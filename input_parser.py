# input_parser.py
# By Sebastian Raaphorst, 2020.

# Parses the input from the astropy tables.

from time_units import *
import numpy as np
from astropy.time import Time as ATime
from astropy.table import Table
from ilp_solver import *
from time import monotonic


def read_tables(time_table: str,
                observation_table: str,
                target_table_metric_visibility_hour_angle: str) -> (TimeSlots, List[Observation]):
    print('Reading tables...')
    start_time = monotonic()

    timetab = Table.read(time_table)
    obstab = Table.read(observation_table)
    targtab_metvisha = Table.read(target_table_metric_visibility_hour_angle)

    # *** PROG / OBS IDS / BANDS ***
    # Get the obs_id of the observations we are considering.
    obs_ids = [row['obs_id'] for row in obstab]

    band_map = {1: Band.Band1, 2: Band.Band2, 3: Band.Band3, 4: Band.Band4}
    bands = {obs_id: band_map[obstab[obs_id == obs_id]['band']] for obs_id in obs_ids}

    # *** TIMING ***
    time_strings = timetab['time']
    time_slot_length = Time(int((ATime(time_strings[1]) - ATime(time_strings[0])).to('minute').round(2).to_value()),
                            unit=TimeUnits.minutes)

    # TODO: This is hard-coded for now. We have to find a way to determine this information from Bryan
    # TODO: using 12 degree twilights. We can do it now with targtab_metvisha's weight_g[sn], but not easily.
    # GS nighttime runs from the specified time_string to the specified time_string, and same for GN.
    # We will make these distinct with no overlap, which is necessary for the ILP algorithm, and will function for the
    # GA algorithm.
    gs_start_time_slot, gs_end_time_slot = 0, 516
    gs_time_slots = gs_end_time_slot - gs_start_time_slot + 1
    gn_start_time_slot, gn_end_time_slot = 264, 927
    gn_time_slots = gn_end_time_slot - gn_start_time_slot + 1

    # Time slots shared in file between GS and GN, starting at gn_start_time_slot and ending at gs_end_time_slot.
    time_slot_overlap = gs_end_time_slot - gn_start_time_slot

    # *** METRIC INFO ***
    # TODO: Should we filter out the 0s?
    # weight_gs should always only be defined between [0:gs_end_time_slot+1].
    # weight_gn should always only be defined between [gn_start_time_slot:-1].

    # Scalars
    # NOTE: weight_gs = metric * visfrac * wha, so we only need weight_gs and weight_gn, along with the range of
    # possible starting positions.
    # metric = {obs_id: obstab[obs_id == obs_id]['metric'] for obs_id in obs_ids}
    # visfrac_gs = {obs_id: obstab[obs_id == obs_id]['visfrac_gs'] for obs_id in obs_ids}
    # visfrac_gn = {obs_id: obstab[obs_id == obs_id]['visfrac_gn'] for obs_id in obs_ids}

    # *** WEIGHTS AND SITE ***
    # TODO: why is wha_gs / wha_gn longer than the number of time_slots for the night?
    weight_gs = {obs_id: targtab_metvisha[idx]['weight_gs'][:gs_time_slots] for idx, obs_id in enumerate(obs_ids)}
    weight_gn = {obs_id: targtab_metvisha[idx]['weight_gn'][gn_start_time_slot:] for idx, obs_id in enumerate(obs_ids)}
    weights = {obs_id: np.append(weight_gs[obs_id], weight_gn[obs_id]) for obs_id in obs_ids}
    # wha_gs = {obs_id: targtab_metvisha[obs_id == obs_id]['wha_gs'] for obs_id in obs_ids}
    # wha_gn = {obs_id: targtab_metvisha[obs_id == obs_id]['wha_gn'] for obs_id in obs_ids}

    # *** OBS LENGTHS ***
    # Get the remaining observation lengths.
    obs_lengths = {row['obs_id']: Time(row['tot_time'] - row['obs_time'], TimeUnits.hours) for row in obstab}

    # The start slots are the dict of lists of indices of the start slots, indexed by obs id.
    # Example: {'GS-2018B-FT-208-18': [15, 16, ... ,86], ...}
    start_slots = {}
    start_slots_gs = {}
    start_slots_gn = {}
    scheduling_unit_length = Time(3)
    units = {obs_id: SchedulingUnits(ceil(obs_lengths[obs_id].mins() / scheduling_unit_length.mins()), 
                                     scheduling_unit_length) for obs_id in obs_ids}
   
    for  obs_id in obs_ids:
        # We don't have observations yet, so this calculation must be repeated.
        #time_slots_needed = ceil(obs_lengths[obs_id].mins() / time_slot_length.mins())
        print(obs_id)
        print(units[obs_id].duration())
        time_slots_needed = ceil(units[obs_id].duration()/time_slot_length.mins())
        # GS time slots
        start_slots_gs[obs_id] = [time_slot_idx for time_slot_idx, weight in enumerate(weight_gs[obs_id]) if weight > 0]

        # GS time slots
        start_slots_gn[obs_id] = [time_slot_idx + gs_time_slots for time_slot_idx, weight in enumerate(weight_gn[obs_id]) if weight > 0]

        # Now drop the slots we cannot start in because we don't have enough time to do so.
        if time_slots_needed > 1:
            start_slots_gs[obs_id] = start_slots_gs[obs_id][:-(time_slots_needed - 1)]
            start_slots_gn[obs_id] = start_slots_gn[obs_id][:-(time_slots_needed - 1)]
        start_slots[obs_id] = start_slots_gs[obs_id] + start_slots_gn[obs_id]
        #units[obs_id] = SchedulingUnits(units_needed, scheduling_unit_length)

    # Create the time slots: in the example data, there should be 173, each of 3 minutes (the granularity).
    time_slots = TimeSlots(time_slot_length, gs_time_slots, gn_time_slots, time_slot_overlap)

    south = {obs_id: len(start_slots_gs[obs_id]) > 0 for obs_id in obs_ids}
    north = {obs_id: len(start_slots_gn[obs_id]) > 0 for obs_id in obs_ids}
    sites = {}
    for obs_id in obs_ids:
        if south[obs_id] and north[obs_id]:
            sites[obs_id] = Site.Both
        elif south[obs_id]:
            sites[obs_id] = Site.GS
        elif north[obs_id]:
            sites[obs_id] = Site.GN
        else:
            sites[obs_id] = None

    # TODO: Hacks to make greedy-max's partial scheduling scheme work with CFHT's genetic algorithm.
    # obs_lengths['GS-2018B-Q-224-34'] is okay: 6 timeslots.
    ## obs_lengths['GS-2018B-Q-207-48'] = Time((43 - 7 + 1) * time_slots.time_slot_length.mins(), TimeUnits.minutes)
    ## obs_lengths['GS-2018B-Q-218-342'] = Time((86 - 44 + 1) * time_slots.time_slot_length.mins(), TimeUnits.minutes)
    ## obs_lengths['GS-2018B-Q-218-363'] = Time((123 - 87 + 1) * time_slots.time_slot_length.mins(), TimeUnits.minutes)
    # obs_lengths['GS-2019A-Q-229-18'] is okay: 5 timeslots.
    # obs_lengths['GS-2018B-Q-112-24'] is okay: some time already used, 13 remaining timeslots.
    # obs_lengths['GS-2018B-Q-112-25'] is okay: some time already used, 15 timeslots.
    # obs_lengths['GS-2018B-Q-112-26'] is okay: some time already used, 16 timeslots.


    distab = [row.lower() for row in obstab['disperser']]
    dispersers = {obs_id: distab[i] for i, obs_id in enumerate(obs_ids) }

    # Create the observations.
    # We parse out the observations that have no start_slot_priorities, i.e. they cannot be started.
    observations = []
    for obs_id in obs_ids:
        # Skip over unschedulable observations: not enough time to schedule them.
        if max(weights[obs_id]) == 0 or (len(start_slots_gs[obs_id]) == 0 and len(start_slots_gn[obs_id]) == 0):
            continue
        observations.append(Observation(obs_id, sites[obs_id], bands[obs_id], obs_lengths[obs_id],
                                        start_slots[obs_id], weights[obs_id], units[obs_id], dispersers[obs_id]))

    print(f"Done reading tables: {monotonic() - start_time} s")
    return time_slots, observations

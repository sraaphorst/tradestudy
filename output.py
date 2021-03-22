# output

from common import *
from typing import Union
from termcolor import colored

def convert_to_schedule(time_slots: TimeSlots,
                        observations: List[Observation],
                        scheduler: Union[None, Scheduling]) -> Union[None, Schedule]:
    """
    Convert scheduling to a schedule.
    """
    schedule = [[],[]]
    for site in {Site.GS, Site.GS}:
        for time_slot_idx, obs_idx in scheduler[site]:
            if len(schedule[site]) > time_slot_idx:
                raise ValueError(f'Observation {obs_idx} illegally scheduled at time slot {time_slot_idx}')
            if len(schedule[site]) < time_slot_idx:
                schedule[site] += [None] * (time_slot_idx - len(schedule[site]))
            obs = observations[obs_idx]
            schedule[site] += [obs_idx] * obs.time_slots_needed(time_slots)
        if len(schedule[site]) < time_slots.num_time_slots_per_site:
            schedule[site] += [None] * (time_slots.num_time_slots_per_site[site] - len(schedule[site]))
        return schedule


def convert_to_scheduling(schedule: Union[None, Schedule]) -> Union[None, Scheduling]:
    """
    Convert a schedule to the easier-to-read scheduling.
    Here is an example schedule:
    [None, None, 4, 4, 4, 2, 2, None, 1, 1, 1, None]
    which means that at time slots:
    0-1: nothing is scheduled.
    2-4: observation 4 is scheduled.
    5-6: observation 2 is scheduled.
    7: nothing is scheduled.
    8-10: observation 1 is scheduled.
    11: nothing is scheduled.

    Convert to:
    [(2, 4), (5, 2), (8, 1)]

    :param schedule: the schedule as solved by the ILP solver or genetic algorithm
    :return: the scheduling equivalent for easier use, or None if there is nothing scheduled
    """
    if schedule is None:
        return None
    scheduling = [[],[]]
    prev_obs_idx = None
    for site in {Site.GS,Site.GN}:
        for time_slot_idx, obs_idx in enumerate(schedule[site]):
            if obs_idx is None:
                continue
            if prev_obs_idx is None or prev_obs_idx != obs_idx:
                prev_obs_idx = obs_idx
                scheduling[site].append((time_slot_idx, obs_idx))

    if len(scheduling) == 0:
        return None

    return scheduling

def calculate_observation_score(site: Site,
                                time_slots: TimeSlots,
                                obs: Observation,
                                initial_slot_idx: int) -> Union[None, float]:
    """
    Given an observation at a site scheduled for a time slot, calculate its score.
    """
    if obs is None:
        return None
    time_slots_needed = obs.time_slots_needed(time_slots)

    score = 0
    overall_time_slot_index = time_slots.get_time_slot(site, initial_slot_idx).idx
    time_slot_length = time_slots.time_slot_length.mins()
    for time_slot_idx in range(time_slots_needed):
        score += obs.weights[time_slot_idx + overall_time_slot_index] * time_slot_length
    return score

def length_of_night_in_mins(site: Site, time_slots: TimeSlots) -> float:
    """
    Returns the length of the night for the given site in minutes.
    """
    return time_slots.num_time_slots_per_site[site] * time_slots.time_slot_length.mins()


def calculate_schedule_score(time_slots: TimeSlots,
                             observations: List[Observation],
                             schedule: Union[None, Schedule]) -> Union[None, float]:
    """
    Given a schedule, e.g [2, 2, 2, 2, None, 3, 3, ...], calculate its score by multiplying the priority
    by the rough integration of the weight.
    """
    if schedule is None:
        return None
    time_slot_length = time_slots.time_slot_length.mins()
    return sum((observations[obs_idx].weights[time_slot_idx] * time_slot_length
                for time_slot_idx, obs_idx in enumerate(schedule)
                if obs_idx is not None)) / length_of_night_in_mins(time_slots)

def calculate_scheduling_score(site: Site, 
                               time_slots: TimeSlots,
                               observations: List[Observation],
                               scheduling: Union[None, Scheduling]) -> Union[None, float]:
    if scheduling is None:
        return None
    
    if site == Site.Both:
        return calculate_scheduling_score(Site.GS, time_slots, observations , scheduling) + \
                calculate_scheduling_score(Site.GN, time_slots, observations, scheduling)
    return sum([calculate_observation_score(site, time_slots, observations[obs_idx], start_slot)
                for start_slot, obs_idx in scheduling[site]]) / length_of_night_in_mins(site, time_slots)


def print_schedule(site: Site, time_slots: TimeSlots, observations: List[Observation], schedule: Schedule) -> None:
    """
    Output the schedule in an easy-to-read format. This should consist of something akin to:

    Site: Gemini South
    Fitness:

    :param time_slots: the time slots
    :param observations: the list of all observations: the observation of obs_idx is at position obs_idx
    :param site: the site, either Site.GS or Site.GN
    :param schedule: the schedule, a list of length time_slots.number_of_time_slots_per_site
    """
    if schedule is None:
        return

    scheduling = convert_to_scheduling(schedule)
    score = calculate_scheduling_score(site, time_slots, observations, scheduling)
    length_of_night = length_of_night_in_mins(site, time_slots)

    # TODO: Now scheduling will have entries of the form (time_slot_index, observation_index)
    # TODO: Here you should iterate over scheduling and output something along the lines of an aligned table like:
    # Site: Gemini South
    # Observation    Band    Start Time    Length    Priority    Hour Angle Priority    Total Priority
    # ...one entry here for each element of scheduling...
    # where:
    # Observation is observation.name
    # Band is observation.band (will be a Band, which is an enum that is a number between 1 and 4)
    # Start Time is the start time of the observation, which in minutes is time_slots.get_time_slot(site, time_slot_index).mins()
    # Length is observation.obs_time.mins(), preferably rounded to two decimal places
    # You should also calculate the end time of the observation, which is the start time + length.
    # Priority is observation.priority, preferably rounded to five decimal places
    # Hour Angle Priority is observation.start_slot_map[time_slot_index], preferably rounded to five decimal places
    # Total Priority is observation.priority * observation.start_slot_map[time_slot_index] * observation.obs_time.mins() / total_time,
    #    preferably rounded to five decimal places
    #
    # The reason to calculate the length is to see if there are any gaps in the schedule and report them.
    # When you are calculating a row, if you see that:
    #    start time of this observation - end time of previous observation > time_slots.time_slot_length, you should
    # put in a row that just says:
    # Gap of x minutes
    # where x = start time of this observation - end time of previous observation.
    pass


def _detailed_scheduling(name: Union[None, str],
                         site:Site,
                         scheduling: Union[None, Scheduling],
                         time_slots: TimeSlots,
                         observations: List[Observation]) -> Union[None, str]:
    """
    Create a string representing the scheduling in detail.
    """
    if scheduling is None:
        return None
    # TODO: Remove this print done :)
    #print(scheduling)
    time_slot_length = time_slots.time_slot_length.mins()

    # We want to indent if there is a name specified.
    line_start = '\n\t' if name is not None else '\n'

    # Header
    data = name if name is not None else ''

    # Now iterate over the scheduling, reporting the observations scheduled and any gaps in the schedule.
    obs_prev_time = 0
    for obs_start_time_slot, obs_idx in scheduling:
        obs = observations[obs_idx]
        obs_start_time = obs_start_time_slot * time_slots.time_slot_length.mins()

        # If we have a gap that is at least as big as a time slot length, report it; otherwise, we don't, as it is
        # considered insignificant since it doesn't represent an unscheduled time slot.
        gap_size = int(obs_start_time - obs_prev_time)
        if gap_size > time_slot_length:
            data += line_start + f'Gap of  {gap_size:>3} min{"s" if gap_size > 1 else ""}'
        data += line_start + f'At time {obs_start_time:>3}: Observation {observations[obs_idx].name:<15}, ' \
                             f'site={Site(observations[obs_idx].site).name:<4}, ' \
                             f'obs_time={round(observations[obs_idx].obs_time.mins(), 2):>5}, ' \
                             f'score={calculate_observation_score(site, time_slots, obs, obs_start_time_slot)}'
        obs_prev_time = obs_start_time + observations[obs_idx].obs_time.mins()

    # Is there a final gap at the end of the schedule?
    gap_size = int(length_of_night_in_mins(site, time_slots) - obs_prev_time)
    if gap_size > 0e-3:
        data += line_start + f'Gap of  {gap_size:>3} min{"s" if gap_size > 1 else ""}'

    return data


def print_schedule2(time_slots: TimeSlots, observations: List[Observation],
                    gn_schedule: Schedule, gs_schedule: Schedule) -> None:
    """
    Given schedules, print the results.
    """

    gn_sched = convert_to_scheduling(gn_schedule)
    gs_sched = convert_to_scheduling(gs_schedule)



    # *** GN ***
    printable_schedule_gn = _detailed_scheduling("Gemini North:", Site.GN, gn_sched, time_slots, observations)
    print(printable_schedule_gn)

    gn_obs = set([obs_idx for obs_idx in gn_schedule if obs_idx is not None])
    gn_usage = sum(observations[obs_idx].obs_time.mins() for obs_idx in gn_obs)
    gn_pct = gn_usage / length_of_night_in_mins(Site.GN, time_slots) * 100
    gn_score = calculate_scheduling_score(Site.GN, time_slots, observations, gn_sched)
    gn_summary = f'\tUsage: {gn_usage}, {gn_pct}%, Score: {gn_score}'
    print(gn_summary + '\n')

    # *** GS ***
    printable_schedule_gs = _detailed_scheduling("Gemini South:", Site.GS, gs_sched, time_slots, observations)
    print(printable_schedule_gs)

    gs_obs = set([obs_idx for obs_idx in gs_schedule if obs_idx is not None])
    gs_usage = sum(observations[obs_idx].obs_time.mins() for obs_idx in gs_obs)
    gs_pct = gs_usage / length_of_night_in_mins(Site.GS, time_slots) * 100
    gs_score = calculate_scheduling_score(Site.GS, time_slots, observations, gs_sched)
    gs_summary = f'\tUsage: {gs_usage} mins, {gs_pct}%, Score: {gs_score}'
    print(gs_summary)

def print_schedule3(time_slots: TimeSlots, observations: List[Observation], schedule: Schedule) -> None:

    sched = convert_to_scheduling(schedule)
    gn_sched = sched[Site.GN]
    gs_sched = sched[Site.GS]
    gn_schedule, gs_schedule  = schedule[0],schedule[1]

    # *** GN ***
    printable_schedule_gn = _detailed_scheduling("Gemini North:", Site.GN, gn_sched, time_slots, observations)
    print(printable_schedule_gn)

    gn_obs = set([obs_idx for obs_idx in gn_schedule if obs_idx is not None])
    gn_usage = sum(observations[obs_idx].obs_time.mins() for obs_idx in gn_obs)
    gn_pct = gn_usage / length_of_night_in_mins(Site.GN,time_slots) * 100
    gn_score = calculate_scheduling_score(Site.GN, time_slots, observations, sched)
    gn_summary = f'\tUsage: {gn_usage}, {gn_pct}%, Score: {gn_score}'
    print(gn_summary + '\n')

    # *** GS ***
    printable_schedule_gs = _detailed_scheduling("Gemini South:", Site.GS, gs_sched, time_slots, observations)
    print(printable_schedule_gs)

    gs_obs = set([obs_idx for obs_idx in gs_schedule if obs_idx is not None])
    gs_usage = sum(observations[obs_idx].obs_time.mins() for obs_idx in gs_obs)
    gs_pct = gs_usage / length_of_night_in_mins(Site.GS,time_slots) * 100
    gs_score = calculate_scheduling_score(Site.GS, time_slots, observations, sched)
    gs_summary = f'\tUsage: {gs_usage} mins, {gs_pct}%, Score: {gs_score}'
    print(gs_summary)
    #obs =  set([obs_idx for obs_idx in schedule if obs_idx is not None])

def in_line_print(schedule):
    sched = []
    print('GS')
    for i in schedule[Site.GS]:
        if i:
            sched.append(colored(i,'red'))
        else:
            sched.append('Free')
    
    print(' '.join(sched))
    sched = []
    print('GN')
    for i in schedule[Site.GN]:
        if i:
            sched.append(colored(i,'red'))
        else:
            sched.append('Free')
    print(' '.join(sched))
    
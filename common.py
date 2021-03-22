# common.py
# By Sebastian Raaphorst, 2020.

# Contains all the code common to both the ILP solver and the genetic algorithm solver.
# We use time slots in both cases. In both cases, time slots can be given a granularity of
# a certain number of minutes, which is the minimal schedule time chunk.
# In the original implementation, for ILPs, the granularity was three minutes.
# For genetic algorithms, it was one minute.

from enum import IntEnum
from typing import List, Union, Dict, Tuple
from math import ceil
from time_units import Time

# The schedule consists of a list corresponding between timeslot indices and observation indices, or None
# if nothing is scheduled at that timeslot index.
# For example, if we have:
# >>> schedule = List(12, 12, 12, None, None, 1, 1, 1, 1, None)
# We would have:
# >>> len(schedule) timeslots, i.e. 10.
# >>> schedule[0:3] = 12, i.e. observation 12 is scheduled for timeslots 0, 1, 2.
# >>> schedule[3:5] = None, i.e. no observations are scheduled for timeslots 3, 4.
# >>> schedule[5:9] = 1, i.e. observation 1 is scheduled for timeslots 5, 6, 7, 8.
# >>> schedule[9:10] = None, i.e. no observations are scheduled for timeslots 10.
Schedule = List[Union[int, None]]

Scheduling = List[Tuple[int, int]]

# The final score of the schedule.
Score = float

# Metrics.
Metric = float

# For observations, a map from available time slot indices to metric
Weights = Dict[int, Metric]


# Define the sites.
class Site(IntEnum):
    """
    The sites (telescopes) available to an observation.
    """
    GN = 0
    GS = 1
    Both = 2


# The band for observations.
class Band(IntEnum):
    """
    The band to which an observation is scheduled.
    """
    Band1 = 1
    Band2 = 2
    Band3 = 3
    Band4 = 4

class TimeSlot:
    counter = 0

    """
    A single representation of a time slot, which is a site, a start time, and a length of time.
    Time can be of any unit, but we typically work with minutes for convenience.
    """
    def __init__(self, site: Site, start_time: Time, length: Time):
        """
        The definition of a time slot available for a site.
        Note that time slots should NOT overlap, i.e. for a site, e.g. GN, if we have time slots:
        1. TimeSlot(GN, Time(0), Time(10)), a TimeSlot(GN, 9, 10) would be illegal as the times intersect at 9.
        The timeslot begins at start_time and lasts for TIMESLOT_LENGTH.

        :param site: the site in question
        :param start_time: the start time of the time slot as an offset from 0, the initial start time
        :param length: the length of the time slot
        """
        if site == Site.Both:
            raise ValueError("TimeSlots must be created for a specific site.")

        self.site = site
        self.start_time = start_time
        self.length = length
        self.idx = TimeSlot.counter
        TimeSlot.counter += 1

    def end_time(self) -> Time:
        """
        Returns the ending time of the time slot for convenience.
        :return: a Time object representing the end of the time slot
        """
        return self.start_time + self.length

    def __hash__(self):
        return hash((self.site, self.start_time))


class TimeSlots:
    """
    A collection of TimeSlot objects.
    """
    # TODO: Change the default value here
    def __init__(self,
                 time_slot_length: Time,
                 gs_time_slots: int,
                 gn_time_slots: int,
                 time_slot_overlap: int):
        """
        Create the collection of time slots, which consist of a collection of TimeSlot
        objects as below for scheduling.

        For each resource, we begin at time 0 and create a time slot.
        We then increment by timeslot_length and create another time slot.
        We continue to create time slots until we have the specified number for each resource.

        :param time_slot_length: the length of the time slots (i.e. the granularity of the schedule)
        :param gs_time_slots: the number of time slots at GS
        :param gn_time_slots: the number of time slots at GN
        :param time_slot_overlap: the overlap in the files between the time slots:
                                  gn_start_time_slots to gs_end_time_slots.
        """
        self.time_slot_length = time_slot_length
        self.num_time_slots_per_site = {Site.GS: gs_time_slots, Site.GN: gn_time_slots}
        self.total_time_slots = gs_time_slots + gn_time_slots
        self._time_slot_overlap = time_slot_overlap

        self._time_slots = []
        for site in [Site.GS, Site.GN]:
            for idx in range(self.num_time_slots_per_site[site]):
                # TODO: Is this calculation correct?
                # TODO: I think so. GN time slots are advanced by the # of GS time slots and then moved back
                # TODO: by the overlap to put them at the correct time.
                time_slot_start = idx + (0 if site == Site.GS else
                                         self.num_time_slots_per_site[Site.GS] - self._time_slot_overlap)
                self._time_slots.append(TimeSlot(site,
                                                 Time(time_slot_start * time_slot_length.mins()),
                                                 Time(time_slot_length.mins())))

    def get_time_slot(self, site: Site, index: int) -> Union[None, TimeSlot]:
        """
        Given a site and an index into its time slots, return the corresponding time slot.
        :param site: the Site
        :param index: the index, in [0, _num_time_slots_per_site[site])
        :return: the TimeSlot, if it exists
        :except: ValueError if the index condition is violated
        """
        if site == Site.Both:
            raise ValueError("get_time_slot requires a specific site")
        elif site == Site.GS:
            return self._time_slots[index] if index < self.num_time_slots_per_site[Site.GS] else None
        else:
            return self._time_slots[self.num_time_slots_per_site[Site.GS] + index] \
                if index < self.num_time_slots_per_site[Site.GN] else None

    def __iter__(self):
        """
        Create an iterator for the time slots.
        :return: an iterator
        """
        return _TimeSlotsIterator(self)


# TODO: Do we need this at all? How do we modify it? To return all the GS slots, and then the GN slots?
class _TimeSlotsIterator:
    """
    Iterator class for TimeSlots.
    """
    def __init__(self, time_slots: TimeSlots):
        self._timeslots = time_slots
        self._index = 0

    def __next__(self) -> TimeSlot:
        """
        Returns the next time slot.
        :return: the next time slot
        :except: StopIteration when the iteration is done
        """
        if self._index < len(self._timeslots._time_slots):
            slot = self._timeslots._time_slots[self._index]
            self._index += 1
            return slot
        raise StopIteration

class SchedulingUnits:
    """
    Basic unit of a schedule. A group of this unit 
    """
    def __init__(self, units: int, length: Time):
        self.units = units
        self.length = length
        self.in_use = self.units 

    def duration(self) -> Time:
        return self.in_use *self.length.mins()
    
    def reduce(self, amount: int) -> bool:
        if self.in_use < amount:
            return False
        
        self.in_use = amount
        return True

    def expand(self) -> None:
        self.in_use = self.units

class Entirety(IntEnum):

    Complete = 0
    Partial = 1


class Observation:
    """
    The basic information that comprises an observation.
    """
    # Keep a static counter of all observations.
    _counter = 0

    def __init__(self, 
                 name: str,
                 site: Site,
                 band: Band,
                 obs_time: Time,
                 start_slots: List[int],
                 weights: Weights,
                 units: SchedulingUnits,
                 disperser: str,
                 allocated_time: Time = None,):
        self.name = name
        self.idx = Observation._counter
        self.site = site
        self.band = band
        self.used_time = Time(0)
        # new atom model
        self.units = units 
        self.disperser = disperser
        self.acq_overhead = Time(0.2) if 'mirror' == self.disperser else Time(0.3)
        self.obs_time = Time(self.acq_overhead.mins() + units.duration())
        self.entirety = Entirety.Complete
        # END new 
        self.allocated_time = obs_time if allocated_time is None else allocated_time
        self.start_slots = start_slots
        self.weights = weights

        Observation._counter += 1

    def time_slots_needed(self, time_slots: TimeSlots) -> int:
        """
        Determine the number of time slots needed (we use ceiling to get this number).
        :param self: the observation
        :return: the number of required time slots
        """
        #return ceil(self.obs_time.mins() / time_slots.time_slot_length.mins())
        return ceil(self.units.duration() / time_slots.time_slot_length.mins())


# __init__.py
# By Sebastian Raaphorst, 2020.

from common import *
from random import choice, randrange, sample, shuffle
from copy import copy
from typing import List, Tuple, Union
from output import calculate_schedule_score, calculate_scheduling_score, convert_to_schedule
from time_units import Time
from itertools import groupby
from termcolor import colored

# A Chromosome maintains a Schedule as defined in common, i.e. a list where the contents of position x are the obs_id
# of the observation scheduled for time_slot x, and None if nothing is scheduled.
# A chromosome is linked to a site. We must handle this carefully for the following reasons:
# 1. For observations that can be scheduled at both, we want them in a chromosome for each site.
# 2. We do not want any observation to be scheduled twice.



class Chromosome:
    def __init__(self, time_slots: TimeSlots, observations: List[Observation]):
        self.time_slots = time_slots
        self.observations = observations
        self.schedule = [[None] * time_slots.num_time_slots_per_site[Site.GN],
                            [None] * time_slots.num_time_slots_per_site[Site.GS]]
        # Pairs of the form (time_slot_idx, obs_idx) where an observation has been scheduled.
        self.scheduling = [[],[]]
        self.on_site = None
        self.partitions = 0

    def __copy__(self):
        c = Chromosome(self.time_slots, self.observations)
        c.schedule = self.schedule.copy()
        c.scheduling = self.scheduling.copy()
        return c

    def on_site(self, site: Site):
        """
        Sets on_site variable to perform different operations to only that site part of of the chromosome 
        """
        if site not in {Site.Both, Site.GS, Site.GN}:
           raise ValueError("Site does not exit.")
        self.on_site = site
    
    def _get_first_gap(self, obs_idx: int, site: int) -> Union[int, None]:
        """
        Given an observation index, if it can be scheduled in this chromosome (the sites are compatible), determine
        the first gap in which it can be scheduled.
        """
        obs = self.observations[obs_idx]
        #print(obs_idx,obs.site,site)
        if obs.site not in {Site.Both, Site.GS, Site.GN}:
            return None

        # We can only schedule between the lower time and the upper time permissible for the chromosome.
        # Get the indices of the minimum and maximum start slots.
        self._min_obs_slot_idx, self._max_obs_slot_idx = None, None
        
        if site == Site.GS:
            #print(obs.start_slots,self.time_slots.num_time_slots_per_site[Site.GS],site, obs.site)
            self._min_obs_slot_idx = min(
                [i for i in obs.start_slots if i < self.time_slots.num_time_slots_per_site[Site.GS]])
            self._max_obs_slot_idx = max(
                [i for i in obs.start_slots if i < self.time_slots.num_time_slots_per_site[Site.GS]])

        elif site == Site.GN:
            #print(obs.start_slots,self.time_slots.num_time_slots_per_site[Site.GS],site, obs.site)
            self._min_obs_slot_idx = min(
                [i for i in obs.start_slots if self.time_slots.num_time_slots_per_site[Site.GS] <= i])
            self._max_obs_slot_idx = max(
                [i for i in obs.start_slots if self.time_slots.num_time_slots_per_site[Site.GS] <= i])
       
        #print(f'min:{self._min_obs_slot_idx} max:{self._max_obs_slot_idx}')
        if obs.entirety == Entirety.Partial:
            # Check if the observation when expanded can be schedule
            og_units = obs.units.in_use
            obs.units.expand()
            slots_needed = obs.time_slots_needed(self.time_slots)
            offset = self.time_slots.num_time_slots_per_site[Site.GS] if site == Site.GN else 0
            unused_time_slots = [time_slot_idx for time_slot_idx, obs_idx in enumerate(self.schedule[site]) if obs_idx is None
                             and time_slot_idx + offset in range(self._min_obs_slot_idx, self._max_obs_slot_idx + 1)]
            for time_slot_idx in unused_time_slots:
                # Check to see if we have slots_needed empty slots starting at time_slot_idx.
                slots_to_check = self.schedule[site][time_slot_idx:(time_slot_idx + slots_needed)]
                if set(slots_to_check) == {None} and len(slots_to_check) == slots_needed:
                    return time_slot_idx
            
            #return to original
            obs.units.reduce(og_units) 
            
        slots_needed = obs.time_slots_needed(self.time_slots)
        

        # Get the sorted indices of the unused time_slots that we can use for scheduling this observation.
        # TODO: WE NEED TO OFFSET TIME_SLOT_IDX BY THE OFFSET, IN THIS CASE
        offset = self.time_slots.num_time_slots_per_site[Site.GS] if site == Site.GN else 0
        
        unused_time_slots = [time_slot_idx for time_slot_idx, obs_idx in enumerate(self.schedule[site]) if obs_idx is None
                             and time_slot_idx + offset in range(self._min_obs_slot_idx, self._max_obs_slot_idx + 1)]
        #print('*** unused slots process ***')
  
        #print([(time_slot_idx,time_slot_idx+offset) for time_slot_idx, obs_idx in enumerate(self.schedule[site]) if obs_idx is None])
        #print('*** END slots process ***')
        #print(f'unused_time_slots: {unused_time_slots}')
        # Now iterate over the unused time slots and determine if this observation can be inserted in the position.
        for time_slot_idx in unused_time_slots:
            # Check to see if we have slots_needed empty slots starting at time_slot_idx.
            slots_to_check = self.schedule[site][time_slot_idx:(time_slot_idx + slots_needed)]

            if set(slots_to_check) == {None} and len(slots_to_check) == slots_needed:
                return time_slot_idx

        # Determine the number of time_slots we need to accommodate this observation.
        
        slots_needed = obs.time_slots_needed(self.time_slots)
        

        # Get the sorted indices of the unused time_slots that we can use for scheduling this observation.
        # TODO: WE NEED TO OFFSET TIME_SLOT_IDX BY THE OFFSET, IN THIS CASE
        offset = self.time_slots.num_time_slots_per_site[Site.GS] if site == Site.GN else 0
        
        unused_time_slots = [time_slot_idx for time_slot_idx, obs_idx in enumerate(self.schedule[site]) if obs_idx is None
                             and time_slot_idx + offset in range(self._min_obs_slot_idx, self._max_obs_slot_idx + 1)]
        #print('*** unused slots process ***')
  
        #print([(time_slot_idx,time_slot_idx+offset) for time_slot_idx, obs_idx in enumerate(self.schedule[site]) if obs_idx is None])
        #print('*** END slots process ***')
        #print(f'unused_time_slots: {unused_time_slots}')
        # Now iterate over the unused time slots and determine if this observation can be inserted in the position.
        for time_slot_idx in unused_time_slots:
            # Check to see if we have slots_needed empty slots starting at time_slot_idx.
            slots_to_check = self.schedule[site][time_slot_idx:(time_slot_idx + slots_needed)]

            if set(slots_to_check) == {None} and len(slots_to_check) == slots_needed:
                return time_slot_idx

        return None

    def _split_observation(self, obs_idx: int, site: int) -> Union[int,None]:
        
        obs = self.observations[obs_idx]
        if obs.obs_time.mins() < 30: #check that for minutes; lower limit 
            return None, None

        slots_needed = obs.time_slots_needed(self.time_slots)
        offset = self.time_slots.num_time_slots_per_site[Site.GS] if site == Site.GN else 0
        time_slot_length = self.time_slots.time_slot_length.mins()
        

        self._min_obs_slot_idx, self._max_obs_slot_idx = None, None
        
        if site == Site.GS:
            self._min_obs_slot_idx = min(
                [i for i in obs.start_slots if i < self.time_slots.num_time_slots_per_site[Site.GS]])
            self._max_obs_slot_idx = max(
                [i for i in obs.start_slots if i < self.time_slots.num_time_slots_per_site[Site.GS]])

        elif site == Site.GN:
            self._min_obs_slot_idx = min(
                [i for i in obs.start_slots if self.time_slots.num_time_slots_per_site[Site.GS] <= i])
            self._max_obs_slot_idx = max(
                [i for i in obs.start_slots if self.time_slots.num_time_slots_per_site[Site.GS] <= i])
        unused_time_slots = [time_slot_idx for time_slot_idx, obs_id in enumerate(self.schedule[site]) if obs_id is None
                             and time_slot_idx + offset in range(self._min_obs_slot_idx, self._max_obs_slot_idx + 1)]       

        if  len(unused_time_slots) < 1:
              return None, None

        unused_time_windows = [] #this might be done with comprehension. 
        for i in range(len(unused_time_slots)):
            if unused_time_slots[i-1]+1 == unused_time_slots[i]:
                unused_time_windows[-1].append(unused_time_slots[i]) 
            else:
                 unused_time_windows.append([unused_time_slots[i]])

        for tw in unused_time_windows:
            #if len(tw)*time_slot_length > 30: #check
            if obs.units.duration() > 30: # check if observation is longer than 30min 
                slots_to_check = self.schedule[site][tw[0]:tw[-1]]
                # the remaining observation also needs to be longer than 30 min
                if set(slots_to_check) == {None} and (slots_needed - len(tw)) > 30:
                    if obs.units.reduce(len(tw)):
                        obs.entirety = Entirety.Partial
                        self.partitions += 1 
                    #remaining_slots = slots_needed - len(tw)
                    #partial_obs = copy(obs)
                    #partial_obs.obs_time = Time(remaining_slots * time_slot_length)
                    #self.observations[obs_idx] = partial_obs
                    #self.partitions += 1
                        return tw[0] , len(tw)
       
        return None, None
        
    def determine_capacity(self, site: int) -> float:
        """
        Determine the amount of time currently used up in this chromosome.
        >= 85% is considered "optimal."
        :return: True if optimal, and False otherwise.
        """
        #return len([i for i in self.schedule if i is not None]) / len(self.schedule) >= 0.85
        return len([i for i in self.schedule[site] if i is not None])

    def determine_fitness(self) -> float:
        """
        Determine the value of the Chromosome, which is just the sum of the metric over its observations multiplied by
        the score for the times chosen.
        :return: the sum of the metric over the observations
        """
        return calculate_scheduling_score(Site.Both, self.time_slots, self.observations, self.scheduling)

    def insert(self, obs_idx, site) -> bool:
        """
        Try to insert observation obs_idx into this chromosome in the earliest possible position. This fails if:
        1. The observation's site is not compatible with this Chromosome's site.
        2. There are no gaps bit enough to accommodate it.
        3. The observation is already in this chromosome.
        Otherwise, it is scheduled.
        :param obs_idx: the index of the observation to try to schedule
        :return: True if we could schedule, and False otherwise
        """
        if obs_idx is None:
            return False # this validation is new, which is odd so I might be a problem when creating the initial popu
        obs = self.observations[obs_idx]
        if site not in {Site.Both, Site.GS, Site.GN}:
            return False
        
        #Check for the obs in both sites so its not duplicated or if 
        for s in {Site.GS, Site.GN}: 
            if obs_idx in self.schedule[s]: # TODO: Must change this
                return False
        
           
              
        # Get the gaps into which we can insert. 
        start_time_slot_idx = self._get_first_gap(obs_idx, site)

        assert(start_time_slot_idx not in [o[0] for o in self.scheduling[site]]) # Check if timeslot is already on used\
        #print(f'obs site: {obs.site.name}, {site.name}, start_ts: {start_time_slot_idx}') 
            
        if start_time_slot_idx is None:
            
            # Check if the observation can be split in partial obs 
            partial_time_slot_idx, slots_needed = self._split_observation(obs_idx, site)
            if partial_time_slot_idx:

                assert(partial_time_slot_idx not in [o[0] for o in self.scheduling[site]]) # Check if timeslot is already on used

                self.scheduling[site].append((partial_time_slot_idx, obs_idx))
                self.scheduling[site].sort()
                self.schedule[site][partial_time_slot_idx:(partial_time_slot_idx 
                                                            + slots_needed)] = [obs_idx]* slots_needed

            return False
        #print(f'{obs_idx} inserted')
        self.scheduling[site].append((start_time_slot_idx, obs_idx))
        self.scheduling[site].sort()

        # Determine the number of time_slots we need to accommodate this observation and modify the schedule array
        # so that it represents obs_idx being scheduled for slots_needed positions starting at start_time_slot_idx.
        slots_needed = obs.time_slots_needed(self.time_slots)
        self.schedule[site][start_time_slot_idx:(start_time_slot_idx + slots_needed)] = [obs_idx] * slots_needed

        return True

    def remove(self, obs_idx: int, site: Site) -> bool:
        """
        Remove the specified observation from the Chromosome.
        :param obs_idx: the observation to remove
        :return: True if it can be removed, and False otherwise.
        """

        if obs_idx in  self.schedule[site]:
            self.schedule[site] = [i if i != obs_idx else None for i in self.schedule[site]]
            self.scheduling[site] = [(i, j) for (i, j) in self.scheduling[site] if j != obs_idx]
            return True
        return False

    def __len__(self):
        """
        Return the number of observations in this schedule by.
        :return:
        """
        return len(self.scheduling[self.on_site])
        #return  len(self.schedule[0])+len(self.schedule[1])

    def __getitem__(self, idx) -> int:
        """
        Return the idxth overvation in this chromosome.
        """
        return self.scheduling[self.on_site][idx][1]
    
    def __str__(self) -> str:
        return f'fitness: {self.determine_fitness():.6} '+ \
                 f'GS cap: {self.determine_capacity(Site.GS)}/{self.time_slots.num_time_slots_per_site[Site.GS]} '+ \
                 f'GN cap: {self.determine_capacity(Site.GN)}/{self.time_slots.num_time_slots_per_site[Site.GN]} '+ \
                 f'\nGS Sched: {self.scheduling[Site.GS]}\nGN Sched: {self.scheduling[Site.GN]}'

    def __contains__(self, idx) -> bool:
        return True if idx in self.schedule[self.on_site] else False

class GeneticAlgorithm:
    def __init__(self, time_slots: TimeSlots, observations: List[Observation], include_greedy_max=False):
        self.time_slots = time_slots
        self.observations = observations
        self.chromosomes = [] #{Site.GS: [], Site.GN: []}
        self.include_greedy_max = include_greedy_max
        self.swaps = 0
        self.mixes = 0
        self.mates = 0
        self.interleaves = 0 
        self.improve_swap = 0
        self.improve_mix = 0
        self.improve_mate = 0
        self.improve_inter = 0   

    def _form_initial_population(self):
        """
        We form the initial population of chromosomes by putting them at the earliest period that we can.
        Every observation is scheduled in one chromosome per site in which it is allowed.
        """
        # Sort the observations by their maximum potential score.
        # sorted_obs_idx_by_score = [obs.idx for obs in sorted(self.observations,
        #                                                      key=lambda x: x.priority * x.obs_time.mins()
        #                                                                    * np.max(list(x.start_slot_map.values())),
        #                                                      reverse=True)]
        sorted_obs_idx_by_score = [obs.idx for obs in sorted(self.observations,
                                                             key=lambda x: max(x.weights) * x.obs_time.mins())]
        #print(f'slots for GN: {self.time_slots.num_time_slots_per_site[Site.GN]},\
        #        slots for GS: {self.time_slots.num_time_slots_per_site[Site.GS]}')
        for obs_idx in sorted_obs_idx_by_score:
            obs = self.observations[obs_idx]
            #print(f'Inserting {obs.idx} slots needed: {obs.time_slots_needed(self.time_slots)}')
            # Determine the sites in which it should be scheduled.
            sites = {Site.GN, Site.GS} if obs.site == Site.Both else {obs.site}
            gs_sched, gn_sched = 0, 0
            for site in sites:
                scheduled = False   
                for chromosome in self.chromosomes:
                    #print(f'chromosome time left for {site.name}: {chromosome.determine_capacity(site)} chromosome fitness: {chromosome.determine_fitness():.6f}')
                    #print(f'score on gs: {calculate_scheduling_score(Site.GS, self.time_slots, self.observations, chromosome.scheduling):.7f} score on gn: {calculate_scheduling_score(Site.GN, self.time_slots, self.observations, chromosome.scheduling):.7f}')
                    if chromosome.insert(obs_idx,site):
                        #print(f'{obs.idx} inserted in existing chromosome at {site.name}')
                        #print(f'chromosome new fitness: {chromosome.determine_fitness():.6f}')
                        scheduled = True
                        break
                if scheduled and site == Site.GS:
                    gs_sched += 1
                if scheduled and site == Site.GN:
                    gn_sched += 1
           
                # Create a new Chromosome for this site and insert it.   
                if not scheduled:
                    c = Chromosome(self.time_slots, self.observations)
                    #print("create a new chromosome")
                    #print(f'{obs.idx} inserted in new chromosome at {site.name}')
                    #print(c.schedule)
                    #print(obs.site,site)
                    if c.insert(obs_idx,site):
                        self.chromosomes.append(c)
                    else:
                        raise ValueError(f'{obs_idx} could not be scheduled at {Site(obs.site).name}')
            self._sort_chromosomes()
            #input()
        #self._print_population()

        if self.include_greedy_max:
        #     # Add Bryan's chromosome to the GA
            print([o.name for o in self.observations])
            gm_scheduling = [[(9,'GN-2019A-FT-101-9'), (168,'GS-2018B-Q-218-342'), 
                              (277,'GN-2018B-Q-238-116'), (329,'GN-2019A-Q-903-71'), 
                              (491,'GS-2018B-Q-224-34')], 
                            [(281,'GN-2018B-Q-228-31'), (318,'GS-2018B-Q-201-28'),
                             (425,'GS-2018B-Q-131-16'), (462,'GS-2018B-Q-207-62'), 
                             (504,'GN-2019A-Q-229-9')]]
         
            for obs in self.observations:
                sites = {Site.GN, Site.GS} if obs.site == Site.Both else {obs.site}
                for site in sites:
                    for i,t in enumerate(gm_scheduling[site]):
                        if t[1] == obs.name:
                           print(obs.idx)
                           gm_scheduling[site][i]= (t[0],obs.idx)
                           print(gm_scheduling[site][i])
            print(gm_scheduling)

        
            g1_chromosome = Chromosome(self.time_slots, self.observations)
            g1_chromosome.schedule = convert_to_schedule(self.time_slots, self.observations, gm_scheduling)
            g1_chromosome.scheduling = gm_scheduling
            self.chromosomes.append(g1_chromosome)
            print(f'Greedy-max chromosome: fitness: {g1_chromosome.determine_fitness()}, '
                   f'score: {calculate_scheduling_score(self.time_slots, self.observations, b1_scheduling)}')
            self._sort_chromosomes()

    def _sort_chromosomes(self):
        self.chromosomes = sorted(self.chromosomes,
                                            key=lambda x: x.determine_fitness(), reverse=True)

    def _single_selection(self) -> int:
        self._sort_chromosomes()
        return choice([n for n, c in enumerate(self.chromosomes)])

    def _pair_selection(self) -> Tuple[int, int]:
        """
        Select two chromosomes: one from the top 25%, and one completely at random.
        :return: the indices of the chromosomes selected.
        """
        self._sort_chromosomes()
        #c1_index = randrange(ceil(len(c1_candidates) / 4))
        c1_idx = choice(range(len(self.chromosomes)))
        c2_idx = None
        while c2_idx is None or c1_idx == c2_idx:
            c2_idx = choice(range(len(self.chromosomes)))

        if self.chromosomes[c1_idx].determine_fitness() > self.chromosomes[c2_idx].determine_fitness():
            return c1_idx, c2_idx
        else:
            return c2_idx, c1_idx

    def _mate(self):
        """
        Mate two chromosomes. This only works if:
        1. They are from the same site.
        2. The timing of the scheduling does not clash with overlaps (overlaps are just dropped).
        If a valid chromosome is found out of the two candidates, pick the higher fitness one and replace the lower
        fitness one in the chromosome list.

        :return: True if mating succeeded, False otherwise.
        """

        c1_index, c2_index = self._pair_selection()
        c1 = self.chromosomes[c1_index]
        c2 = self.chromosomes[c2_index]

        c3 = Chromosome(self.time_slots, self.observations)
        c4 = Chromosome(self.time_slots, self.observations)
        for site in {Site.GN,Site.GS}:
            
            c1.on_site, c2.on_site, c3.on_site, c4.on_site = site, site, site, site
            
            # Pick a crossover point. We want some of each chromosome, so pick between [1, len-1].
            # If either is too short, we can't mate.
            if len(c1) <= 1 or len(c2) <= 1:
                return False


            # Pick a point from the scheduling from c1 and c2.
            c1_point = randrange(1, len(c1))
            c2_point = randrange(1, len(c2))

            for i in range(c1_point):
                c3.insert(c1[i], site)
            for i in range(c2_point, len(c2)):
                if c2[i] not in set(c3.schedule[site]):
                    c3.insert(c2[i], site)

            for i in range(c2_point):
                c4.insert(c2[i], site)
            for i in range(c1_point, len(c1)):
                if c1[i] not in set(c4.schedule[site]): 
                    c4.insert(c1[i], site)
            
        # Reset site variable
        c1.on_site, c2.on_site, c3.on_site, c4.on_site = None, None, None, None 
        
        #print('MATING PROCESS')
        #print('Chromosomes to interleave:')
        #print(f'c1 {c1}')
        #print(f'c2 {c2}')
        #print('New breed chromosomes :')
        #print(f'c3 {c3}')
        #print(f'c4 {c4}')
        self.mates += 1

        #If we have improvement in one of the matings, then replace the lower-valued chromosome.
        max_c = c3 if c3.determine_fitness() > c4.determine_fitness() else c4
        if max_c.determine_fitness() > c2.determine_fitness() and not self._contains(max_c):
            self.chromosomes[c2_index] = max_c
            self._sort_chromosomes()
            self.improve_mate+=1
            return True

        return False

    def _contains(self, c: Chromosome):
        for c2 in self.chromosomes:
            if c.scheduling == c2.scheduling:
                return True
        return False

    def _interleave(self):
        """
        Perform the interleave operation between chromosomes.
        """
        c1_index, c2_index = self._pair_selection()
        c1 = self.chromosomes[c1_index]
        c2 = self.chromosomes[c2_index]
        
        c3 = Chromosome(self.time_slots, self.observations)
        c4 = Chromosome(self.time_slots, self.observations)
        # Interleave to produce the chromosomes.
        for site in {Site.GN,Site.GS}:
            c1.on_site, c2.on_site, c3.on_site, c4.on_site = site, site, site, site

            for i in range(min(len(c1), len(c2))):
                c3.insert(c1[i] if i % 2 == 0 and c1[i] not in c2 else c2[i], site)
                c4.insert(c2[i] if i % 2 == 0 and c2[i] not in c1 else c1[i], site)
        
        c1.on_site, c2.on_site, c3.on_site, c4.on_site = None, None, None, None

        #print('INTERLEAVE PROCESS')
        #print('Chromosomes to interleave:')
        #print(f'c1 {c1}')
        #print(f'c2 {c2}')
        #print('New breed chromosomes :')
        #print(f'c3 {c3}')
        #print(f'c4 {c4}')
        self.interleaves += 1
        
        # If we have improvement in one of the crossovers, then replace the lower-valued chromosome.
        max_c = c3 if c3.determine_fitness() > c4.determine_fitness() else c4
        if max_c.determine_fitness() > c2.determine_fitness() and not self._contains(max_c):
            self.chromosomes[c2_index] = max_c
            self._sort_chromosomes()
            self.improve_inter += 1
            return True

        return False

    def _mutation_swap(self):
        """
        Swap two observations in the chromosome.
        """
        c_idx = self._single_selection()
        c = self.chromosomes[c_idx]

        new_c = copy(c)
        for site in {Site.GN,Site.GS}:

            c.on_site, new_c.on_site = site, site
            
            if len(c) < 2:
                return False

            # Sample two observations to swap.
            # This only works if the re-add switches the order.
            pos1, pos2 = sample(range(len(c)), 2)
            pos1, pos2 = (pos1, pos2) if pos1 > pos2 else (pos2, pos1)
            #print(f'swaping pos {pos1} {c[pos1]} for pos {pos2} {c[pos2]}')
            new_c.remove(c[pos1], site)
            #print(f'after removing {c[pos1]} {new_c}')
            new_c.remove(c[pos2], site)
            #print(f'after removing {c[pos2]} {new_c}')
            new_c.insert(c[pos2], site)
            #print(f'after inserting {c[pos2]} {new_c}')
            new_c.insert(c[pos1], site)
            #print(f'after inserting {c[pos1]} {new_c}')
 
            # if new_c.scheduling == c.scheduling:
            #     return False

        c.on_site, new_c.on_site = None, None
        #print("SWAP MUTATOR")
        
        #print(f'old chromosome {c}')
        #print(f'new chromosome {new_c}')
        self.swaps += 1
        #input()
        if new_c.determine_fitness() > c.determine_fitness() and not self._contains(new_c):
            self.chromosomes[c_idx] = new_c
            self._sort_chromosomes()
            self.improve_swap += 1
            return True

        return False

    def _mutation_mix(self) -> bool:
        """
        Try to replace a random number of observations in a randomly selected chromosome.
        """
        c_idx = self._single_selection()
        c = self.chromosomes[c_idx]
        new_c = copy(c)

        for site in {Site.GN,Site.GS}:
            c.on_site, new_c.on_site = site, site
            #print(len(c))
             
            if len(c) <= 1:
                return False
            n = randrange(1, len(c))
            
            # Pick n random observation indices from c to drop.
            obs_idx_to_drop = sorted(sample([obs_idx for _, obs_idx in c.scheduling[site]], n), reverse=True)
            for obs_idx in obs_idx_to_drop:
                new_c.remove(obs_idx, site)

            # Pick n random observation indices to try to insert.
            candidates = [o.idx for o in self.observations if o.site in {site, Site.Both}]
            #print(len(candidates),n)
            #obs_idx_to_add = sample(range(len(candidates)), min(len(candidates), n))
            obs_idx_to_add = sample(candidates,n)
            #print([self.observations[obss].site for obss in obs_idx_to_add])
            #print(site)
            for obs_idx in obs_idx_to_add:
               new_c.insert(obs_idx, site)
        
        c.on_site, new_c.on_site = None, None
        #print("MIX MUTATOR")
        #print(f'old chromosome {c}')
        #print(f'new chromosome {new_c}')
        self.mixes += 1
        if new_c.determine_fitness() > c.determine_fitness() and not self._contains(new_c):
            self.chromosomes[c_idx] = new_c
            self._sort_chromosomes()
            self.improve_mix +=1
            return True

        return False

    def _shuffle(self) -> bool:
        """
        Reorder the observations in a chromosome by adding them to a new chromosome in a random order
        and then seeing if this does better than the original chromosome.
        """
        c_idx = self._single_selection()
        c = self.chromosomes[c_idx]

        if len(c) <= 1:
            return False

        shuffled_obs_idxs = [obs_idx for _, obs_idx in c.scheduling]
        shuffle(shuffled_obs_idxs)

        new_c = Chromosome(self.time_slots, self.observations)
        for obs_idx in shuffled_obs_idxs:
            c.insert(obs_idx)

        if new_c.determine_fitness() > c.determine_fitness() and not self._contains(new_c):
            self.chromosomes[c_idx] = new_c
            self._sort_chromosomes()
            return True

        return False

    def _print_best_fitness(self, sites: List[Site], i: int = None) -> None:
        c = self.chromosomes[0]
        for site in sites:
            print(f"Best fitness for {site.name}{f' iteration {i}' if i is not None else ''}: {c.determine_fitness()} {c.scheduling[site]}")

    def _print_chromosomes(self) -> None:
        for c in self.chromosomes:
            
            sched = []
            print(f'Fitness: {c.determine_fitness()}')
            print('GS')
            print(f'len: {len(c.schedule[Site.GS])}')
            for i in c.schedule[Site.GS]:
                if i:
                    sched.append(colored(i,'grey','on_white'))
                else:
                    sched.append(colored('F','white','on_white'))
            t = colored(',','white','on_white')
            print(f'{t}'.join(sched))
            sched = []
            print('GN')
            print(f'len: {len(c.schedule[Site.GN])}')
            for i in c.schedule[Site.GN]:
                if i:
                    sched.append(colored(i,'grey','on_white'))
                else:
                    sched.append(colored('F','white','on_white'))
            t = colored(',','white','on_white')
            print(f'{t}'.join(sched))
            input()

    def _print_population(self) -> None:

        for i, c in enumerate(self.chromosomes):
            print(f'chromosome {i} fitness: {c.determine_fitness():.6}'+
                 f'GS cap: {c.determine_capacity(Site.GS)}/{self.time_slots.num_time_slots_per_site[Site.GS]} '+
                 f'GN cap: {c.determine_capacity(Site.GN)}/{self.time_slots.num_time_slots_per_site[Site.GN]} '+
                 f'Partials obs: {c.partitions}')
         
            print(f'{c.scheduling}')
        input()

    def _run(self, max_iterations_without_improvement) -> Union[None, Chromosome]:
        """
        The meat of the run algorithm. We do this twice: once to get a chromosome for each site as described in the
        run method.
        """

        

        self._sort_chromosomes()
        #best_c_gs = None if len(self.chromosomes[Site.GS]) == 0 or Site.GS not in sites else copy(self.chromosomes[Site.GS][0])
        #best_c_gn = None if len(self.chromosomes[Site.GN]) == 0 or Site.GN not in sites else copy(self.chromosomes[Site.GN][0])
        best_c = None if len(self.chromosomes) == 0 else copy(self.chromosomes[0])
        sites = {Site.GS, Site.GN}
        # Count the chromosomes at each site.
        #gs_chromosomes = len(self.chromosomes[Site.GS])
        #gn_chromosomes = len(self.chromosomes[Site.GN])
        n_chromosomes = len(self.chromosomes)
        # Perform all of the operations
        counter = 0
        while counter < max_iterations_without_improvement:
            # print('*** START ITERATION ***')
            # print(f'GS chromosomes: {len(self.chromosomes[Site.GS])}')
            # print(f'GN chromosomes: {len(self.chromosomes[Site.GN])}')
            
            #pass
            
            self._mate()
            self._interleave()
            self._mutation_swap()
            self._mutation_mix()
            #self._shuffle()
            #input()
            # print(f'GS chromosomes: {len([c for c in self.chromosomes if c.site == Site.GS])}')
            # print(f'GN chromosomes: {len([c for c in self.chromosomes if c.site == Site.GN])}')
            # print('*** DONE ITERATION ***')

            # See if we have a better best chromosome.
            #new_best_c_gs = None if len(self.chromosomes[Site.GS]) == 0 else copy(self.chromosomes[Site.GS][0])
            #new_best_c_gn = None if len(self.chromosomes[Site.GN]) == 0 else copy(self.chromosomes[Site.GN][0])
            new_best = None if len(self.chromosomes) == 0 else copy(self.chromosomes[0])

            improvement = False
            #if new_best_c_gs is not None and best_c_gs is not None and new_best_c_gs.determine_fitness() > best_c_gs.determine_fitness():
            #    best_c_gs = copy(self.chromosomes[Site.GS][0])
            #    improvement = True
            #if new_best_c_gn is not None and best_c_gn is not None and new_best_c_gn.determine_fitness() > best_c_gn.determine_fitness():
            #    best_c_gn = copy(self.chromosomes[Site.GN][0])
            #    improvement = True

            if new_best is not None and best_c is not None and new_best.determine_fitness() > best_c.determine_fitness():
                best_c = copy(self.chromosomes[0])
                improvement = True

            if improvement:
                self._print_best_fitness(sites, counter)
                counter = 0
            else:
                counter += 1

        # Pick the best chromosome and return it.
        #if best_c_gs is not None and (best_c_gn is None or best_c_gs.determine_fitness() > best_c_gn.determine_fitness()):
        #    return Site.GS, copy(best_c_gs)
        #elif best_c_gn is not None:
        #    return Site.GN, copy(best_c_gn)
        #else:
        #    return None
        
        if best_c is not None:
            return copy(best_c)
        else:
            None

    def run(self, max_iterations_without_improvement=100) -> Schedule:
        """
        Run the genetic algorithm prototype and return the best chromosomes for GN and GS.
        There is a danger that an observation that can be scheduled at both sites will be.
        We cannot allow this, so pick the highest fitness chromosome for the site it represents.
        Change all of its observations to the relevant site.
        Drop all of the observations from the remaining chromosomes.
        Drop all of the chromosomes for the site we covered.
        Iterate more to get a good result for the other site.
        :return: two chromosomes, one for GN and one for GS
        """
        self._form_initial_population()
        self._print_population()
        #ites = ([Site.GS] if len(self.chromosomes[Site.GS]) > 0 else 0) + \
        #        ([Site.GN] if len(self.chromosomes[Site.GN]) > 0 else 0)

        # Check that everything is scheduled.
        obs_in_gs_chroms = len(set([obs_idx for c in self.chromosomes for _, obs_idx in c.scheduling[Site.GS]]))
        obs_in_gn_chroms = len(set([obs_idx for c in self.chromosomes for _, obs_idx in c.scheduling[Site.GN]]))
        gs_obs = len([obs for obs in self.observations if obs.site == Site.GS])
        gn_obs = len([obs for obs in self.observations if obs.site == Site.GN])
        gb_obs = len([obs for obs in self.observations if obs.site == Site.Both])

        #print(f'obs_in_gs_chroms={obs_in_gs_chroms}, gs_obs={gs_obs}, both={gb_obs}, total={gs_obs + gb_obs}')
        #print(f'obs_in_gn_chroms={obs_in_gn_chroms}, gn_obs={gn_obs}, both={gb_obs}, total={gn_obs + gb_obs}')
        #assert(obs_in_gs_chroms == gs_obs + gb_obs)
        #assert(obs_in_gn_chroms == gn_obs + gb_obs)

        # print(f'GS obs: {gs_obs}, GN obs: {gn_obs}, Both: {gb_obs}')
        # print('\nGS Initial Chromosomes')
        # for n, c in enumerate(self.chromosomes[Site.GS]):
        #     print(f'{n}: {c.determine_fitness()} {len([i for i in c.schedule if i is not None])} {c.scheduling}')
        #     # print(f'\t{c.schedule}')
        # print('GN Initial chromosomes:')
        # for n, c in enumerate(self.chromosomes[Site.GN]):
        #     print(f'{n}: {c.determine_fitness()} {len([i for i in c.schedule if i is not None])} {c.scheduling}')
        #     # print(f'\t{c.schedule}')

        #best_c_gn = None
        best_c = self._run(max_iterations_without_improvement)
        
        print('Total operations done'+ 
              f'\nMate\nDone: {self.mates}\tImprove: {self.improve_mate}'+
              f'\nInterleave\nDone: {self.interleaves}\tImprove: {self.improve_inter}'+ 
              f'\nMutator Swap\nDone: {self.swaps}\tImprove: {self.improve_swap}'+
              f'\nMutator Mix\nDone: {self.mixes}\tImprove: {self.improve_mix}')
        #best_c_gs = None

        #results = self._run(max_iterations_without_improvement)
        #if results is None:
        #    return None, None

        if best_c is None:
            return None
        #best_site, best_c = results
        #if best_site == Site.GN:#
        #    best_c_gn = best_#c
        #if best_site == Site.#GS:
        #    best_c_gs = best_#c#

        # print('\nGS Final Chromosomes')
        # for n, c in enumerate(self.chromosomes[Site.GS]):
        #     print(f'{n}: {c.determine_fitness()} {len([i for i in c.schedule if i is not None])} {c.scheduling}')
        #     # print(f'\t{c.schedule}')
        # print('GN Final chromosomes:')
        # for n, c in enumerate(self.chromosomes[Site.GN]):
        #     print(f'{n}: {c.determine_fitness()} {len([i for i in c.schedule if i is not None])} {c.scheduling}')
        #     # print(f'\t{c.schedule}')

        # obs_in_gs_chroms = len(set([obs_idx for c in self.chromosomes[Site.GS] for _, obs_idx in c.scheduling]))
        # obs_in_gn_chroms = len(set([obs_idx for c in self.chromosomes[Site.GN] for _, obs_idx in c.scheduling]))
        # gs_obs = len([obs for obs in self.observations if obs.site == Site.GS])
        # gn_obs = len([obs for obs in self.observations if obs.site == Site.GN])
        # gb_obs = len([obs for obs in self.observations if obs.site == Site.Both])

        # print(f'obs_in_gs_chroms={obs_in_gs_chroms}, gs_obs={gs_obs}, both={gb_obs}, total={gs_obs + gb_obs}')
        # print(f'obs_in_gn_chroms={obs_in_gn_chroms}, gn_obs={gn_obs}, both={gb_obs}, total={gn_obs + gb_obs}')

        # Drop all the observations from the other site that have been scheduled for this site.
        #other_site = Site.GN if best_site == Site.GS else Site.GS
        #for _, obs_idx in best_c.scheduling:
        #    for c in self.chromosomes[other_site]:
        #        c.remove(obs_idx)

        # Remove any blank chromosomes.
        self.chromosomes = [o for o in self.chromosomes if len(o.scheduling) > 0]

        # print(f'**** NEW CHROMOSOMES ****')
        # for n, c in enumerate([c for c in self.chromosomes[other_site]]):
        #     print(f'{n}: {len([i for i in c.schedule if i is not None])} {c.scheduling}')
        #     # print(f'\t{c.schedule}')

        # Now repeat the process if chromosomes are left. All that is left are chromosomes from the other site.
        #if len(self.chromosomes) > 0:
        #    results = self._run(max_iterations_without_improvement)
        #    if results is not None:
        #        best_site, best_c = results
        #        if best_site == Site.GS:
        #            best_c_gs = best_c
        #        else:
        #            best_c_gn = best_c

        # If either is still None, return an empty schedule.
        #best_gs = best_c_gs.schedule if best_c_gs is not None else [None] * self.time_slots.num_time_slots_per_site[Site.GS]
        #best_gn = best_c_gn.schedule if best_c_gn is not None else [None] * self.time_slots.num_time_slots_per_site[Site.GN]
        best_c = best_c.schedule if best_c is not None else [None] * (self.time_slots.num_time_slots_per_site[Site.GN]+
                                                                        self.time_slots.num_time_slots_per_site[Site.GS])
        return best_c

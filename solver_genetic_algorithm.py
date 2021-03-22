# solver_generic_algorithm.py
# By Sebastian Raaphorst, 2020.

from input_parser import read_tables
from ga_solver import *
import output
from time import monotonic

if __name__ == '__main__':
    tabdir = './data/'
    prefix = ''
    suffix = '_gngs_20201116'
    time_slots, observations = read_tables(f'{tabdir}{prefix}timetab{suffix}.fits',
                                           f'{tabdir}{prefix}obstab{suffix}.fits',
                                           f'{tabdir}{prefix}targtab_metvisha{suffix}.fits')

    start_time = monotonic()
    num_runs = 10
    tot_score = 0
    best_score = 0
    best_schedule = None
    schedule = None
    for i in range(num_runs):
        ga = GeneticAlgorithm(time_slots, observations, include_greedy_max=False)
        schedule = ga.run()
        gs_score = calculate_scheduling_score(Site.GS, time_slots, observations,
                                              output.convert_to_scheduling(schedule))
        gn_score = calculate_scheduling_score(Site.GN, time_slots, observations,
                                              output.convert_to_scheduling(schedule))
        score = gn_score + gs_score
        # We have to schedule these pairwise, because otherwise we may get schedules where an observation is
        # scheduled at both sites: thus, for scoring, use the sum of the scores of the schedules to determine
        # the best schedule. We may move to a single pool of longer chromosomes similar to how we do things for ILP
        # instead of the way we do things now for GA with a pool for GS and GN, although this will making handling
        # Site.both observations more challenging. It will, on the other hand, simplify AND / OR.
        print(f'Final score: {score}')
        input()
        #print(f'Final GN score: {gn_score}')
        if score > best_score:
            best_schedule = schedule
            best_score = score
            #print(f'Best schedule now: {best_schedule}')
            print(f'Best schedule now:')
            output.in_line_print(best_schedule)
            input()
        tot_score += score

    total_time = monotonic() - start_time
    print('\n\n*** RESULTS ***')
    print(f'Number of runs: {num_runs}')
    print(f'Average score: {tot_score / num_runs}')
    #print(f'Average GN score: {gn_tot_score / num_runs}')
    #print(f'Average total score: {(gn_tot_score + gs_tot_score) / num_runs}')
    print(f'Time taken: {total_time} s')
    
    print(f'Average time per run: {total_time / num_runs} s')
    
    #print(schedule)
    #print(best_schedule)
    #print('************************************************')

    # *** DELETE ***
    # gn_scheduling = output.convert_to_scheduling(gn_schedule)
    # print(f'GN fitness: {output.calculate_score(time_slots, observations, gn_scheduling)}')
    # print(f'GN: {gn_schedule}')
    # 
    # gs_scheduling = output.convert_to_scheduling(gs_schedule)
    # print(f'GS fitness: {output.calculate_score(time_slots, observations, gs_scheduling)}')
    # print(f'GS: {gs_schedule}')
    # *** END DELETE ***

    # Once print_schedule is implemented, delete the indicated area.
    print('\n\n')
    output.print_schedule3(time_slots, observations, best_schedule)
    
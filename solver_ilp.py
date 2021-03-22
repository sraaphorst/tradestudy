# solver_ilp.py
# By Sebastian Raaphorst, 2020.

from input_parser import read_tables
from ilp_solver import *
import output
from time import monotonic

if __name__ == '__main__':
    tabdir = './data/'
    prefix = ''
    suffix = '_gngs_20201116'
    time_slots, observations = read_tables(f'{tabdir}{prefix}timetab{suffix}.fits',
                                           f'{tabdir}{prefix}obstab{suffix}.fits',
                                           f'{tabdir}{prefix}targtab_metvisha{suffix}.fits')

    # Run the solver.
    start_time = monotonic()
    gs_schedule, gn_schedule = ilp_scheduler(time_slots, observations)
    total_time = monotonic() - start_time
    schedule = [gn_schedule, gs_schedule]

    # *** DELETE ***
    gn_scheduling = output.convert_to_scheduling(Site.GN, schedule)
    print(f'GN fitness: {output.calculate_scheduling_score(Site.GN, time_slots, observations, gn_scheduling)}')
    print(f'GN: {gn_schedule}')

    gs_scheduling = output.convert_to_scheduling(Site.GS, schedule)
    print(f'GS fitness: {output.calculate_scheduling_score(Site.GS, time_slots, observations, gs_scheduling)}')
    print(f'GS: {gs_schedule}')
    # *** END DELETE ***

    # Once print_schedule is implemented, delete the indicated area.
    print('\n\n*** RESULTS ***')
    print(f'Time taken: {total_time} s')
    output.print_schedule3(time_slots, observations, schedule)


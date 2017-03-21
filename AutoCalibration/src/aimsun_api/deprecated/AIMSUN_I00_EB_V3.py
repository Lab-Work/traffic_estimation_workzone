__author__ = 'Carlos and Yanning'
"""
This script is used to validate calibration for a synthetic work zone I00. It can also be used for visualizing how parameters
influencing the calibration errors, and sensitivity analysis.
"""


from AIMSUNFUNCTIONS_V3 import *
from copy import deepcopy
import numpy as np


# ======================================================================================
# Configurations in this block
# ======================================================================================

# The communication between the optimization solver and AIMSUN is via files:
sim_para_file = 'C:/sim_com/sim_paras.txt'  # the parameter file
sim_val_file = 'C:/sim_com/sim_val.txt'     # the objective value file
sim_sol_file = 'C:/sim_com/sim_sol.txt'     # the final solution file. Create this to stop AIMSUN from simulation

# number of replications and associated seeds
g_num_rep = 10
seed_list = [3503, 23798, 28860, 12358, 6370, 14452, 21488, 7893, 586, 30447]

# g_num_rep = 20
# seed_list = [3503, 23798, 28860, 12358, 6370, 14452, 21488, 7893, 586, 30447,
#              3645, 22234, 25345, 12442, 6232, 14623, 21525, 7092, 523, 30194]

# This is the file that saves the evaluated parameters, so no need to re-simulate
# opt_step_file = 'c:/Users/TrafficControl/Google Drive/AIMSUN_carlos/Workzone_calibration/I00_EB/congflow/Logs/test_synthetic_with_error_model/20160516_125457_opt_log.csv'
# opt_step_file = 'c:/Users/TrafficControl/Google Drive/AIMSUN_carlos/Workzone_calibration/I00_EB/congflow/Logs/20160513_170644_opt_log.csv'
# opt_step_file = 'c:/Users/TrafficControl/Google Drive/AIMSUN_carlos/Workzone_calibration/I00_EB/congflow/Logs/20160602_151254_opt_log_1_1254.csv'
opt_step_file = None

# -----------------------------------------------------------------------
# Change the following to be consistent with the ang network file.
all_det_strs =   ['EB1','EB2','EB3','EB4']

# Make sure OptQuest ends first. This program ends once it found sim_sol.txt file
# This value is used for calculating the computation time
g_max_iter = 4000

# Define the way to compute the objective function
# obj_fun[0]*RMS_speed + obj_fun[1]*RMS_count
obj_fun = [1, 0]

# detector used for validation
det_used = all_det_strs
det_used_weight = {'EB1':1.0, 'EB2':1.0, 'EB3':1.0, 'EB4':1.0}

# The following noise model is used to generate validation data.
speed_noise_model = {}
speed_noise_model['EB1'] = [-2, 2]  # [bias, std]
speed_noise_model['EB2'] = [+2, 2]  # [bias, std]
speed_noise_model['EB3'] = [-2, 2]  # [bias, std]
speed_noise_model['EB4'] = [+2, 2]  # [bias, std]

count_noise_model = {}
count_noise_model['EB1'] = [0, 20] # std: 20 cnt/5min = 240 count/hr = 10% if 2400 flow
count_noise_model['EB2'] = [0, 20] # std: 20 cnt/5min = 240 count/hr = 10% if 2400 flow
count_noise_model['EB3'] = [0, 20] # std: 20 cnt/5min = 240 count/hr = 10% if 2400 flow
count_noise_model['EB4'] = [0, 20] # std: 20 cnt/5min = 240 count/hr = 10% if 2400 flow


# Files used in this file
folder_name = 'congflow'
# this is the traffic state names that should be used in .ang
# demand name: traffic_state + '_demand'
# scenario name: traffic_state
# experiment name: traffic_state + '_exp'
traffic_state = 'congflow'

parsFilePath = 'c:/Users/TrafficControl/Dropbox/DesktopServer/Workzone_autocalibration/I00_EB/'+ folder_name +'/demand_data/paras.txt'
flowsFilePath = 'c:/Users/TrafficControl/Dropbox/DesktopServer/Workzone_autocalibration/I00_EB/'+ folder_name +'/demand_data/flows.txt'
turnsFilePath = 'c:/Users/TrafficControl/Dropbox/DesktopServer/Workzone_autocalibration/I00_EB/'+ folder_name +'/demand_data/turns.txt'

presetParaFile = 'c:/Users/TrafficControl/Dropbox/DesktopServer/Workzone_autocalibration/I00_EB/'+ folder_name +'/preset_paras.txt'
validFilePath = 'c:/Users/TrafficControl/Dropbox/DesktopServer/Workzone_autocalibration/I00_EB/'+ folder_name +'/validation_data/'  # in current folder
logger_path = 'c:/Users/TrafficControl/Dropbox/DesktopServer/Workzone_autocalibration/I00_EB/'+ folder_name +'/Logs/'


# the following two time strings are used to import the validation data in the selected time period
# They should be the exact start and end time of the simulated period (e.g. fist 15 cong states:  05/01/2015 15:30 ~ 05/01/2015 16:45)
# start_time_str = '05/01/2015 15:30'
# end_time_str = '05/01/2015 18:00'
start_time_str = None
end_time_str = None

main_entrance_id = 329
main_entrance_discount = 1.0



# Put one line of description of the simulation
description = 'This is the synthetic calibration of the {0} state on I00.\n'.format(traffic_state)

#============================================================================
# set up default parameters, and true parameters if known
# Note, here the default parameters only focuses on the values to be calibrated; others will be static in preset_paras
simulate_default = False     # simulate default parameters
simulate_user = True        # simulate the user set parameters
simulate_true = True         # simulate the true parameters, save the true result
simulate_newseed = False     # simulate the newseed with the obtained optimal parameters or the default parameters in the
                            # synthetic validataion.

# ------------------------------------------------------------------
default_paras = {}

default_paras['car_maxAccel'] = [3, 0.2, 2.6, 3.4]
default_paras['car_reactionTime'] = [0.8, 1.2, 1.6, 1]
default_paras['car_minDist'] = [1, 0.3, 0.5, 1.5]       # range in OptQuest [0.5, 10]
default_paras['car_sensitivityFactor'] = [1, 0, 1, 1]
default_paras['truck_maxAccel'] = [1, 0.5, 0.6, 1.8]
default_paras['truck_reactionTime'] = [0.8, 1.3, 1.7, 1]
default_paras['truck_minDist'] = [1.5, 0.5, 1.0, 2.5]   # range in OptQuest [1, 15]
default_paras['truck_sensitivityFactor'] = [1, 0, 1, 1]
#
default_paras['car_speedAcceptance'] = [1.1, 0.1, 0.85, 1.3]
default_paras['truck_speedAcceptance'] = [1.05, 0.1, 0.85, 1.1]

# ------------------------------------------------------------------
user_paras = {}

user_paras['car_speedAcceptance'] = [1.29, 0.1, 0.85, 1.3]
user_paras['truck_maxAccel'] = [1.18, 0.5, 0.6, 1.8]
user_paras['car_sensitivityFactor'] = [0.88, 0, 0.8, 0.8]
user_paras['truck_sensitivityFactor'] = [1.1, 0, 0.85, 0.85]
user_paras['car_reactionTime'] = [1.0, 1.2, 1.6, 1]
user_paras['truck_reactionTime'] = [1.0, 1.3, 1.7, 1]

user_paras['car_maxAccel'] = [3.0, 0.2, 2.6, 3.4]
user_paras['car_minDist'] = [1.0, 0.3, 0.5, 1.5]
user_paras['truck_minDist'] = [1.5, 0.5, 1.0, 2.5]
user_paras['truck_speedAcceptance'] = [1.05, 0.1, 0.85, 1.1]

# ------------------------------------------------------------------
true_paras = {}
# changed from the default paras, all steps are 0.01, except for reactionTime is 0.1
# In toatal 4.2*10^12 combinations
"""
Additional constraints on the calibrated parameters:
1. mean car_speedAcceptance >= mean default_truck_speedAcceptance
2. car_sensitivityFactor <= truck_sensitivityFactor ===> Cars more aggressive.
3. change the simulation step from 0.8 to 0.2, hence finer grid
"""
true_paras['car_speedAcceptance'] = [1.08, 0.1, 0.85, 1.3]
true_paras['truck_maxAccel'] = [0.9, 0.5, 0.6, 1.8]     # range: [0.6, 1.8], # 121
true_paras['car_sensitivityFactor'] = [0.95, 0, 1.05, 1.05]     # range: [0.9, 1.1], # 21
true_paras['truck_sensitivityFactor'] = [1.03, 0, 1.03, 1.03]   # range: [0.9, 1.1], # 21
true_paras['car_reactionTime'] = [0.8, 1.2, 1.6, 1]     # range: {0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6}
true_paras['truck_reactionTime'] = [1.2, 1.3, 1.7, 1]   # range: {0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6}

# same as default, not changed
true_paras['car_maxAccel'] = [3.0, 0.2, 2.6, 3.4]   # range: [2.6, 3.4], # 81
true_paras['car_minDist'] = [1.0, 0.3, 0.5, 1.5]   # range: [0.5, 1.5], # 101
true_paras['truck_minDist'] = [1.5, 0.5, 1.0, 2.5]     # range: [1.0, 2.5], # 151
true_paras['truck_speedAcceptance'] = [1.05, 0.1, 0.85, 1.1]

# true_paras = {}
# true_paras['car_speedAcceptance'] = [1.05, 0.1, 0.85, 1.3]  # changed
# true_paras['car_maxAccel'] = [3.2, 0.2, 2.6, 3.4]           # changed
# true_paras['car_reactionTime'] = [0.8, 1.2, 1.6, 1]
# true_paras['car_minDist'] = [2, 0.3, 0.5, 3]       # changed
# true_paras['car_sensitivityFactor'] = [0.9, 0, 0.9, 0.9]    # changed
# true_paras['truck_speedAcceptance'] = [0.9, 0.1, 0.85, 1.1] # changed
# true_paras['truck_maxAccel'] = [1.2, 0.5, 0.6, 1.8]         # changed
# true_paras['truck_reactionTime'] = [0.8, 1.3, 1.7, 1]
# true_paras['truck_minDist'] = [3, 0.5, 1.0, 5]   # changed
# true_paras['truck_sensitivityFactor'] = [0.95, 0, 0.95, 0.95]


# ======================================================================================
# End of configurations in this block
# ======================================================================================


# Load a network using ANGConsole
def main(argv):

    # a timeout flag. If not paras or solution paras are obtained in 60 s. Then stop AIMSUN.
    timeout_flag = False
    timeout_counter = 0

    solved_flag = False

    # keep track of the current best parameter and its associated objective
    best_para = deepcopy(default_paras)
    best_obj = 0

    # keep track of start time:
    start_time = datetime.now()

    # Start the cmd output logger
    start_cmd_logger(logger_path, start_time)
    start_opt_logger(logger_path, start_time)

    # print python version
    # print '\nPython interpreter version: {0} \n'.format(sys.version)

    print_cmd('\n\n========================================================================')
    print_cmd('=========================AIMSUN Autocalibration=========================')
    print_cmd('========================================================================')
    print_cmd('Auto calibration started at {0}\n'.format(datetime.now().strftime("%Y-%m-%d %H:%M:%S")))

    # write notes and parameters to cmd_log
    cmd_logger_header(description,
                      g_max_iter, g_num_rep, seed_list,
                      obj_fun,
                      det_used,
                      main_entrance_id, main_entrance_discount,
                      [('default', default_paras), ('true', true_paras)])

    print_cmd('Sensor noise models are:\n---- speed: {0}\n---- count: {1}\n'.format(speed_noise_model, count_noise_model))

    if (len(argv) < 2):
        print_cmd('Usage: aconsole.exe -script SCRIPT ANG_FILE')
        return -1
    else:
        # Start the Aimsun console
        console = ANGConsole()
        if console.open(argv[1]):

            print_cmd('\nAimsun opening {0} ...\n'.format(argv[1]))

            #========================================================
            # Read demand
            # Get the current Aimsun  model
            model = GKSystem.getSystem().getActiveModel()

            demand = read_demand_from_file(model, traffic_state + '_demand',
                                           parsFilePath, flowsFilePath, turnsFilePath,
                                           main_entrance_id, main_entrance_discount)
            print_cmd('Finished setting Demand. \n')

            #=====================================================================================
            # Setup model
            #=====================================================================================
            # create scenario
            # if exists, then just load
            scenario = setup_scenario(model, traffic_state, demand)

            # create experiment
            experiment = setup_experiment(model, traffic_state + '_exp', scenario)

            # create replications
            avg_result = setup_replication(model, experiment, g_num_rep, seed_list)

            # create simulator
            simulator = create_simulator(model)
            # plugin is a module which can compute the average for the GKExperimentResult object
            plugin = GKSystem.getSystem().getPlugin( "GGetram" )


            #=====================================================================================
            # validation data for computing the objective
            #=====================================================================================
            # ----read from file
            # valid_data = read_validation_data(start_time_str, end_time_str, det_used, validFilePath)
            # ----or generate

            set_new_paras(model, experiment, true_paras)
            # save
            console.save("visual_I00_EB_valid.ang")
            print_cmd('\nvisual_I00_EB_valid.ang saved.')

            valid_result = simulate_rep_from_paras(model, experiment, true_paras,
                                                  simulator, avg_result, plugin,
                                                   None, det_used_weight,
                                                   'true')
            clean_meas_data = valid_result[1]
            valid_data = clean_meas_data
            # valid_data = add_noise_to_data( clean_meas_data, speed_noise_model, count_noise_model)

            # we keep the single loop run time for evaluating the time needed for simulation

            # save the complete validation data used for visualization
            # all_valid_data = read_validation_data(start_time_str, end_time_str, all_det_strs, validFilePath)
            save_solution_data(true_paras, clean_meas_data, logger_path, start_time, 'clean')
            save_solution_data(true_paras, valid_data, logger_path, start_time, 'valid')


            #=====================================================================================
            # load previous solution
            #=====================================================================================
            if opt_step_file is not None:
                solution_list = load_previous_opt_solutions(opt_step_file)
                print_cmd('\nLoaded previous OptQuest solutions:\n---- {0}'.format(opt_step_file))

                # get the current best parameters
                previous_obj_tuple = np.array( solution_list[1] )
                previous_obj_val = obj_fun[0]*previous_obj_tuple[:,0] + obj_fun[1]*previous_obj_tuple[:,1]

                idx = np.argmin( previous_obj_val )

                best_para = deepcopy(solution_list[0][idx])
                best_obj = deepcopy(solution_list[1][idx])

                print_cmd('\nBest obj from loaded parameters: {0}\n'.format( obj_fun[0]*best_obj[0] +
                                                                           obj_fun[1]*best_obj[1] ))

            else:
                solution_list = None


            #=====================================================================================
            # generate the result using default values
            #=====================================================================================
            if simulate_default is True:

                set_new_paras(model, experiment, default_paras)
                # test
                console.save("visual_I00_EB_default.ang")
                print_cmd('\nvisual_I00_EB_default.ang saved.')

                default_result = simulate_rep_from_paras(model, experiment, default_paras,
                                                         simulator, avg_result, plugin,
                                                         valid_data, det_used_weight,
                                                         'default')

                save_solution_data(default_paras, default_result[1], logger_path, start_time, 'default')

                #========================================================
                # For printing result
                # compute the objective function value if using the true parameters
                # true_obj_value = valid_result[0]
                default_obj_val = default_result[0]

                best_obj = deepcopy(default_obj_val)

            else:
                default_obj_val = (-1, -1)


            #=====================================================================================
            # generate the result using default values
            #=====================================================================================
            if simulate_user is True:

                set_new_paras(model, experiment, user_paras)
                # test
                console.save("visual_I00_EB_user.ang")
                print_cmd('\nvisual_I00_EB_user.ang saved.')

                user_result = simulate_rep_from_paras(model, experiment, user_paras,
                                                         simulator, avg_result, plugin,
                                                         valid_data, det_used_weight,
                                                         'user')

                save_solution_data(user_paras, user_result[1], logger_path, start_time, 'user')

                #========================================================
                # For printing result
                # compute the objective function value if using the true parameters
                # true_obj_value = valid_result[0]
                user_obj_val = user_result[0]

            else:
                user_obj_val = (-1, -1)


            #=====================================================================================
            # save the true simulation result:
            #=====================================================================================
            if simulate_true is True:

                result = simulate_rep_from_paras(model, experiment, true_paras,
                                                 simulator, avg_result, plugin,
                                                 valid_data, det_used_weight, 'iteration True')
                true_obj_val = result[0]
            else:
                true_obj_val = (0, 0)



            #=====================================================================================
            # Simulating parameters starts here
            #=====================================================================================
            opt_solution = []

            iter_counter = 1
            simulated_iter_counter = 1

            # g_max_iter + 1 is to make sure OptQuest stops first
            while iter_counter <= g_max_iter + 1 and solved_flag is False:

                run_time = datetime.now() - start_time
                expected_finish_time = datetime.now() + run_time*(g_max_iter + 3 - iter_counter)/simulated_iter_counter

                print_cmd('\n--------------------Iteration {0}------------------------'.format(iter_counter))

                # print out run time information
                print_cmd('Started simulation at: {0}'.format(start_time.strftime("%Y-%m-%d %H:%M:%S")))
                print_cmd('Simulation run time  : {0}'.format(str(run_time)))
                print_cmd('Expected to finish at: {0}\n'.format(expected_finish_time.strftime("%Y-%m-%d %H:%M:%S")))

                #=====================================================================================
                # Read parameters
                #=====================================================================================
                # first check if solution has been found
                paras = read_optquest_solution(sim_sol_file)

                if paras is None:
                    # solution not found; keep going
                    # read initial/new parameters from OptQuest
                    paras = read_optquest_paras(sim_para_file)

                    # if paras is still None, then probably a solution has been written
                    if paras is None:

                        if timeout_flag is False:
                            if timeout_counter >= 3:
                                timeout_flag = True
                                break
                            else:
                                timeout_counter += 1
                                continue

                    else:
                        # wait a few cycles and got new paras, reset timeout
                        timeout_counter = 0
                        timeout_flag = False

                else:
                    # mark as solved
                    print_cmd('\n------------------Optimal solution----------------------')
                    solved_flag = True

                #=====================================================================================
                # Simulate and get the objective function value
                #=====================================================================================
                # first check previous solutions
                tmp_obj_val = try_get_obj_from_previous_solutions(solution_list, paras)

                if tmp_obj_val is None:

                    # the true parameters are the default.
                    print_cmd('Reset default paras:')
                    set_new_paras(model, experiment, true_paras)
                    # overwrite a subset of the parameters from paras
                    print_cmd('Overwrite paras:')
                    set_new_paras(model, experiment, paras)

                    # no solution, then simulate
                    result = simulate_rep_from_paras(model, experiment, paras,
                                                         simulator, avg_result, plugin,
                                                         valid_data, det_used_weight,
                                                         'iteration {0}'.format(iter_counter))

                    # evaluate the objective function value
                    obj_value = result[0]

                    # truly simulated number
                    simulated_iter_counter += 1

                else:
                    obj_value = tmp_obj_val
                    print_cmd('Got obj value {0} from previous solutions'.format(tmp_obj_val))

                #=====================================================================================
                # send simulation value to OptQuest
                #=====================================================================================
                write_simval(obj_fun[0]*obj_value[0] + obj_fun[1]*obj_value[1], sim_val_file)


                #=====================================================================================
                # update the current best value
                #=====================================================================================
                if obj_fun[0]*obj_value[0] + obj_fun[1]*obj_value[1] <= \
                    obj_fun[0]*best_obj[0] + obj_fun[1]*best_obj[1] :
                    # found better parameters
                    best_para = deepcopy(paras)
                    best_obj = deepcopy(obj_value)

                #=====================================================================================
                # Log the result
                log_opt_step(opt_solution, iter_counter, paras, obj_value)

                print_opt_steps(opt_solution, obj_fun, [('default', default_obj_val), ('true', true_obj_val)])

                iter_counter += 1

            # while loop ends here
            #========================================================

            # save logger to file
            stop_opt_log()

            # save the optimal simulation result:
            if solved_flag is True:
                optimal_paras = paras

                set_new_paras(model, experiment, optimal_paras)
                # save the ang file
                console.save("visual_I80_EB_optimal.ang")
                print_cmd('\nvisual_I80_EB_optimal.ang saved.')

                # re-simulate, in case no data available
                result = simulate_rep_from_paras(model, experiment, optimal_paras,
                                                         simulator, avg_result, plugin,
                                                         valid_data, det_used_weight,
                                                         'iteration Optimal')

                optimal_obj_val = result[0]
                optimal_data = result[1]
                save_solution_data(optimal_paras, optimal_data, logger_path, start_time, 'optimal')

                #========================================================
                # see how the seed affects the optimal solution
                if simulate_newseed is True:
                    set_random_seed(avg_result)
                    newseed_result = simulate_rep_from_paras(model, experiment, true_paras,
                                                             simulator, avg_result, plugin,
                                                             valid_data, det_used_weight,
                                                             'newseed')
                    save_solution_data(true_paras, newseed_result[1], logger_path, start_time, 'newseed')
                    newseed_obj_val = newseed_result[0]
                else:
                    newseed_obj_val = (-1, -1)

            end_time = datetime.now()

            # if timeout, then print out the optimal parameter so far
            if timeout_flag is True:
                set_random_seed(avg_result)
                newseed_result = simulate_rep_from_paras(model, experiment, true_paras,
                                                             simulator, avg_result, plugin,
                                                             valid_data, det_used_weight,
                                                             'newseed')
                save_solution_data(true_paras, newseed_result[1], logger_path, start_time, 'newseed')
                newseed_obj_val = newseed_result[0]

                print_cmd(  '\n\n==================== Calibration timeout ===================================')
                print_cmd(  '===========================================================================')

                print_results([ ('default', default_paras),
                                ('true', true_paras),
                                ('newseed', true_paras),
                                ('optimal', best_para)],
                              [ ('default', default_obj_val),
                                ('true', [0,0]),
                                ('newseed', newseed_obj_val),
                                ('optimal', best_obj)])

            if solved_flag is True:

                print_cmd(  '\n\n==================== Calibration finished ===================================')
                print_cmd(  '===========================================================================')

                print_results([ ('optimal', optimal_paras),
                                ('default', default_paras),
                                ('true', true_paras)],
                              [ ('optimal',optimal_obj_val),
                                ('default', default_obj_val),
                                ('true', true_obj_val),
                                ('newseed', newseed_obj_val)])


            print_cmd(  '\nCalibration started at:    {0}'.format(start_time.strftime("%Y-%m-%d %H:%M:%S")))
            print_cmd(  '              ended at:      {0}'.format(end_time.strftime("%Y-%m-%d %H:%M:%S")))
            print_cmd(  '         took in total:      {0}'.format(str(end_time - start_time)))
            print_cmd(  '===========================================================================')


            print_cmd(  '\nAimsun is now closing...\n')

            # stop logger and save into file
            stop_cmd_logger()

            console.close()

        else:
            console.getLog().addError("Could not open")


if __name__ == "__main__":
    sys.exit(main(sys.argv))
















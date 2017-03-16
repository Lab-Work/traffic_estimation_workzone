__author__ = 'Carlos and Yanning'

# This code is used for automating AIMSUN and interacting with the optimization solver.
# All functions are in AIMSUNFUNCTIONS_V3


from AIMSUNFUNCTIONS_V3 import *
from copy import deepcopy

# ======================================================================================
# Configurations in this block
# ======================================================================================

# ---------- File directories ----------
# Use the following file for communication with the optimization program
sim_para_file = 'C:/thread2_sim_com/sim_paras.txt'
sim_val_file = 'C:/thread2_sim_com/sim_val.txt'
sim_sol_file = 'C:/thread2_sim_com/sim_sol.txt'

# names that should be used in .ang
#  - demand name: traffic_state + '_demand'
#  - scenario name: traffic_state
#  - experiment name: traffic_state + '_exp'
traffic_state = 'congflow'
validFilePath = 'C:/Users/TrafficControl/Dropbox/DesktopServer/Workzone_autocalibration/I80_EB/validation_data/'
logger_path = 'C:/Users/TrafficControl/Dropbox/DesktopServer/Workzone_autocalibration/I80_EB/Logs/thread2/'

# the file that saved previous simulation result for each optimization step.
opt_step_file = 'C:/Users/TrafficControl/Dropbox/DesktopServer/Workzone_autocalibration/I80_EB/Logs/thread2_previous_sim_result.csv'
# opt_step_file = None

# ---------- Specify demand ----------
# Demand can be specified in two ways: traffic states or OD matrix.
#   - For traffic states: use function read_demand_from_file with three files as inputs
#   - For OD matrix: use built in script in AIMSUN and load via GUI. Then in this script,
#     use function load_demand_from_ang() with the demand name imported in AIMSUN.
# deprecated
main_entrance_id = 21216
main_entrance_discount = 1.0

# ---------- Configure seeds ----------
# number of replications. We are not calibrating a work zone. We are REPLICATING a work zone.
# Hence, one replication is sufficient.
g_num_rep = 1
seed_list = [3503]
# seed_list = [33, 23798, 28860, 12358, 6370, 14452, 21488, 7893, 586, 30447]

# make sure OptQuest ends first. This program ends once it found sim_sol.txt file
# this value is used for calculating the computation time
g_max_iter = 1000

# ---------- Detector list for computing the error ----------
# EB3 was used for generating the main inflow, hence do not use as the validation data.
# smart_det_strs = ['EB5','EB7','EB9','EB12']
# radar_det_strs = ['EB4','EB6','EB8','EB10','EB14','EB15','EB16']
all_det_strs =   ['EB4','EB5','EB6','EB7','EB8','EB9','EB10', 'EB11', 'EB12','EB14','EB15','EB16']

# detector used for validation and their associated weights
# detector EB14 is right before the work zone, i.e., where merging occurs, hence less weight.
det_used = ['EB4','EB5','EB6','EB7','EB8','EB9','EB10', 'EB11', 'EB12', 'EB14']
det_used_weight = {'EB4':1.0, 'EB5':1.0, 'EB6':1.0, 'EB7':1.0, 'EB8':1.0,
                   'EB9':1.0, 'EB10':1.0, 'EB11':1.0, 'EB12':1.0, 'EB14':0.5}

# Define the way to compute the objective function
# obj_fun[0]*RMS_speed + obj_fun[1]*RMS_count
obj_fun = [1, 0]

# ------------ For processing the validation data ------------
# the following two time strings are used to import the validation data in the selected time period
# They should be the exact start and end time of the simulated period (e.g. fist 15 cong states:  05/01/2015 15:30 ~ 05/01/2015 16:45)
start_time_str = '05/01/2015 15:30'
end_time_str = '05/01/2015 18:00'

# Put one line of description of the simulation
description = 'This is the calibration of the {0} state. Simulate time from {1} to {2}\n'.format(traffic_state,
                                                                                                 start_time_str,
                                                                                                 end_time_str) + \
              'The traffic demand is manually tuned.'


# ======================================================================================
# The parameters for simulation
# ======================================================================================
# set up the simulation flag
simulate_default = True
simulate_user = False
simulate_newseed = True

# ------------ Default parameters ------------
# set up default parameters
default_paras = OrderedDict()
default_paras['car_speedAcceptance'] = [1.1, 0.1, 0.85, 1.3]
default_paras['truck_maxAccel'] = [1, 0.5, 0.6, 1.8]
default_paras['car_sensitivityFactor'] = [1, 0, 1, 1]
default_paras['truck_sensitivityFactor'] = [1, 0, 1, 1]
default_paras['car_reactionTime'] = [0.8, 1.2, 1.6, 1]
default_paras['truck_reactionTime'] = [0.8, 1.3, 1.7, 1]
default_paras['car_minHeadway'] = [0, 0, 0, 0]
default_paras['truck_minHeadway'] = [0, 0, 0, 0]

default_paras['car_maxAccel'] = [3, 0.2, 2.6, 3.4]
default_paras['car_minDist'] = [1, 0.3, 0.5, 1.5]       # range in OptQuest [0.5, 10]
default_paras['truck_speedAcceptance'] = [1.05, 0.1, 0.85, 1.1]
default_paras['truck_minDist'] = [1.5, 0.5, 1.0, 2.5]   # range in OptQuest [1, 15]


# ------------ User parameters ------------
# specify a set of parameters the user would like to simulate and get the simulation result
# This is the optimal parameters we calibrated using the fine tuned parameters
user_paras = OrderedDict()
user_paras['car_speedAcceptance'] = [0.96, 0.1, 0.85, 1.3]  # calibrated value from freeflow
user_paras['truck_maxAccel'] = [1.58, 0.5, 0.6, 1.8]
user_paras['car_sensitivityFactor'] = [0.8, 0, 0.8, 0.8]
user_paras['truck_sensitivityFactor'] = [0.85, 0, 0.85, 0.85]
user_paras['car_reactionTime'] = [0.8, 1.2, 1.6, 1]
user_paras['truck_reactionTime'] = [0.8, 1.3, 1.7, 1]
user_paras['car_minHeadway'] = [1, 0, 1, 1]
user_paras['truck_minHeadway'] = [2, 0, 2, 2]

user_paras['car_maxAccel'] = [3.02, 0.2, 2.6, 3.4]
user_paras['car_minDist'] = [3.4, 0.3, 0.5, 6]     # range in OptQuest [0.3, 6]
user_paras['truck_speedAcceptance'] = [0.99, 0.1, 0.85, 1.1]
user_paras['truck_minDist'] = [2.8, 0.5, 1, 10]


# deprecated
# parsFilePath = 'c:/Users/TrafficControl/Google Drive/AIMSUN_carlos/Workzone_calibration/I80_EB/'+ folder_name +'/demand_data/paras.txt'
# flowsFilePath = 'c:/Users/TrafficControl/Google Drive/AIMSUN_carlos/Workzone_calibration/I80_EB/'+ folder_name +'/demand_data/flows.txt'
# turnsFilePath = 'c:/Users/TrafficControl/Google Drive/AIMSUN_carlos/Workzone_calibration/I80_EB/'+ folder_name +'/demand_data/turns.txt'
# presetParaFile = 'c:/Users/TrafficControl/Google Drive/AIMSUN_carlos/Workzone_calibration/I80_EB/'+ folder_name +'/preset_paras.txt'


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
    print_cmd('==============================Thread 2==================================')
    print_cmd('Auto calibration started at {0}\n'.format(datetime.now().strftime("%Y-%m-%d %H:%M:%S")))

    # write notes and parameters to cmd_log
    cmd_logger_header(description,
                      g_max_iter, g_num_rep, seed_list,
                      obj_fun,
                      det_used,
                      main_entrance_id, main_entrance_discount,
                      [('default', default_paras), ('user', user_paras)])

    #========================================================
    # validation data for computing the objective
    # ----read from file
    valid_data = read_validation_data(start_time_str, end_time_str, det_used, validFilePath)

    # save the complete validation data used for visualization
    all_valid_data = read_validation_data(start_time_str, end_time_str, all_det_strs, validFilePath)
    save_solution_data(None, all_valid_data, logger_path, start_time, 'valid')


    #========================================================
    # load previous solution
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


    #========================================================
    # start AIMSUN simulation
    if (len(argv) < 2):
        print_cmd('Usage: aconsole.exe -script SCRIPT ANG_FILE')
        return -1
    else:
        # Start the Aimsun console
        console = ANGConsole()
        if console.open(argv[1]):

            print_cmd('\nAimsun opening {0} ...\n'.format(argv[1]))

            #========================================================
            # Set up AIMSUN simulation
            #========================================================
            # Get the current Aimsun  model
            model = GKSystem.getSystem().getActiveModel()

            demand = load_demand_from_ang(model, 'congflow_demand')

            # Setup model
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

            #========================================================
            # generate the result using default values
            #========================================================
            if simulate_default is True:

                set_new_paras(model, experiment, default_paras)
                # test
                console.save("thread2_calib_I80_EB_default.ang")
                print_cmd('\nthread2_calib_I80_EB_default.ang saved.')

                default_result = simulate_rep_from_paras(model, experiment, default_paras,
                                                         simulator, avg_result, plugin,
                                                         valid_data, det_used_weight,
                                                         'default')

                save_solution_data(default_paras, default_result[1], logger_path, start_time, 'default')

                # --------------------------------------------------
                # For printing result
                # compute the objective function value if using the true parameters
                # true_obj_value = valid_result[0]
                default_obj_val = default_result[0]

                best_obj = deepcopy(default_obj_val)

            else:
                default_obj_val = (-1, -1)

            #========================================================
            # generate the result using user specified values
            #========================================================
            if simulate_user is True:

                set_new_paras(model, experiment, user_paras)
                # test
                console.save("thread2_calib_I80_EB_user.ang")
                print_cmd('\nthread2_calib_I80_EB_user.ang saved.')

                user_result = simulate_rep_from_paras(model, experiment, user_paras,
                                                         simulator, avg_result, plugin,
                                                         valid_data, det_used_weight,
                                                         'user')

                save_solution_data(user_paras, user_result[1], logger_path, start_time, 'user')

                # --------------------------------------------------
                # For printing result
                # compute the objective function value if using the true parameters
                # true_obj_value = valid_result[0]
                user_obj_val = user_result[0]

            else:
                user_obj_val = (-1, -1)

            #========================================================
            # Iterating parameters
            #========================================================
            opt_solution = []

            iter_counter = 1
            simulated_iter_counter = 1
            # g_max_iter + 1 is to make sure OptQuest stops first
            while iter_counter <= g_max_iter + 1 and solved_flag is False:

                run_time = datetime.now() - start_time
                expected_finish_time = datetime.now() + run_time*(g_max_iter + 3 - iter_counter)/simulated_iter_counter

                print_cmd('\n-------------T2-------Iteration {0}------------------------'.format(iter_counter))

                # print out run time information
                print_cmd('Started simulation at: {0}'.format(start_time.strftime("%Y-%m-%d %H:%M:%S")))
                print_cmd('Simulation run time  : {0}'.format(str(run_time)))
                print_cmd('Expected to finish at: {0}\n'.format(expected_finish_time.strftime("%Y-%m-%d %H:%M:%S")))

                #==================================================
                # Read parameters
                #==================================================
                paras = read_optquest_solution(sim_sol_file)

                if paras is None:
                    # solution not found; keep going
                    # read initial/new parameters from OptQuest
                    paras = read_optquest_paras(sim_para_file)

                    # if paras is still None, then probably a solution has been written
                    if paras is None:
                        if timeout_flag is False:
                            if timeout_counter >=3:
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
                    print_cmd('\n------------------Optimal solution----------------------'.format(iter_counter))
                    solved_flag = True


                #==================================================
                # Simulate and get the objective function value
                #==================================================
                # first check previous solutions
                tmp_obj_val = try_get_obj_from_previous_solutions(solution_list,paras)

                if tmp_obj_val is None:

                    # set the default parameters, and then overwrite. The reason is the new paras may not be
                    # complete, hence need to make sure the other paras are in default.
                    # the true parameters are the default.
                    print_cmd('Reset default paras:')
                    set_new_paras(model, experiment, default_paras)
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

                    # console.save("test_I80_EB_middle.ang")
                    # print_cmd('\ntest_I80_EB_middle.ang saved.')

                else:
                    # [RMS_speed, RMS_count]; just follow the tradition, though we are not using count
                    obj_value = tmp_obj_val
                    print_cmd('Got obj value {0} from previous solutions'.format(tmp_obj_val))

                #==================================================
                # send simulation value to OptQuest
                #==================================================
                write_simval(obj_fun[0]*obj_value[0] + obj_fun[1]*obj_value[1], sim_val_file)

                #==================================================
                # update the current best value
                #==================================================
                if obj_fun[0]*obj_value[0] + obj_fun[1]*obj_value[1] <= \
                    obj_fun[0]*best_obj[0] + obj_fun[1]*best_obj[1] :
                    # found better parameters
                    best_para = deepcopy(paras)
                    best_obj = deepcopy(obj_value)

                #==================================================
                # Log the result
                #==================================================
                log_opt_step(opt_solution, iter_counter, paras, obj_value)

                print_opt_steps(opt_solution, obj_fun, [('default', default_obj_val)])

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
                console.save("thread2_calib_I80_EB_optimal.ang")
                print_cmd('\nthread2_calib_I80_EB_optimal.ang saved.')

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
                    newseed_result = simulate_rep_from_paras(model, experiment, optimal_paras,
                                                             simulator, avg_result, plugin,
                                                             valid_data, det_used_weight,
                                                             'newseed')
                    save_solution_data(optimal_paras, newseed_result[1], logger_path, start_time, 'newseed')
                    newseed_obj_val = newseed_result[0]
                else:
                    newseed_obj_val = (-1, -1)

            end_time = datetime.now()


            # if timeout, then print out the optimal parameter so far
            if timeout_flag is True:
                set_random_seed(avg_result)
                newseed_result = simulate_rep_from_paras(model, experiment, best_para,
                                                             simulator, avg_result, plugin,
                                                             valid_data, det_used_weight,
                                                             'newseed')
                save_solution_data(best_para, newseed_result[1], logger_path, start_time, 'newseed')
                newseed_obj_val = newseed_result[0]

                print_cmd(  '\n\n==================== Calibration timeout ===================================')
                print_cmd(  '===========================================================================')

                print_results([ ('default', default_paras),
                                ('optimal', best_para),
                                ('newseed', best_para)],
                              [ ('default', default_obj_val),
                                ('optimal', best_obj),
                                ('newseed', newseed_obj_val)])

            if solved_flag is True:

                print_results([ ('optimal', optimal_paras),
                                ('default', default_paras)],
                              [ ('optimal',optimal_obj_val),
                                ('default', default_obj_val),
                                ('newseed', newseed_obj_val)])
            else:

                print_cmd(  '\n\n====================Calibration timeout===================================')
                print_cmd(  '===========================================================================')

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

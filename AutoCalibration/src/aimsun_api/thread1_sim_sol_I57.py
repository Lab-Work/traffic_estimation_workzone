__author__ = 'Carlos and Yanning'

# This code is used for automating AIMSUN and interacting with the optimization solver.
# All functions are in AIMSUNFUNCTIONS_V3


from AIMSUNFUNCTIONS_V3 import *
from copy import deepcopy

# ======================================================================================
# Configurations in this block
# ======================================================================================

# Something that does not have to be modified===============================
# simulation para file, and simulation obj value file path
sim_para_file = 'C:/thread1_sim_com/sim_paras.txt'
sim_val_file = 'C:/thread1_sim_com/sim_val.txt'
sim_sol_file = 'C:/thread1_sim_com/sim_sol.txt'

# this is the traffic state names that should be used in .ang
# demand name: traffic_state + '_demand'
# scenario name: traffic_state
# experiment name: traffic_state + '_exp'
traffic_state = 'congflow'
validFilePath = 'C:/Users/TrafficControl/Dropbox/DesktopServer/Workzone_autocalibration' + \
                '/I57_SB/validation_data_15mph/'  # in current folder
logger_path = 'C:/Users/TrafficControl/Dropbox/DesktopServer/Workzone_autocalibration' + \
              '/I57_SB/Logs/thread1'
# opt_step_file = 'c:/Users/TrafficControl/Google Drive/AIMSUN_carlos/Workzone_calibration/I80_EB/congflow/Logs/20150723_150924/20150723_150924_opt_log.csv'
opt_step_file = None


# ---------- Specify demand ----------
main_entrance_id = 51413
main_entrance_discount = 1.0


# ---------- Configure seeds ----------
# number of replications. We are not calibrating a work zone. We are REPLICATING a work zone.
# Hence, one replication is sufficient.
g_num_rep = 1
seed_list = [23798]
# seed_list = [3503, 23798, 28860, 12358, 6370, 14452, 21488, 7893, 586, 30447]

# make sure OptQuest ends first. This program ends once it found sim_sol.txt file
# this value is used for calculating the computation time
g_max_iter = 800


# ---------- Detector list for computing the error ----------
# EB3 was used for generating the main inflow, hence do not use as the validation data.
all_det_strs =   ['SB1','SB2','SB3','SB4','SB5','SB6','SB7', 'SB8', 'SB9']
# detector used for validation
det_used = all_det_strs
det_used_weight = {'SB1':1.0, 'SB2':1.0, 'SB3':1.0,
                   'SB4':1.0, 'SB5':1.0, 'SB6':1.0,
                   'SB7':1.0, 'SB8':1.0, 'SB9':1.0}

# Define the way to compute the objective function
# obj_fun[0]*RMS_speed + obj_fun[1]*RMS_count
obj_fun = [1, 0]


# ------------ For processing the validation data ------------
# the following two time strings are used to import the validation data in the selected time period
# They should be the exact start and end time of the simulated period (e.g. fist 15 cong states:  05/01/2015 15:30 ~ 05/01/2015 16:45)
start_time_str = '11/26/2014 15:30'
end_time_str = '11/26/2014 18:00'


# Put one line of description of the simulation
description = 'This is the calibration of the {0} state ib I57_SB. Simulate time from {1} to {2}\n'.format(traffic_state, start_time_str, end_time_str) + \
              'Calibrate with missing speed replaced by 15 mph. Demand is manually calibrated.'


# deprecated
# parsFilePath = 'c:/Users/TrafficControl/Google Drive/AIMSUN_carlos/Workzone_calibration/I57_SB/'+ folder_name +'/demand_data/paras.txt'
# flowsFilePath = 'c:/Users/TrafficControl/Google Drive/AIMSUN_carlos/Workzone_calibration/I57_SB/'+ folder_name +'/demand_data/flows.txt'
# turnsFilePath = 'c:/Users/TrafficControl/Google Drive/AIMSUN_carlos/Workzone_calibration/I57_SB/'+ folder_name +'/demand_data/turns.txt'
# presetParaFile = 'c:/Users/TrafficControl/Google Drive/AIMSUN_carlos/Workzone_calibration/I57_SB/'+ folder_name +'/preset_paras.txt'


# ======================================================================================
# The parameters for simulation
# ======================================================================================
# set up the simulation flag
simulate_default = False
simulate_newseed = False
simulate_user = False

#============================================================================
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
user_paras = {}
user_paras['car_speedAcceptance'] = [1.1, 0.1, 0.85, 1.3]  # calibrated value from freeflow
user_paras['truck_maxAccel'] = [1.68, 0.5, 0.6, 1.8]
user_paras['car_sensitivityFactor'] = [0.8, 0, 0.8, 0.8]
user_paras['truck_sensitivityFactor'] = [0.85, 0, 0.85, 0.85]
user_paras['car_reactionTime'] = [0.8, 1.2, 1.6, 1]
user_paras['truck_reactionTime'] = [0.8, 1.3, 1.7, 1]


user_paras['car_maxAccel'] = [3.0, 0.2, 2.6, 3.4]
user_paras['car_minDist'] = [4, 0.3, 0.5, 10]     # range in OptQuest [0.5, 10]
user_paras['truck_speedAcceptance'] = [0.87, 0.1, 0.85, 1.1]
user_paras['truck_minDist'] = [15, 0.5, 1.0, 15]

true_paras = None

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
    print_cmd('=============================== Thread 1 ===============================')
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
                console.save("thread1_calib_I57_SB_default.ang")
                print_cmd('\nthread1_calib_I57_SB_default.ang saved.')

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
                print_cmd('Reset default paras:')
                set_new_paras(model, experiment, default_paras)
                print_cmd('Overwrite user defined paras:')
                set_new_paras(model, experiment, user_paras)
                # test
                console.save("thread1_calib_I57_SB_user.ang")
                print_cmd('\nthread1_calib_I57_SB_user.ang saved.')

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
            # Simulate optimal
            #========================================================
            paras = read_optquest_solution(sim_sol_file)
            if paras is None:
                raise Exception('No optimal solution found')
            else:
                print_cmd('Reset default paras:')
                set_new_paras(model, experiment, default_paras)
                print_cmd('Overwrite with optimal paras:')
                set_new_paras(model, experiment, paras)
                # save the ang file
                console.save("thread1_calib_I57_SB_optimal.ang")
                print_cmd('\nthread1_calib_I57_SB_optimal.ang saved.')

                # re-simulate, in case no data available
                result = simulate_rep_from_paras(model, experiment, paras,
                                                         simulator, avg_result, plugin,
                                                         valid_data, det_used_weight,
                                                         'iteration Optimal')

                optimal_obj_val = result[0]
                optimal_data = result[1]
                save_solution_data(paras, optimal_data, logger_path, start_time, 'optimal')

                #========================================================
                # see how the seed affects the optimal solution
                if simulate_newseed is True:
                    set_random_seed(avg_result)
                    newseed_result = simulate_rep_from_paras(model, experiment, paras,
                                                             simulator, avg_result, plugin,
                                                             valid_data, det_used_weight,
                                                             'newseed')
                    save_solution_data(paras, newseed_result[1], logger_path, start_time, 'newseed')
                    newseed_obj_val = newseed_result[0]
                else:
                    newseed_obj_val = (-1, -1)


            # save logger to file
            stop_opt_log()


            print_cmd( 'Finished simulating optimal parameters' )
            print_cmd(  '\nAimsun is now closing...\n')

            # stop logger and save into file
            stop_cmd_logger()

            console.close()

        else:
            console.getLog().addError("Could not open")


if __name__ == "__main__":
    sys.exit(main(sys.argv))

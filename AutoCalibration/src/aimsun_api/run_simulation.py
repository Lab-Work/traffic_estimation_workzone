from AIMSUN_API import *
from copy import deepcopy

"""
This script is used to start AIMSUN simulation for autocalibration.
It communicates with optimizer NOMAD via txt file.
    - Once received parameters, it simulates the parameters in AIMSUN, and write the objective value in file.
    - NOMAD proposes new parameters to simulate based on the last objective value.

Usage:
    - cd to the current folder containing this function
    - run: ./run_simulation.py [workzone] [thread]
        [workzone]: 'I80' or 'I57' or 'I00'
        [thread]: 'thread1' or 'thread2'; Current version of AIMSUN only allows running two instances simultaneously.
    - start NOMAD optimizer.
"""

__author__ = 'Yanning Li'


def main(argv):

    # =================================================
    # Read and parse the configuration file for this work zone and thread
    workzone = argv[1]
    thread = argv[2]
    config = parse_config(workzone, thread)

    # =================================================
    # a timeout flag. If not paras or solution paras are obtained in 60 s. Then stop AIMSUN.
    timeout_flag = False
    timeout_counter = 0

    solved_flag = False

    # keep track of the current best parameter and its associated objective
    default_paras = load_paras('../data/I80/default_paras.txt')
    best_para = deepcopy(default_paras)
    best_obj = 0

    # keep track of start time:
    start_time = datetime.now()

    # Start the cmd output logger
    start_cmd_logger(config['logger_path'], start_time)
    start_opt_logger(config['logger_path'], start_time)

    # print python version
    # print '\nPython interpreter version: {0} \n'.format(sys.version)

    print_cmd('\n\n========================================================================')
    print_cmd('=========================AIMSUN Autocalibration=========================')
    print_cmd('=============================== Thread 1 ===============================')
    print_cmd('Auto calibration started at {0}\n'.format(datetime.now().strftime("%Y-%m-%d %H:%M:%S")))

    # write notes and parameters to cmd_log
    cmd_logger_header('AIMSUN calibration for {0} {1}'.format(workzone, thread),
                      config['g_max_iter'], config['g_num_rep'], config['seed_list'],
                      config['obj_fun'], config['det_used'],
                      config['main_entrance_id'], config['main_entrance_discount'],
                      [('default', default_paras)])

    #========================================================
    # validation data for computing the objective
    # ----read from file
    valid_data = read_validation_data(config['start_time_str'], config['end_time_str'],
                                      config['det_used'], config['validFilePath'])

    # save the complete validation data used for visualization
    all_valid_data = read_validation_data(config['start_time_str'], config['end_time_str'],
                                          config['all_det_strs'], config['validFilePath'])
    save_solution_data(None, all_valid_data, config['logger_path'], start_time, 'valid')


    #========================================================
    # load previous solution
    if config['opt_step_file'] is not None:
        solution_list = load_previous_opt_solutions(config['opt_step_file'])
        print_cmd('\nLoaded previous OptQuest solutions:\n---- {0}'.format(config['opt_step_file']))

        # get the current best parameters
        previous_obj_tuple = np.array( solution_list[1] )
        previous_obj_val = config['obj_fun'][0]*previous_obj_tuple[:,0] + config['obj_fun'][1]*previous_obj_tuple[:,1]

        idx = np.argmin( previous_obj_val )

        best_para = deepcopy(solution_list[0][idx])
        best_obj = deepcopy(solution_list[1][idx])

        print_cmd('\nBest obj from loaded parameters: {0}\n'.format( config['obj_fun'][0]*best_obj[0] +
                                                                   config['obj_fun'][1]*best_obj[1] ))

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
            scenario = setup_scenario(model, config['traffic_state'], demand)

            # create experiment
            experiment = setup_experiment(model, config['traffic_state'] + '_exp', scenario)

            # create replications
            avg_result = setup_replication(model, experiment, config['g_num_rep'], config['seed_list'])

            # create simulator
            simulator = create_simulator(model)
            # plugin is a module which can compute the average for the GKExperimentResult object
            plugin = GKSystem.getSystem().getPlugin( "GGetram" )

            #========================================================
            # generate the result using default values
            #========================================================
            if config['simulate_default'] is True:

                set_new_paras(model, experiment, default_paras)
                # test
                console.save("thread1_calib_I80_EB_default.ang")
                print_cmd('\nthread1_calib_I80_EB_default.ang saved.')

                default_result = simulate_rep_from_paras(model, experiment, default_paras,
                                                         simulator, avg_result, plugin,
                                                         valid_data, config['det_used_weight'],
                                                         'default')

                save_solution_data(default_paras, default_result[1], config['logger_path'], start_time, 'default')

                # --------------------------------------------------
                # For printing result
                # compute the objective function value if using the true parameters
                # true_obj_value = valid_result[0]
                default_obj_val = default_result[0]

                best_obj = deepcopy(default_obj_val)

            else:
                default_obj_val = (-1, -1)

            #========================================================
            # Iterating parameters
            #========================================================
            opt_solution = []

            iter_counter = 1
            simulated_iter_counter = 1
            # g_max_iter + 1 is to make sure OptQuest stops first
            while iter_counter <= config['g_max_iter'] + 1 and solved_flag is False:

                run_time = datetime.now() - start_time
                expected_finish_time = datetime.now() + run_time*(config['g_max_iter'] + 3 - iter_counter)/simulated_iter_counter

                print_cmd('\n-----------T1-------Iteration {0}------------------------'.format(iter_counter))

                # print out run time information
                print_cmd('Started simulation at: {0}'.format(start_time.strftime("%Y-%m-%d %H:%M:%S")))
                print_cmd('Simulation run time  : {0}'.format(str(run_time)))
                print_cmd('Expected to finish at: {0}\n'.format(expected_finish_time.strftime("%Y-%m-%d %H:%M:%S")))

                #==================================================
                # Read parameters
                #==================================================
                paras = read_optquest_solution(config['sim_sol_file'])

                if paras is None:
                    # solution not found; keep going
                    # read initial/new parameters from OptQuest
                    paras = read_optquest_paras(config['sim_para_file'])

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
                                                         valid_data, config['det_used_weight'],
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
                write_simval(config['obj_fun'][0]*obj_value[0] + config['obj_fun'][1]*obj_value[1], config['sim_val_file'])

                #==================================================
                # update the current best value
                #==================================================
                if config['obj_fun'][0]*obj_value[0] + config['obj_fun'][1]*obj_value[1] <= \
                    config['obj_fun'][0]*best_obj[0] + config['obj_fun'][1]*best_obj[1] :
                    # found better parameters
                    best_para = deepcopy(paras)
                    best_obj = deepcopy(obj_value)

                #==================================================
                # Log the result
                #==================================================
                log_opt_step(opt_solution, iter_counter, paras, obj_value)

                print_opt_steps(opt_solution, config['obj_fun'], [('default', default_obj_val)])

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
                console.save("thread1_calib_I80_EB_optimal.ang")
                print_cmd('\nthread1_calib_I80_EB_optimal.ang saved.')

                # re-simulate, in case no data available
                result = simulate_rep_from_paras(model, experiment, optimal_paras,
                                                         simulator, avg_result, plugin,
                                                         valid_data, config['det_used_weight'],
                                                         'iteration Optimal')

                optimal_obj_val = result[0]
                optimal_data = result[1]
                save_solution_data(optimal_paras, optimal_data, config['logger_path'], start_time, 'optimal')

                #========================================================
                # see how the seed affects the optimal solution
                if config['simulate_newseed'] is True:
                    set_random_seed(avg_result)
                    newseed_result = simulate_rep_from_paras(model, experiment, optimal_paras,
                                                             simulator, avg_result, plugin,
                                                             valid_data, config['det_used_weight'],
                                                             'newseed')
                    save_solution_data(optimal_paras, newseed_result[1], config['logger_path'], start_time, 'newseed')
                    newseed_obj_val = newseed_result[0]
                else:
                    newseed_obj_val = (-1, -1)

            end_time = datetime.now()


            # if timeout, then print out the optimal parameter so far
            if timeout_flag is True:
                set_random_seed(avg_result)
                newseed_result = simulate_rep_from_paras(model, experiment, best_para,
                                                             simulator, avg_result, plugin,
                                                             valid_data, config['det_used_weight'],
                                                             'newseed')
                save_solution_data(best_para, newseed_result[1], config['logger_path'], start_time, 'newseed')
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


def parse_config(workzone, thread):

    config = {}

    config_f = '../{0}_configuration.txt'.format(workzone)

    # parse the file
    with open(config_f) as f:
        for line in f:
            if len(line) == 0 or line[0] == '#':
                continue

            items = line.strip().split(':')

            # three types of values: string, number, and list
            if items[0] == 'com_dir' or items[0] == 'traffic_state' or  items[0] == 'validFilePath' \
                or items[0] == 'logger_path' or items[0] == 'start_time_str' \
                or items[0] == 'end_time_str':
                # for all string valued entries
                config[items[0]] = items[1]

            elif items[0] == 'main_entrance_id' or items[0] == 'g_max_iter' or items[0] == 'g_num_rep':
                # for all int valued entries
                config[items[0]] = int(items[1])

            elif items[0] == 'main_entrance_discount':
                # for all float valued centries
                config[items[0]] = float(items[1])

            elif items[0] == 'seed_list':
                # list of int
                config[items[0]] = [int(i) for i in items[1].split(',')]

            elif items[0] == 'obj_fun':
                # list of float
                config[items[0]] = [float(i) for i in items[1].split(',')]

            elif items[0] == 'all_det_strs' or items[0] == 'det_used':
                # list of strings
                config[items[0]] = items[1].split(',')

            elif items[0] == 'det_used_weight':
                # a dict
                config[items[0]] = {}
                for entry in items[1].split((';')):
                    config[items[0]][entry.split(',')[0]] = float(entry.split(',')[1])

            elif items[0] == 'simulate_default' or items[0] == 'simulate_user' or items[0] == 'simulate_newseed' \
                    or items[0] == 'opt_step_file':
                if items[1] == 'True':
                    config[items[0]] = True
                elif items[1] == 'False':
                    config[items[0]] = False
                elif items[1] == 'None':
                    config[items[0]] = None
                else:
                    raise Exception('simulate_xxx can only take these values: True, False')

            else:
                raise Exception('Unrecognized configuration entry')

    # specify the folder depending on the thread specification
    config['sim_para_file'] = config['com_dir'] + thread + 'sim_paras.txt'
    config['sim_val_file'] = config['com_dir'] + thread + 'sim_val.txt'
    config['sim_sol_file'] = config['com_dir'] + thread + 'sim_sol.txt'
    config['validFilePath'] += '{0}/validation_data/'.format(workzone)
    config['logger_path'] += '{0}/Logs/{1}/'.format(workzone, thread)

    return config


def load_paras(parafile):
    """
    A universal paras reader assuming the file exist.
    :param parafile: the file to read
    :return: a dict, key-value
    """
    paras = OrderedDict()

    f = open(parafile, 'r')

    for line in f:
        line = line.strip()
        items = line.split(',')

        # first item is the key
        paras[items[0]] = []

        for i in range(1, len(items)):
            paras[items[0]].append(float(items[i]))

    f.close()

    return paras


if __name__ == "__main__":
    sys.exit(main(sys.argv))


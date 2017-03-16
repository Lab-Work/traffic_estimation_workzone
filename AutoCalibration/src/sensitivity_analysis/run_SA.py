__author__ = 'Yanning Li'
"""
This script conducts the sensitivity analysis of twelve parameters and narrow down to fewer parameters for calibration.

Background:
Sensitivity analysis is a method for determining the influence of each factor to the system response. Ideally, an ANOVA
should be performed to obtain statistic significance of each factor and the interaction between factors. However, ANOVA
requires a large number of simulations. E.g., if 12 factors, 4 levels of value each, then in total 4^12 = 16,777,216.
Assuming each iteration of the simulation takes 10 s, it will take 5.3 years to finish the simulation. Realistically,
the maximum number of factors should be 8, which takes roughly 10 days to finish ANOVA test.

Alternatively, we can roughly infer the sensitivity by investigating one factor at a time while keeping other factors
fixed as default. This relies on an assumption: the factors are independent. Though this assumption is not completely
justified, it provides a faster way to infer the sensitivity of factors.

Setup of the experiment.
1. Twelve parameters are investigated by X-Y plot (change one parameter at a time), including the follow for two types
of vehicles:
- speed acceptance. (mean, std, min, max) NOTE: the desired speed must be set unrealistically high in AIMSUN ang.
        The reason is the vehicle speed is determined by min( desired_speed, speed_limit*speed_acceptance ).
        Set the desired_speed to be high for the speed_acceptance to the active one.
- maximum acceleration. (mean, std, min, max), This parameter for Gipps car-following model.
- reaction time. scalar. This is the time it takes a driver to react to speed changes in the preceding vehicle.
        It is used in the car-following model
- minimum distance. (mean, std, min, max) This is the distance, in metres, that a vehicle keeps between itself and the
        preceding vehicle when stopped. Smaller distance => higher capacity.
- sensitivity factor. scalar. In the deceleration component of the car-following model, the follower makes an estimation
        of the deceleration of the leader using the sensitivity factor.
        If lower than 1, the vehicle underestimate the deceleration of the leader => more aggressive.
        If higher than 1, then less aggressive.
- coefficient of variance. std/mean. We assume a single parameter denotes the std/mean for three parameters with scale:
            Car
            + speed_acceptance, std = mean*cov/3.3;
            + accel, std = mean*cov/4.5;
            + minDist, std = mean*cov
            Car COV , default, 0.3, range [0.1, 0.5]
            Truck
            + speed_acceptance, std = mean*cov/3.5;
            + accel, std = mean*cov/0.66;
            + minDist, std = mean*cov
            Truck COV , default, 0.33, range [0.1, 0.5]

2. The initial parameter values are default.

3. Make sure the maximum desired speed for car and truck are unrealistically high. See speed accpetance factor for reasons.

4. Make sure the vehicle generation is set as const, hence generating the flow exact as specified.

5. Fix all other parameters as default, and change one parameter to see the response.

6. Make sure AIMSUN has ten replications replication and set new seeds in every run.

"""


# This code does a sensitivity analysis:
# Given the optimal parameters, and the bounds of each parameter,
# this code then change one parameters at each round, while keeping the other parameters
# at the optimal values.

import sys
import os
from os.path import exists
import numpy as np
import matplotlib.pyplot as plt
import csv
import time
from datetime import datetime
import copy


# here is the address of some files to communicate with AIMSUN
sim_para_file = 'C:/sim_com/sim_paras.txt'
sim_val_file = 'C:/sim_com/sim_val.txt'
sim_sol_file = 'C:/sim_com/sim_sol.txt' # write this file to terminate AIMSUN

folder_name = 'congflow'

logger_path = 'C:/Users/TrafficControl/Google Drive/AIMSUN_carlos/Workzone_calibration/I00_EB/'+ folder_name +'/Logs/'

# set the initial value as the default.
default_paras = {}
default_paras['car_speedAcceptance'] = [1.1, 0.1, 0.85, 1.3] # range: [0.85, 1.3]
default_paras['car_maxAccel'] = [3, 0.2, 2.6, 3.4]      # range: [2.6, 3.4]
default_paras['car_reactionTime'] = [0.8, 1.2, 1.6, 1]  # range: [0.5, 2]
default_paras['car_minDist'] = [1, 0.3, 0.5, 1.5]       # range: [0.5, 1.5]
default_paras['car_sensitivityFactor'] = [1, 0, 1, 1]   # range: [0.5, 1.5]
default_paras['car_cov'] = 0.3                          # range: [0.1, 0.5]

default_paras['truck_speedAcceptance'] = [1.05, 0.1, 0.85, 1.1] # range: [0.85, 1.1]
default_paras['truck_maxAccel'] = [1, 0.5, 0.6, 1.8]        # range: [0.6, 1.8]
default_paras['truck_reactionTime'] = [0.8, 1.3, 1.7, 1]    # range: [0.5, 2]
default_paras['truck_minDist'] = [1.5, 0.5, 1.0, 2.5]       # range: [1.0, 2.5]
default_paras['truck_sensitivityFactor'] = [1, 0, 1, 1]     # range: [0.5, 1.5]
default_paras['truck_cov'] = 0.33                           # range: [0.1, 0.5]

# The fine grid
bounds = {}
# min, max, step
bounds['car_speedAcceptance'] = [0.85, 1.3, 0.01]   # 46
bounds['car_maxAccel'] = [2.6, 3.4, 0.01]           # 81
bounds['car_reactionTime'] = [0.4, 2.0, 0.2]        # 3
bounds['car_minDist'] = [0.5, 1.5, 0.01]             # 101
bounds['car_sensitivityFactor'] = [0.5, 1.5, 0.01]  # 101
bounds['car_cov'] = [0.1, 0.5, 0.01]

bounds['truck_speedAcceptance'] = [0.85, 1.1, 0.01]     # 26
bounds['truck_maxAccel'] = [0.6, 1.8, 0.01]             # 121
bounds['truck_reactionTime'] = [0.4, 2.0, 0.2]          # 3
bounds['truck_minDist'] = [1, 2.5, 0.01]                # 151
bounds['truck_sensitivityFactor'] = [0.5, 1.5, 0.01]    # 101
bounds['truck_cov'] = [0.15, 0.15, 0.01]


# The coarse grid just to debug code
# bounds = {}
# # min, max, step
# bounds['car_speedAcceptance'] = [0.85, 1.3, 0.45]   # 46
# bounds['car_maxAccel'] = [2.6, 2.6, 1]           # 81
# bounds['car_reactionTime'] = [0.8, 0.8, 1]        # 3
# bounds['car_minDist'] = [0.5, 0.5, 1]             # 101
# bounds['car_sensitivityFactor'] = [0.5, 0.5, 1]  # 101
# bounds['car_cov'] = [0.1, 0.1, 1]
#
# bounds['truck_speedAcceptance'] = [0.85, 0.85, 1]     # 26
# bounds['truck_maxAccel'] = [0.6, 1.8, 1.2]             # 121
# bounds['truck_reactionTime'] = [0.8, 0.8, 1]          # 3
# bounds['truck_minDist'] = [1, 1, 1]                # 151
# bounds['truck_sensitivityFactor'] = [0.5, 0.5, 1]    # 101
# bounds['truck_cov'] = [0.1, 0.1, 1]





def main(argv):

    start_time = datetime.now()
    start_cmd_logger(logger_path, start_time)

    # the following variables are used to estimate run time
    iter_counter = 1
    timeout_flag = False

    paras_to_analyze = ['car_reactionTime', 'truck_reactionTime']


    # Compute total number of different parameters. This value is used to compute the simulation time.
    g_total_iter = 0
    for key in paras_to_analyze:
        g_total_iter += (bounds[key][1] - bounds[key][0])/bounds[key][2] + 1
    g_total_iter = int(g_total_iter)

    # start logger
    cmd_logger_header(paras_to_analyze, start_time)

    analysis_logger = {}

    print_cmd('\n=========Total {0} sets of paras =============\n'.format(g_total_iter))

    for para_name in paras_to_analyze:

        print_cmd('\n===================== Analyzing {0}  =====================\n'.format(para_name))

        start_value = bounds[para_name][0]
        end_value = bounds[para_name][1]
        step_size = bounds[para_name][2]

        # the first list is the new value,
        # the second list is the objective function value
        analysis_logger[para_name] = [[],[]]

        try_value = start_value
        # +0.005 is to make sure the end_value is included. Otherwise, there may be numerical error
        while try_value <= end_value+0.005:

            run_time = datetime.now() - start_time
            expected_finish_time = datetime.now() + \
                                   run_time*(  g_total_iter - iter_counter)/iter_counter

            print_cmd('\n--------------------Set {0} / {1}------------------------'.format(iter_counter, g_total_iter))

            # print out run time information
            print_cmd('Started simulation at: {0}'.format(start_time.strftime("%Y-%m-%d %H:%M:%S")))
            print_cmd('Simulation run time  : {0}'.format(str(run_time)))
            print_cmd('Expected to finish at: {0}\n'.format(expected_finish_time.strftime("%Y-%m-%d %H:%M:%S")))


            try_value = np.round(try_value, decimals=2)
            paras = generate_new_para(para_name, try_value)

            print_cmd('\n------ Analyzing {0}: {1:.2f}\n'.format(para_name, try_value))

            # for each parameters, simulate 10 iterations
            for i in range(0,10):

                # write paras for AIMSUN to simulate
                write_paras(paras, sim_para_file)

                # wait for teh objective function value
                # we wait in the function
                obj_val = read_obj_value(sim_val_file)

                # If obj_val is none, then it means the AIMSUN simulator stopped,
                # Break and save the result
                if obj_val is None:
                    timeout_flag = True
                    break

                analysis_logger[para_name][0].append(try_value)
                analysis_logger[para_name][1].append(obj_val)

            if timeout_flag is True:
                break

            try_value += step_size

            iter_counter += 1

        # save the result into logger. It will be overwritten once new parameter is analyzed
        write_analysis_steps(analysis_logger, logger_path, start_time)

        if timeout_flag is True:
            break

    # all paras done, then write the sim_sol to stop AIMSUN
    write_paras(default_paras, sim_sol_file)

    stop_cmd_logger()



# this function replace the value entry in the start paras by the new set value
def generate_new_para(para_name, new_value):

    new_para = copy.deepcopy(default_paras)

    # remove the keys for the cov since they are computed and set in other parameters.
    if 'car_cov' in new_para: del new_para['car_cov']
    if 'truck_cov' in new_para: del new_para['truck_cov']

    if para_name == 'car_speedAcceptance' \
            or para_name == 'car_maxAccel' \
            or para_name == 'car_reactionTime'\
            or para_name == 'car_minDist'\
            or para_name == 'truck_speedAcceptance'\
            or para_name == 'truck_maxAccel'\
            or para_name == 'truck_reactionTime'\
            or para_name == 'truck_minDist':
        new_para[para_name][0] = new_value

    elif para_name == 'car_sensitivityFactor'\
            or para_name == 'truck_sensitivityFactor':
        new_para[para_name][0] = new_value
        new_para[para_name][2] = new_value
        new_para[para_name][3] = new_value

    elif para_name == 'car_cov':
        # change the parameter to the maxAccel, minDist, and speedAcceptance
        new_para['car_speedAcceptance'][1] = np.round(new_para['car_speedAcceptance'][0]*new_value/3.3, decimals=2)
        new_para['car_maxAccel'][1] = np.round(new_para['car_maxAccel'][0]*new_value/4.5, decimals=2)
        new_para['car_minDist'][1] = np.round(new_para['car_minDist'][0]*new_value, decimals=2)

    elif para_name == 'truck_cov':
        # change the parameter to the maxAccel, minDist, and speedAcceptance
        new_para['truck_speedAcceptance'][1] = np.round(new_para['truck_speedAcceptance'][0]*new_value/3.5, decimals=2)
        new_para['truck_maxAccel'][1] = np.round(new_para['truck_maxAccel'][0]*new_value/0.666, decimals=2)
        new_para['truck_minDist'][1] = np.round(new_para['truck_minDist'][0]*new_value, decimals=2)

    else:
        print_cmd('ERROR: unrecognized para name: {0}\n'.format(para_name))

    return new_para



# write the new parameters
def write_paras(paras, para_file):

    f = open(para_file, 'w')

    for key in paras.keys():

        if key != 'car_cov' and key!= 'truck_cov':
            # Those two are already converted to the std in other parameters
            f.write(key)

            for val in paras[key]:
                f.write( ',' + str(val))

            f.write('\n')

    f.close()

    print_paras(paras)




# read the new objective value
def read_obj_value(obj_val_file):

    wait_time = 0
    timeout = 0
    print_cmd('waiting for objective value...\n')
    while not exists(obj_val_file):
        time.sleep(0.1)   # sleep 1 second
        wait_time += 1
        timeout += 1
        if wait_time >= 600: # print one every 1 min
            print_cmd('waiting for objective value...\n')
            wait_time = 0

        if timeout >= 3000: # waited 5 min
            print_cmd('Timeout Error.\n')
            return None

    if exists(obj_val_file):
        f = open(obj_val_file, 'r')

        line = f.readline()
        tup = line.split(',')
        obj_val = float(tup[0])

        f.close()

        os.remove(obj_val_file)

        print_cmd('Have read objective value: {0}\n'.format(obj_val))

        return obj_val




# print the paras just set
def print_paras(paras):

    print_cmd('Set new paras:')
    for key in paras.keys():
        if isinstance(paras[key], list):
            print_cmd('---- {0}: {1}'.format(key, ', '.join( '{0:.2f}'.format(val) for val in  paras[key]) ) )
        else:
            print_cmd('---- {0}: {1}'.format(key, paras[key] ) )


# start logger for the cmd output
def start_cmd_logger(path_str, start_time):

    global cmd_logger

    file_name = start_time.strftime("%Y%m%d_%H%M%S") + '_SA_cmd_log'
    cmd_logger = open(path_str + file_name + '.txt', 'wb')

    # return cmd_logger


# command logger header
# para_list =[ (name, paras) ]
def cmd_logger_header(paras_list, start_time):

    print_cmd('\n================================================================')
    print_cmd('=======================Sensitivity Analysis=====================')
    print_cmd('Start time: {0}'.format(start_time.strftime("%Y-%m-%d %H:%M:%S")))

    print_cmd('\nAnalyzing parameters:')
    for para_name in paras_list:
        print_cmd('-- {0}: {1}'.format(para_name, bounds[para_name]))




# print out on the cmd, and save the cmd output to a file
# If log_file is None, will only print
def print_cmd(line_str):

    if cmd_logger is None:
        print line_str
    else:
        # save every time
        # Hence even if the file is not correctly closed, those lines are still saved.
        cmd_logger.write(line_str + '\n')
        print line_str


# stop logger for the cmd output
def stop_cmd_logger():

    cmd_logger.close()



# save the analysis logger
# analysis_logger[para_name] = [[para_val], [obj_val]]
def write_analysis_steps(analysis_logger, path_name, start_time):

    file_name = start_time.strftime("%Y%m%d_%H%M%S") + '_SA_step_log'
    f = open(path_name + file_name + '.csv', 'wb')
    writer = csv.writer(f)

    for key in analysis_logger.keys():

        l_val = []
        l_val.append(key + '_val')
        for val in analysis_logger[key][0]:
            l_val.append(str(val))

        writer.writerow(l_val)

        l_obj = []
        l_obj.append(key + '_obj')
        for obj in analysis_logger[key][1]:
            l_obj.append(str(obj))

        writer.writerow(l_obj)

    f.close()




if __name__ == "__main__":
    sys.exit(main(sys.argv))








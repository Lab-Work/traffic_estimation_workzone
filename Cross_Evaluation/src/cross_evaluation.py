import bisect
import os
import random
import sys
import time
from collections import OrderedDict
from copy import deepcopy
from os.path import exists
import matplotlib.colors as colors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.offsetbox import AnchoredText
from scipy import optimize

from EnKF_AN import EnKF_AN
from FISB_workzone import FISB_workzone
from Interpolation import Interpolation
from Virtual_Sensors import Virtual_Sensors

np.seterr(invalid='ignore')

__author__ = 'Yanning Li'
"""
This class is used to:
    - Run estimators using a variety of sensor network configurations and algorithms
    - Calibrate the fundamental diagram
    - Visualize the true states
    - Visualize the evaluation results
    - Compute the estimation errors
    - Compare and visualize the estimation errors
"""


class CrossEval:
    """
    This class defines a cross evaluation environment which builds different configurations of the sensor network,
    exacts the corresponding virtual sensor data, and runs different algorithms.
    The output is the quantitative measures of the performance, as well as the visualization.
    """

    # ================================================================================ #
    # Public functions to be called.
    # ================================================================================ #
    def __init__(self, workzone='I57', log_dir=None, data_dir=None, grid_res=(5, 200), replications='all'):
        """
        The constructor of the cross evaluation class
        :param workzone: the work zone string
        :param log_dir: the top directory of the logs. Subfolders Virtual_sensor_data/ and Result/ will be created under
                    this folder
        :param data_dir: the top level directory of the data, contains /Trajectory_data/, /Virtual_sensor_data/,/Result/
        :param grid_res: tuple in (second, meter), the resolution of the discretization grid. It discretize to run the
                    estimation algorithms and assess its performance.
        :param replications: the list of replications to select from.
        :return:
        """

        # separately work on two workzones
        self.workzone = workzone
        self.grid_res = grid_res
        # load the topology of the work zone
        # self.sections and self.fwy_sec_order are used to convert the relative location on sections to the
        # absolute location
        self.time_grid, self.space_grid, self.workzone_topo, self.__aimsun_start_dur_step = \
            self.__load_workzone(log_dir + self.workzone + '_topology.txt')

        # the replications to estimate on
        if replications == 'all':
            self.active_replications = self.workzone_topo['replications']
        else:
            # check if there is any illegal replication ids
            for rep in replications:
                if rep not in self.workzone_topo['replications']:
                    raise Exception('Error:unknown replications {0}.'.format(rep))
            self.active_replications = replications

        # ===========================================================================================
        # set up the directories for the files and data
        # this file contains a list of the virtual sensor data previously generated in data_dir
        self.__log_dir = log_dir
        self.__sensor_data_log = OrderedDict()
        self.__sensor_data_to_fetch = OrderedDict()
        self.__result_log = OrderedDict()
        self.__sensor_data_dir = OrderedDict()
        self.__est_result_dir = OrderedDict()
        self.__truestate_result_dir = OrderedDict()
        for rep in self.active_replications:

            if sys.platform == 'darwin' or sys.platform == 'linux2' or sys.platform == 'linux':
                # Mac directory or linux directory
                self.__sensor_data_log[rep] = log_dir + 'Virtual_sensor_data/{0}_rep{1}_generated.txt'.format(
                    self.workzone, rep)
                self.__sensor_data_to_fetch[rep] = log_dir + 'Virtual_sensor_data/{0}_rep{1}_to_generate.txt'.format(
                    self.workzone, rep)
                self.__result_log[rep] = log_dir + 'Estimation_results/{0}_rep{1}_res{2}s{3}m_generated.txt'.format(
                    self.workzone, rep,
                    self.grid_res[0],
                    self.grid_res[1])

                self.__sensor_data_dir[rep] = data_dir + 'Virtual_sensor_data/{0}/rep{1}/'.format(self.workzone, rep)
                self.__truestate_result_dir[rep] = data_dir + 'Estimation_results/{0}/rep{1}/truestate/'.format(
                    self.workzone, rep)
                self.__est_result_dir[rep] = data_dir + 'Estimation_results/{0}/rep{1}/{2}s{3}m/'.format(self.workzone,
                                                                                                         rep,
                                                                                                         self.grid_res[
                                                                                                             0],
                                                                                                         self.grid_res[
                                                                                                             1])

            elif sys.platform == 'win32':
                self.__sensor_data_log[rep] = log_dir + 'Virtual_sensor_data\\{0}_rep{1}_generated.txt'.format(
                    self.workzone, rep)
                self.__sensor_data_to_fetch[rep] = log_dir + 'Virtual_sensor_data\\{0}_rep{1}_to_generate.txt'.format(
                    self.workzone, rep)
                self.__result_log[rep] = log_dir + 'Estimation_results\\{0}_rep{1}_res{2}s{3}m_generated.txt'.format(
                    self.workzone, rep,
                    self.grid_res[0],
                    self.grid_res[1])

                self.__sensor_data_dir[rep] = data_dir + 'Virtual_sensor_data\\{0}\\rep{1}\\'.format(self.workzone, rep)
                self.__truestate_result_dir[rep] = data_dir + 'Estimation_results\\{0}\\rep{1}\\truestate\\'.format(
                    self.workzone, rep)
                self.__est_result_dir[rep] = data_dir + 'Estimation_results\\{0}\\rep{1}\\{2}s{3}m\\'.format(
                    self.workzone, rep,
                    self.grid_res[0],
                    self.grid_res[1])

            else:
                raise Exception('Error: run "import sys; sys.platform" to check the platform. ')

            # make directories
            if not exists(self.__sensor_data_dir[rep]):
                os.makedirs(self.__sensor_data_dir[rep])
            if not exists(self.__truestate_result_dir[rep]):
                os.makedirs(self.__truestate_result_dir[rep])
            if not exists(self.__est_result_dir[rep]):
                os.makedirs(self.__est_result_dir[rep])

        # ===========================================================================================
        # keep a list of all configurations that need to run or evaluate cross.
        # [config_name]['sensors'][sensor_id][sensor_paras]
        #              ['algorithms'][algorithm_paras]
        self.all_config = OrderedDict()

        # ===========================================================================================
        # create a virtual sensor class which will be used to fetch the sensor data and generate the true states
        self.__data_generator = Virtual_Sensors(self.workzone, log_dir, data_dir,
                                                self.workzone_topo['sections'], self.workzone_topo['fwy_sec_order'],
                                                self.active_replications,
                                                self.__aimsun_start_dur_step,
                                                [self.space_grid[0], self.space_grid[-1]])

    def __load_workzone(self, topo_file):
        """
        This function loads the topology of the network
        :param topo_file: the file path
        :return: t_grid, x_grid, topology of the network including the location of ramps, freeway sec order and
            replications
        """
        start_dur_step = None
        t_grid = None
        x_grid = None

        workzone_topo = {'sections': OrderedDict(),
                         'loc_onramp': None,
                         'loc_offramp': None,
                         'fwy_sec_order': None,
                         'replications': None}

        f_topo = open(topo_file, 'r')

        for line in f_topo:

            if line[0] == '#':
                continue

            if 'sec_id' in line:
                items = line.strip().split(';')
                sec_id = int(items[0].split(':')[1])
                workzone_topo['sections'][sec_id] = {}

                # add length
                workzone_topo['sections'][sec_id][items[1].split(':')[0]] = float(items[1].split(':')[1])

                if 'upstream' in line:
                    workzone_topo['sections'][sec_id]['connections'] = {}
                    for i in range(2, len(items)):
                        entry = items[i].split(':')
                        if entry[0] == 'upstream':
                            workzone_topo['sections'][sec_id]['connections']['upstream'] = [int(j) for j in
                                                                                            entry[1].split(',')]
                        elif entry[0] == 'downstream':
                            workzone_topo['sections'][sec_id]['connections']['downstream'] = [int(j) for j in
                                                                                              entry[1].split(',')]

            elif 'fwy_sec_order' in line:
                items = line.strip().split(':')
                workzone_topo['fwy_sec_order'] = [int(i) for i in items[1].split(',')]

            elif 'loc_start_end' in line:
                items = line.strip().split(':')
                start = float(items[1].split(',')[0])
                end = float(items[1].split(',')[1])

                # get x_grid
                x_grid = np.arange(end, start, -self.grid_res[1])
                if x_grid[-1] != start:
                    x_grid = np.concatenate([x_grid, np.array([start])])
                x_grid = x_grid[::-1]

                # round every digit to two decimals
                for i in range(0, len(x_grid)):
                    x_grid[i] = round(x_grid[i], 2)

                x_grid = x_grid.tolist()

            elif 'time_start_dur_step' in line:
                start_dur_step = [float(i) for i in line.split(':')[1].split(',')]
                t_grid = np.arange(0, start_dur_step[1], self.grid_res[0])
                if t_grid[-1] != start_dur_step[1]:
                    t_grid = np.concatenate([t_grid, np.array([start_dur_step[1]])])

                # round every digit to two decimals
                for i in range(0, len(t_grid)):
                    t_grid[i] = round(t_grid[i], 2)

                t_grid = t_grid.tolist()

            elif 'loc_onramp' in line:
                items = line.strip().split(':')
                workzone_topo['loc_onramp'] = sorted([float(i) for i in items[1].split(',')])
                workzone_topo['cell_onramp'] = [bisect.bisect(x_grid, i) - 1 for i in workzone_topo['loc_onramp']]

            elif 'loc_offramp' in line:
                items = line.strip().split(':')
                workzone_topo['loc_offramp'] = sorted([float(i) for i in items[1].split(',')])
                workzone_topo['cell_offramp'] = [bisect.bisect(x_grid, i) - 1 for i in workzone_topo['loc_offramp']]

            elif 'replications' in line:
                items = line.strip().split(':')
                workzone_topo['replications'] = [int(i) for i in items[1].split(',')]

        f_topo.close()

        # the following entries must be set in the topology file
        if t_grid is None or x_grid is None or workzone_topo['replications'] is None:
            raise Exception('Error: loc_start_end, simulation_duration, replications must be set in workzone topology.')

        # print('t_grid: {0}'.format(t_grid))
        # print('x_grid: {0}'.format(x_grid))

        print('workzone topo {0}'.format(workzone_topo))

        return t_grid, x_grid, workzone_topo, start_dur_step

    def load_config(self, config_file):
        """
        This function loads the configuration file including the sensor network specification and the
        algorithm specification.
        The configuration file can be generated from a separate script or manually configured in Excel.
        :param config_file: the file name
        :return:
        """

        # append configuration file to the all config log
        # config_log = self.__log_dir + '{0}_configurations_log.txt'.format(self.workzone)
        # file_config_log = open(config_log, 'a+')

        f_config = open(config_file, 'r')

        config_id = None
        for line in f_config:

            line = line.strip()

            # append line to all configuration file
            # file_config_log.write( line +'\n' )

            if len(line) == 0:
                # finished parsing this configuration
                config_id = None

            elif line[0] == '#':
                # skip all the comment
                continue

            else:

                items = line.split(';')
                line_att = items[0].split(':')[0].strip()

                if line_att == 'config_id':
                    # new configuration block
                    config_id = items[0].split(':')[1].strip()
                    self.all_config[config_id] = OrderedDict()

                elif line_att == 'sensor':
                    # directly save the attributes into a sensor key-value store
                    if items[1].split(':')[0] != 'id':
                        raise Exception('Error: The id of the sensor must come first in the definition of each sensor.')
                    else:
                        sensor_id = items[1].split(':')[1]
                        if 'sensors' not in self.all_config[config_id]:
                            self.all_config[config_id]['sensors'] = OrderedDict()

                        self.all_config[config_id]['sensors'][sensor_id] = OrderedDict()

                    for i in range(2, len(items)):
                        category = items[i].split(':')[0].strip()
                        # print('category {0}: {1}'.format(category, items[i]))
                        value = self.__assertType(items[i].split(':')[1].strip())
                        self.all_config[config_id]['sensors'][sensor_id][category] = value


                elif line_att == 'algorithm':

                    # directly save the attributes into a sensor key-value store
                    if items[1].split(':')[0] != 'id':
                        raise Exception('The id of the algorithm must come first in the definition of the algorithm.')
                    else:
                        alg_id = items[1].split(':')[1]
                        if 'algorithms' not in self.all_config[config_id]:
                            self.all_config[config_id]['algorithms'] = OrderedDict()
                        self.all_config[config_id]['algorithms'][alg_id] = OrderedDict()

                    for i in range(2, len(items)):
                        category = items[i].split(':')[0].strip()

                        value = self.__assertType(items[i].split(':')[1].strip())

                        self.all_config[config_id]['algorithms'][alg_id][category] = value

        f_config.close()
        # file_config_log.close()

        # post processing the loaded sensor:
        # - convert the relative location to absolute location
        # - find the cell of the sensors

        for config_id in self.all_config.keys():

            for s_id in self.all_config[config_id]['sensors'].keys():
                try:
                    self.all_config[config_id]['sensors'][s_id]['loc'] = \
                        self.__get_abs_loc(self.all_config[config_id]['sensors'][s_id]['section'],
                                           self.all_config[config_id]['sensors'][s_id]['distance'])

                    loc = self.all_config[config_id]['sensors'][s_id]['loc']
                    if type(loc) is float or type(loc) is int:
                        if loc in self.space_grid:
                            # on the freeway
                            self.all_config[config_id]['sensors'][s_id]['cell'] = self.space_grid.index(loc)
                        elif 'loc_onramp' in self.workzone_topo.keys() and self.workzone_topo['loc_onramp'] is not None \
                                and loc in self.workzone_topo['loc_onramp']:
                            # on the onramp
                            self.all_config[config_id]['sensors'][s_id]['cell'] = bisect.bisect(self.space_grid,
                                                                                                loc) - 1
                        elif 'loc_offramp' in self.workzone_topo.keys() and self.workzone_topo[
                            'loc_offramp'] is not None \
                                and loc in self.workzone_topo['loc_offramp']:
                            self.all_config[config_id]['sensors'][s_id]['cell'] = bisect.bisect(self.space_grid,
                                                                                                loc) - 1
                        else:
                            print('Location of sensor:{0}'.format(loc))
                            raise Exception('Error: flow sensors must locate on the grid or on/off ramps.')
                    elif type(loc) is tuple:
                        loc1, loc2 = loc
                        self.all_config[config_id]['sensors'][s_id]['cell'] = (self.space_grid.index(loc1),
                                                                               self.space_grid.index(loc2))

                except KeyError:
                    raise Exception(
                        'Error: section and distance must be specified for each sensor in the configuration file.')

        # all configurations loaded
        print('Status: All configurations loaded from file {0}'.format(config_file))

    def just_fetch(self):
        """
        This function directly calls Virtual_sensors class to generate all the virual sensor data.
        The sensor_to_be_generate file can be generated using the generate_sensors.py
        :return:
        """
        self.__data_generator.generate_virtual_sensors_data()

    def fetch_data(self):
        """
        This function fetches all the data needed to run all the configurations.
        This function merges the network specifications of all the configurations. Then it only request to virtual_sensor
        class to generate the data that has not been previously generated and saved.
        :return:
        """

        # check if data has been generated for each replication

        for rep in self.active_replications:
            # read generated sensor data record from self.sensor_data_record
            record_list = []
            if not exists(self.__sensor_data_log[rep]):
                print('Warning: No log of previously generated sensor data was found.')
            else:
                with open(self.__sensor_data_log[rep], 'r') as f_record:
                    for line in f_record:

                        if len(line) == 0 or line[0] == '#':
                            continue

                        # the first tuple is the sensor id
                        items = line.split(';')
                        line_att = items[0].split(':')[0].strip()
                        if line_att != 'id':
                            raise Exception('Error: The first entry of each sensor definition must be its id.')
                        else:
                            # sensor id can be a string
                            record_list.append(items[0].split(':')[1].strip())
                f_record.close()

                print(
                    'Status: The data for the following sensors has been previously generated: {0}'.format(record_list))

            # first remove the to_generate.txt in the previous round
            try:
                os.remove(self.__sensor_data_to_fetch[rep])
            except OSError:
                pass

            # check if the sensors in all the configurations have been generated before. If not, push to stack for generating data.
            f_to_fetch = None
            sensors_to_fetch = []
            for config in self.all_config.keys():
                for sensor_id in self.all_config[config]['sensors'].keys():
                    # find out new sensors with no previously generated data
                    if sensor_id not in record_list and sensor_id not in sensors_to_fetch:

                        sensors_to_fetch.append(sensor_id)
                        # build a string to write
                        if f_to_fetch is None:
                            f_to_fetch = open(self.__sensor_data_to_fetch[rep], 'w+')

                        # convert the sensor dictionary to a string and write to file
                        line = self.__sensorDictToString(sensor_id, self.all_config[config]['sensors'][sensor_id])
                        f_to_fetch.write('{0}\n'.format(line))

            # save to file
            if f_to_fetch is not None:
                f_to_fetch.close()

            if len(sensors_to_fetch) != 0:
                print('Status: To fetch data for replication {0}: {1}.'.format(rep, sensors_to_fetch))
            else:
                print('Status: All sensor data available for replication {0}.'.format(rep))

        # request virtual_sensor class to generate those data.
        self.__data_generator.generate_virtual_sensors_data()

    def generate_true_state(self, cus_grid, queue_threshold):
        """
        This function generates true states for a specific time-space resolution.
        The grid is used for:
            - Run the estimators. The resolution of the estimates may be different for the same sensor network and
             algorithm depending on the grid.
            - Compute the performance metric. The true states is extracted based on this grid from the trajectory data
        It first checks if the true states have been previously generated. If not, then generate.
        :param cus_grid: time x space resolution (s, m), specify this to generate the true states on
            this customized grid.
        :param queue_threshold: the speed in m/s, under which will be considered as congested.
        :return: trues states in metric units, m/s, m, s, saved in corresponding folder.
        """
        if cus_grid is None:
            cus_grid = self.grid_res

        # generate the true state for each replication
        for rep in self.active_replications:

            # check if the true states for a particular grid have been generated previously
            if not exists(self.__truestate_result_dir[rep] + self.__true_state_name(cus_grid)):
                print('Status: Generating true states for Replication {0} at grid {1}s x {2}m...'.format(rep,
                                                                                                         cus_grid[0],
                                                                                                         cus_grid[1]))
                self.__data_generator.generate_true_states_data(cus_grid, rep, queue_threshold)
            else:
                print(
                    'Status: True states for Replication {0} at grid {1}s x {2}m were generated.'.format(rep,
                                                                                                         cus_grid[0],
                                                                                                         cus_grid[1]))

    def run_estimators(self):
        """
        This function runs the algorithms specified in each configuration.
        :return: The estimation results, _speed.txt, _queue.txt, _traveltime.txt are saved in corresponding folders.
        """

        # run the estimators on each replication dataset
        for rep in self.active_replications:

            # read generated result record from self.result_record
            record_list = []
            if not exists(self.__result_log[rep]):
                print('Warning: No record of previously generated result was found.')
            else:
                with open(self.__result_log[rep], 'r') as f_log:
                    # find out the id of each configuration
                    for line in f_log:
                        if len(line) == 0 or line[0] == '#':
                            continue
                        else:
                            items = line.split('&')
                            record_list.append((items[0].strip(), items[1].strip()))
                f_log.close()

                print('Status: Previously generated result for replication {0}: {1}'.format(rep, record_list))

            for config_id in self.all_config.keys():

                # a configuration may potentially have multiple algorithms
                for alg_id in self.all_config[config_id]['algorithms'].keys():

                    if (config_id, alg_id) not in record_list:
                        # run estimators for this configuration and algorithm

                        # get the data for this configuration
                        flow_data, speed_data, traveltime_data, sensors = self.__prepare_data_for_estimation(config_id,
                                                                                                             rep)

                        if self.all_config[config_id]['algorithms'][alg_id]['type'] == 'interpolation':

                            print('Status: running estimator {0} for {1}...'.format(alg_id, config_id))
                            time0 = time.time()
                            # run estimator
                            estimator = self.__run_interpolation_estimator(flow_data, speed_data, traveltime_data,
                                                                           sensors,
                                                                           self.all_config[config_id]['algorithms'][
                                                                               alg_id]['interpolation_option'],
                                                                           self.all_config[config_id]['algorithms'][
                                                                               alg_id]['queue_threshold'],
                                                                           self.all_config[config_id]['algorithms'][
                                                                               alg_id]['missing_data'])
                            time1 = time.time()
                            print('      --- Took {0:.2f} s.'.format(time1 - time0))

                        elif self.all_config[config_id]['algorithms'][alg_id]['type'] == 'enkf_TFD' or \
                                        self.all_config[config_id]['algorithms'][alg_id]['type'] == 'enkf_QLFD' or \
                                        self.all_config[config_id]['algorithms'][alg_id]['type'] == 'enkf_AN':

                            print('Status: running estimator {0} for {1}'.format(alg_id, config_id))

                            time0 = time.time()

                            estimator = self.__run_EnKF_estimator(flow_data, speed_data, traveltime_data, sensors,
                                                                  self.all_config[config_id]['algorithms'][alg_id])
                            time1 = time.time()
                            print('      --- Took {0:.2f} s.'.format(time1 - time0))

                        elif self.all_config[config_id]['algorithms'][alg_id]['type'] == 'fisb':

                            print('Status: running estimator {0} for {1}'.format(alg_id, config_id))

                            time0 = time.time()
                            estimator = self.__run_FISB_estimator(flow_data, speed_data, traveltime_data, sensors,
                                                                  self.all_config[config_id]['algorithms'][alg_id])
                            time1 = time.time()
                            print('      --- Took {0:.2f} s.'.format(time1 - time0))

                        else:
                            raise Exception('Error: other algorithms not supported yet.')

                        # save the generated result into corresponding result_dir, and update result_record_file
                        result_name_prefix = self.__est_result_dir[rep] + '{0}_{1}'.format(config_id, alg_id)
                        # save the transposed density and speed: num_steps x num_cells on the finest grid
                        np.savetxt(result_name_prefix + "_density.txt", estimator.est_density.T, delimiter=",")
                        np.savetxt(result_name_prefix + "_speed.txt", estimator.est_speed.T, delimiter=",")

                        # Use a standard method for computing the queueu and travel time.
                        self.compute_queue_for_scenario(rep, config_id, alg_id,
                                                        self.all_config[config_id]['algorithms'][alg_id][
                                                            'queue_threshold'])
                        self.compute_traveltime_for_scenario(rep, config_id, alg_id)

                        # update the configuration_log
                        with open(self.__result_log[rep], 'a+') as f_log:
                            f_log.write('{0}&{1}\n'.format(config_id, alg_id))


    # ================================================================================ #
    # private functions used for running estimators
    # ================================================================================ #

    def __prepare_data_for_estimation(self, config_id, rep):
        """
        This function prepares the data set for the estimators. It reads the virtual sensor data from files and save in
        a dictionary for flow data, speed data, and travel time data.
            - Note, the timestamp in the virtual sensor data should be shifted by +30s, and the last entry
            should be removed. These are cleaned in this function.
            - The virtual sensor data files are saved in corresponding folders with columns:
                time (s), speed (kph), count
              This function converts the speed to m/s, and count to flow (veh/s). All algorithms internally use m, s,
              and m/s to do the computation.
        :param config_id: the configuration id
        :param rep: the replication id
        :return: flow_data, speed_data, traveltime_data,
                each is a dict with keys: 'time':[timstamps (s)], 'sensors':[sensor_ids (str)], 'data':2d array,
                flow_data['data'][index(sensor_id)][index(timestamp)] = data collected at timestamp by sensor_id
        """
        flow_data = {'time': [], 'sensors': [], 'data': []}
        speed_data = {'time': [], 'sensors': [], 'data': []}
        traveltime_data = {'time': [], 'sensors': [], 'data': []}
        sensors = {}

        for s_id in self.all_config[config_id]['sensors'].keys():

            data_file_name = self.__sensor_data_dir[rep] + '{0}.txt'.format(s_id)

            if not exists(data_file_name):
                raise Exception('Error: The data file for {0} does not exist.'.format(s_id))

            if self.all_config[config_id]['sensors'][s_id]['type'] == 'radar' or \
                            self.all_config[config_id]['sensors'][s_id]['type'] == 'rtms' or \
                            self.all_config[config_id]['sensors'][s_id]['type'] == 'icone' or \
                            self.all_config[config_id]['sensors'][s_id]['type'] == 'vol_gt':
                # if the file is RTMS, radar, or icone sensor, format:
                # timestamp, speed, flow
                flow_data['sensors'].append(s_id)
                speed_data['sensors'].append(s_id)

                # read data for the sensor from file
                sensor_data = np.genfromtxt(data_file_name, delimiter=',')

                # replace inf, -1, with nan
                # nan_idx = (sensor_data == np.inf) | (sensor_data < 0)
                # sensor_data[nan_idx] = np.nan
                sensor_data = self.__clean_data(sensor_data)

                # The current data format is shifted by 30 second
                dt_shift = self.all_config[config_id]['sensors'][s_id]['aggregation_sec']
                sensor_data[:, 0] += dt_shift

                # only get the data within the investigated time domain.
                time_index_max = np.searchsorted(sensor_data[:, 0], self.time_grid[-1], side='right')

                flow_data['time'] = sensor_data[0:time_index_max, 0].tolist()
                speed_data['time'] = sensor_data[0:time_index_max, 0].tolist()

                # change units speed: kph->m/s; count -> flow veh/s
                speed_data['data'].append((sensor_data[0:time_index_max, 1] * 1000.0 / 3600).tolist())
                flow_data['data'].append((sensor_data[0:time_index_max, 2] * 1.0 /
                                          self.all_config[config_id]['sensors'][s_id]['aggregation_sec']).tolist())

                # add the information to the sensor network
                sensors[s_id] = {}
                sensors[s_id]['loc'] = self.__get_abs_loc(self.all_config[config_id]['sensors'][s_id]['section'],
                                                          self.all_config[config_id]['sensors'][s_id]['distance'])
                sensors[s_id]['aggregation_sec'] = self.all_config[config_id]['sensors'][s_id]['aggregation_sec']
                sensors[s_id]['flow_std'] = self.all_config[config_id]['sensors'][s_id]['flow_std_specs']
                sensors[s_id]['speed_std'] = self.all_config[config_id]['sensors'][s_id]['speed_std_specs']

            elif self.all_config[config_id]['sensors'][s_id]['type'] == 'bluetooth' or \
                            self.all_config[config_id]['sensors'][s_id]['type'] == 'tt_gt':
                # if bluetooth sensor, then format:
                # time stamp, travel time
                traveltime_data['sensors'].append(s_id)

                # read data for the sensor from file
                sensor_data = np.genfromtxt(data_file_name, delimiter=',')

                # replace inf, -1, with nan
                # nan_idx = (sensor_data == np.inf) | (sensor_data < 0)
                # sensor_data[nan_idx] = np.nan
                sensor_data = self.__clean_data(sensor_data)

                # The current data format is shifted by 30 second
                dt_shift = self.all_config[config_id]['sensors'][s_id]['aggregation_sec']
                sensor_data[:, 0] += dt_shift

                # only get the data for the time under investigated
                time_index_max = np.searchsorted(sensor_data[:, 0], self.time_grid[-1], side='right')

                traveltime_data['time'] = sensor_data[0:time_index_max, 0].tolist()

                print('Status: prepared Sensor {0} data time from {1} to {2}'.format(s_id, sensor_data[0, 0],
                                                                                     sensor_data[-1, 0]))

                traveltime_data['data'].append(sensor_data[0:time_index_max, 1].tolist())

                # register information to the sensors
                sensors[s_id] = {}
                sensors[s_id]['loc'] = (self.__get_abs_loc(self.all_config[config_id]['sensors'][s_id]['section_up'],
                                                           self.all_config[config_id]['sensors'][s_id]['distance_up']),
                                        self.__get_abs_loc(self.all_config[config_id]['sensors'][s_id]['section_down'],
                                                           self.all_config[config_id]['sensors'][s_id][
                                                               'distance_down']))
                sensors[s_id]['traveltime_std'] = self.all_config[config_id]['sensors'][s_id]['traveltime_std_specs']

        if len(flow_data['sensors']) == 0:
            flow_data = None
        else:
            flow_data['data'] = np.matrix(flow_data['data'])

        if len(speed_data['sensors']) == 0:
            speed_data = None
        else:
            speed_data['data'] = np.matrix(speed_data['data'])

        if len(traveltime_data['sensors']) == 0:
            traveltime_data = None
        else:
            traveltime_data['data'] = np.matrix(traveltime_data['data'])

        return flow_data, speed_data, traveltime_data, sensors

    def __clean_data(self, data):
        """
        This function cleans the data.
            - It replaces the inf value by the last measurement. NOTE inf value means no measurment of a vehicle
            - It replaces -1 value by nan
        :param data: a 2d matrix, with each row as a new measurment: e.g. time, speed, flow
        :return: cleaned data, process data in the same matrix
        """
        # if no row or only one row, then no need to clean
        if data is None or data.shape[0] <= 1:
            return data

        # clean inf
        last_row = data[1, :]
        for row in range(1, data.shape[0]):
            # replace inf values by last value
            is_inf = np.isinf(data[row, :])
            data[row, is_inf] = last_row[is_inf]
            last_row = data[row, :]

        # clean -1
        data[data < 0] = np.nan

        return data

    def __run_interpolation_estimator(self, flow_data, speed_data, traveltime_data, sensors, option,
                                      queue_threshold=17.88,
                                      missing_data='fill'):
        """
        This function creates an instance of the interpolation estimator, and runs the estimator
        :param flow_data: the flow data prepared
        :param speed_data: the speed data prepared
        :param traveltime_data: the travel time data prepared
        :param sensors: list of sensor ids used for this estimator
        :param option: 'constant' or 'linear'
        :param queue_threshold: speed (m/s), the threshold for determining the queue
        :missing_data: 'fill', 'nan', or 'no_update', different ways for dealing with the missing data
        :return: estimation results are saved in csv in corresponding folders
        """

        estimator = Interpolation(self.time_grid, self.space_grid, self.workzone_topo['loc_onramp'],
                                  self.workzone_topo['loc_offramp'],
                                  sensors, queue_threshold, missing_data)

        estimator.set_meas_data(flow_data=flow_data, speed_data=speed_data, BT_data=traveltime_data)

        if option == 'constant':
            estimator.constant_interpolation()
        elif option == 'linear':
            estimator.linear_interpolation()
        else:
            raise Exception('Error: unrecognized interpolation method.')

        return estimator

    def __run_FISB_estimator(self, flow_data, speed_data, traveltime_data, sensors, alg_para):
        """
        This function creates an instance of the FISB estimator and runs the estimator
        :param flow_data: the flow data prepared by function __prepare_data_for_estimation
        :param speed_data: the speed data prepared
        :param traveltime_data: the travel time data prepared
        :param sensors: the list of sensors ids
        :param alg_para: the parameters for the fisb algorithm
        :return: estimation results are saved in csv in corresponding folders
        """
        omg_c = alg_para['omg_c']
        omg_f = alg_para['omg_f']
        sigma_fac = alg_para['sigma_fac']
        tau_fac = alg_para['tau_fac']
        v_crit = alg_para['v_crit']
        delta_v = alg_para['delta_v']
        lag = alg_para['lag']
        time_diff_cutoff = alg_para['time_diff_cutoff']
        queue_threshold = alg_para['queue_threshold']

        fisb = FISB_workzone(self.time_grid, self.space_grid, sensors,
                             self.workzone_topo['loc_onramp'], self.workzone_topo['loc_offramp'],
                             queue_threshold,
                             omg_c, omg_f, sigma_fac, tau_fac,
                             v_crit, delta_v, lag, time_diff_cutoff)
        fisb.set_meas_data(flow_data=flow_data, speed_data=speed_data, BT_data=traveltime_data)
        fisb.generate_estimates()

        return fisb

    def __run_EnKF_estimator(self, flow_data, speed_data, traveltime_data, sensors, alg_para):
        """
        :param flow_data: the flow data prepared by function __prepare_data_for_estimation
        :param speed_data: the speed data prepared
        :param traveltime_data: the travel time data prepared
        :param sensors: the list of sensors ids
        :param alg_para: the parameters for the fisb algorithm
        :return: estimation results are saved in csv in corresponding folders
        """

        if alg_para['type'] == 'enkf_AN':

            vm = alg_para['vm']
            beta = alg_para['beta']
            rhoc = alg_para['rhoc']
            wc = alg_para['wc']

            # the standard deviation of each parameter
            # Change the std to convariance assuming they are independent
            std_model_noise = OrderedDict()
            std_model_noise['cell'] = alg_para['std_cell']
            std_model_noise['oncell'] = alg_para['std_oncell']
            std_model_noise['offcell'] = alg_para['std_offcell']
            std_model_noise['qin'] = alg_para['std_qin']
            std_model_noise['qout'] = alg_para['std_qout']

            # print('std: {0}'.format([Q['vm'], Q['beta'], Q['rhoc'], Q['wc']]))

            num_ensembles = alg_para['num_ensembles']
            queue_threshold = alg_para['queue_threshold']

            init_rho = alg_para['init_rho']
            init_qin = alg_para['init_qin']
            init_qout = alg_para['init_qout']

            enkf = EnKF_AN(self.time_grid, self.space_grid, sensors,
                           self.workzone_topo['loc_onramp'], self.workzone_topo['loc_offramp'],
                           vm, beta, rhoc, wc,
                           num_ensembles,
                           std_model_noise, queue_threshold,
                           init_rho, init_qin, init_qout)

            enkf.set_meas_data(flow_data=flow_data, speed_data=speed_data, BT_data=traveltime_data)
            enkf.run_batch_filter()

            return enkf

        else:
            raise Exception('Deprecated or invalid enkf algorithm type: {0}'.format(alg_para['type']))

    # deprecated
    def convert_speed_kph2mps(self, rep, config_id, alg_id):
        """
        The FISB estimation result was originally saved in kph, convert to m/s
        :param rep: the repoiation id, int
        :param config_id: the configuraiotn id, str
        :param alg_id: the algorithm id, str
        :return:
        """
        # get the estimated result from this scenario
        est_speed_file = self.__est_result_dir[rep] + '{0}_{1}_speed.txt'.format(config_id, alg_id)
        est_speed = np.genfromtxt(est_speed_file, delimiter=',')
        speed_data = np.matrix(est_speed) / 3.6

        # save the updated speed file
        np.savetxt(est_speed_file, speed_data, delimiter=",")

    # ================================================================================ #
    # functions used for analyzing and visualizing the estimation results
    # ================================================================================ #
    # deprecated
    def final_compare_speed_error(self, replications, config_ids, alg_ids, grid, norm='L1',
                                  areas=('all', 'aroundqueue'), queue_buffer=800,
                                  x_axis='config', xlabel='Number of sensors', xticklabels=[],
                                  unit='metric', plot_style='line', title=None,
                                  save_fig=False, save_fig_name='speed_num_sensors.pdf',
                                  plot_cost_eff=False, cost_eff_title=None, cost_eff_num_sensors=None):
        """
        This function is the final function which compares the speed estimation errod and visualize them
        It will generate three plots. Each plots consists of several lines. Each line represent one algorithm. The x-axis
        is the configurations, while the y-axis is the error for the speed, queue, and travel time.
        :return:
        """

        # ==================================
        # get the I57 result
        # I57_speed_const, I57_queue_const, I57_traveltime_const = \
        #     self.__compute_L1_norm_for_scenario_against_grid(rep, 'configI57', 'constFILL', grid)
        # I57_speed_linear, I57_queue_linear, I57_traveltime_linear = self.__compute_L1_norm_for_scenario(rep, 'configI57', 'linearFILL')
        # I57_speed_enkf, I57_queue_enkf, I57_traveltime_enkf = self.__compute_L1_norm_for_scenario(rep, 'configI57', 'enkf50')

        if x_axis == 'alg' and plot_style == 'bar':
            raise Exception('bar plot across configurations is not supported.')

        # ==================================
        # compute the errors
        err_speed = OrderedDict()

        if x_axis == 'config':
            for alg_id in alg_ids:
                err_speed[alg_id] = OrderedDict()

                for area in areas:
                    err_speed[alg_id][area] = []

                    for config_id in config_ids:
                        # different types of error metric
                        L_speed_reps = []
                        for rep in replications:
                            L_speed = self.compute_speed_error_for_scenario_on_grid(rep, config_id, alg_id, grid,
                                                                                    norm=norm, area=area,
                                                                                    queue_buffer=queue_buffer)
                            L_speed_reps.append(L_speed)
                        # compute the mean error across replications
                        L_speed_mean = np.sum(L_speed_reps) / len(replications)
                        # print('appended {0}: {1}'.format( config_id ,L_speed_mean))
                        err_speed[alg_id][area].append(L_speed_mean)

                        # print('{0}+{1}: {2}'.format(alg_id, area, err_speed[alg_id][area]))

                print('Status: finished computing errors for {0}'.format(alg_id))

        elif x_axis == 'alg':
            for config_id in config_ids:
                err_speed[config_id] = OrderedDict()

                for area in areas:
                    err_speed[config_id][area] = []

                    for alg_id in alg_ids:
                        # different types of error metric
                        L_speed_reps = []
                        for rep in replications:
                            L_speed = self.compute_speed_error_for_scenario_on_grid(rep, config_id, alg_id, grid,
                                                                                    norm=norm, area=area,
                                                                                    queue_buffer=queue_buffer)
                            L_speed_reps.append(L_speed)
                        # compute the mean error across replications
                        L_speed_mean = np.sum(L_speed_reps) / len(replications)
                        err_speed[config_id][area].append(L_speed_mean)

        # ==================================
        # visualize the speed
        fig = plt.figure(figsize=(15, 10), dpi=100)
        ax = fig.add_axes([0.1, 0.15, 0.82, 0.75])
        ax2 = ax.twinx()

        if unit == 'imperial':
            # I57_speed_const = self.__metric2imperial(I57_speed_const, 'speed')
            # I57_speed_linear = self.__metric2imperial(I57_speed_linear, 'speed')
            # I57_speed_enkf = self.__metric2imperial(I57_speed_enkf, 'speed')
            unit_str = 'mph'
            for key in err_speed.keys():
                for a in err_speed[key].keys():
                    err_speed[key][a] = self.__metric2imperial(np.array(err_speed[key][a]), 'speed')
        elif unit == 'metric':
            unit_str = 'm/s'

        # print('VerMac speed: const:{0}; linear: {1}; enkf: {2}'.format(I57_speed_const, I57_speed_linear, I57_speed_enkf))
        # ax.scatter([800.0, 800.0], [I57_speed_linear, I57_speed_enkf], marker='o' )

        # setup bar plot
        if plot_style == 'bar':
            loc = np.arange(len(xticklabels)).astype(float)
            width = 0.2
            loc_interpolation = loc - width * 1.5
            loc_filtering = loc - width * 0.5
            loc_enkf = loc + width * 0.5

        for key in err_speed.keys():

            for mae_area in err_speed[key].keys():

                # set labels, linestyle, marker, color
                if x_axis == 'config':
                    if 'linear' in key:
                        col = 'r'
                        linestyle = ':'
                        if mae_area == 'aroundqueue':
                            label = 'Interpolation - around queue'
                            marker = '^'
                        else:
                            label = 'Interpolation'
                            marker = 'o'

                        if plot_style == 'bar':
                            loc = loc_interpolation
                            rect_interpolation = ax.bar(loc, err_speed[key][mae_area], width, color=col)

                    elif 'fisb' in key:
                        col = 'deepskyblue'
                        linestyle = '--'
                        if mae_area == 'aroundqueue':
                            label = 'Kernel filtering - around queue'
                            marker = '^'
                        else:
                            label = 'Kernel filtering'
                            marker = 'o'
                        if plot_style == 'bar':
                            loc = loc_filtering
                            rect_filtering = ax.bar(loc, err_speed[key][mae_area], width, color=col)

                    elif 'enkf' in key:
                        col = 'g'
                        linestyle = '-'
                        if mae_area == 'aroundqueue':
                            label = 'Kalman filter - around queue'
                            marker = '^'
                        else:
                            label = 'Kalman filter'
                            marker = 'o'

                        if plot_style == 'bar':
                            loc = loc_enkf
                            rect_enkf = ax.bar(loc, err_speed[key][mae_area], width, color=col)

                elif x_axis == 'alg':
                    label = key.split('_')[-1]
                else:
                    raise Exception('unrecognized x_axis')

                # ==============================================
                # plot the line
                if plot_style == 'line':

                    if mae_area == 'all':
                        ax.plot(err_speed[key][mae_area], label=label,
                                linewidth=2, marker=marker,
                                linestyle=linestyle, color=col,
                                markersize=10, fillstyle='full')
                    else:
                        ax2.plot(err_speed[key][mae_area], label=label,
                                 linewidth=2, marker=marker,
                                 linestyle=linestyle, color=col,
                                 markersize=10, fillstyle='full')
        if norm == 'L1':
            option_str = 'MAE'
        elif norm == 'L2':
            option_str = 'RMSE'
        else:
            raise Exception('Unrecognized norm')

        if title is None:
            ax.set_title('{0} norm of speed error', fontsize=32)
        else:
            ax.set_title(title, fontsize=32)

        ax.set_xlabel('{0}'.format(xlabel), fontsize=30)
        if plot_style == 'line':
            ax.set_xlim([0, len(xticklabels) - 1])
        elif plot_style == 'bar':
            ax.set_xlim([-0.5, len(xticklabels) - 0.5])
        x_ticks = np.arange(0, len(xticklabels))
        x_ticklabels = xticklabels

        ax.set_xticks(x_ticks)
        ax.set_xticklabels(x_ticklabels, fontsize=28)
        ax.tick_params(axis='both', which='major', labelsize=28)
        ax2.tick_params(axis='both', which='major', labelsize=28)

        ax.set_ylabel('Speed error ({0})'.format(unit_str), fontsize=28)
        ax.set_ylim([0, 25])
        ax2.set_ylabel('Error around queue ({0})'.format(unit_str), fontsize=28)
        ax2.set_ylim([5, 35])

        #        ax.grid(True)

        if plot_style == 'line':
            ax.legend(prop={'size': 20}, loc=2)
            ax2.legend(prop={'size': 20}, loc=6)
        elif plot_style == 'bar':
            ax.legend((rect_interpolation[0], rect_filtering[0], rect_enkf[0]),
                      ('Interpolaiton', 'Kernel filtering', 'Kalman filter'), prop={'size': 20}, loc='best')

        if save_fig is True:
            plt.savefig('{0}'.format(save_fig_name))
            plt.clf()
            plt.close()
        else:
            plt.draw()

        # ==================================
        # plot the cost-effectiveness: how much more improvement in average by adding each sensor
        if plot_cost_eff is True:
            fig = plt.figure(figsize=(15, 6), dpi=100)
            ax = fig.add_axes([0.1, 0.15, 0.82, 0.75])

            # print('VerMac speed: const:{0}; linear: {1}; enkf: {2}'.format(I57_speed_const, I57_speed_linear, I57_speed_enkf))
            # ax.scatter([800.0, 800.0], [I57_speed_linear, I57_speed_enkf], marker='o' )
            if cost_eff_num_sensors is None:
                print('Warning: Please specify cost_eff_num_sensors for cost-effectiveness plot')
                num_sensors = np.arange(0, len(xticklabels))
            else:
                num_sensors = np.array([float(i) for i in cost_eff_num_sensors])

            for key in err_speed.keys():

                if x_axis == 'config':
                    if 'linear' in key:
                        label = 'Interpolation'
                    elif 'fisb' in key:
                        label = 'Kernel filtering'
                    elif 'enkf' in key:
                        label = 'Kalman filter'
                elif x_axis == 'alg':
                    label = key.split('_')[-1]

                eff_per_sensor = []
                eff_per_sensor.append(100)
                for num in range(1, len(num_sensors)):
                    rel_improvment = -100.0 * (err_speed[key][num] - err_speed[key][num - 1]) / err_speed[key][num - 1]
                    eff_per_sensor.append(rel_improvment / (num_sensors[num] - num_sensors[num - 1]))

                ax.plot(eff_per_sensor, label=label, linewidth=2)

            if norm == 'L1':
                option_str = 'MAE'
            elif norm == 'L2':
                option_str = 'RMSE'

            if cost_eff_title is None:
                ax.set_title('Cost-effectiveness of speed {0}'.format(option_str), fontsize=32)
            else:
                ax.set_title(cost_eff_title, fontsize=32)

            ax.set_xlabel('{0}'.format(xlabel), fontsize=30)
            ax.set_xlim([0, len(xticklabels) - 1])

            x_ticks = np.arange(0, len(xticklabels))
            x_ticklabels = xticklabels

            ax.set_xticks(x_ticks)
            ax.set_xticklabels(x_ticklabels, fontsize=28)
            ax.tick_params(axis='both', which='major', labelsize=28)

            ax.set_ylabel('Improvement per sensor (%)', fontsize=30)
            ax.grid(True)

            plt.legend(prop={'size': 20})

            if save_fig is True:
                plt.savefig('cost_eff_{0}'.format(save_fig_name))
                plt.clf()
                plt.close()
            else:
                plt.draw()

    # deprecated
    def final_compare_queuelength_and_tt_error(self, replications, config_ids, alg_ids, grid,
                                               norm='L1', quantities=['queue', 'tt'],
                                               x_axis='config', xlabel='Number of sensors', xticklabels=[],
                                               unit='metric', plot_style='line', title=None,
                                               save_fig=False, save_fig_name='queue_num_sensors.pdf',
                                               plot_cost_eff=False, cost_eff_title=None,
                                               cost_eff_num_sensors=None):
        """
        This function computes the errors for the speed, queue, travel time.
        It will generate three plots. Each plots consists of several lines. Each line represent one algorithm. The x-axis
        is the configurations, while the y-axis is the error for the speed, queue, and travel time.
        :return:
        """

        # ==================================
        # get the I57 result
        # I57_speed_const, I57_queue_const, I57_traveltime_const = \
        #     self.__compute_L1_norm_for_scenario_against_grid(rep, 'configI57', 'constFILL', grid)
        # I57_speed_linear, I57_queue_linear, I57_traveltime_linear = self.__compute_L1_norm_for_scenario(rep, 'configI57', 'linearFILL')
        # I57_speed_enkf, I57_queue_enkf, I57_traveltime_enkf = self.__compute_L1_norm_for_scenario(rep, 'configI57', 'enkf50')

        if x_axis == 'alg' and plot_style == 'bar':
            raise Exception('bar plot across configurations is not supported.')

        # ==================================
        # compute the errors
        err_estimate = OrderedDict()

        if x_axis == 'config':
            for alg_id in alg_ids:
                err_estimate[alg_id] = OrderedDict()

                for quan in quantities:
                    err_estimate[alg_id][quan] = []

                    for config_id in config_ids:
                        L_est_reps = []
                        for rep in replications:
                            if quan == 'queue':
                                L_queue = \
                                    self.compute_queuelength_error_for_scenario_on_grid(rep,
                                                                                        config_id,
                                                                                        alg_id,
                                                                                        grid, norm)
                                L_est_reps.append(L_queue)
                            elif quan == 'tt':
                                L_tt = self.compute_traveltime_error_for_scenario_on_grid(rep, config_id, alg_id,
                                                                                          grid, norm)
                                L_est_reps.append(L_tt)
                            else:
                                raise Exception('unrecognized quantity')

                        # compute the mean queue error
                        L_est_mean = np.sum(L_est_reps) / len(replications)
                        err_estimate[alg_id][quan].append(L_est_mean)

                print('Status: finished computing errors for {0}'.format(alg_id))

        elif x_axis == 'alg':
            for config_id in config_ids:

                err_estimate[config_id] = OrderedDict()

                for quan in quantities:
                    err_estimate[config_id][quan] = []
                    for alg_id in alg_ids:
                        L_est_reps = []
                        for rep in replications:
                            if quan == 'queue':
                                L_queue = \
                                    self.compute_queuelength_error_for_scenario_on_grid(rep, config_id, alg_id, grid,
                                                                                        norm)
                                L_est_reps.append(L_queue)
                            elif quan == 'tt':
                                L_tt = self.compute_traveltime_error_for_scenario_on_grid(rep, config_id, alg_id,
                                                                                          grid, norm)
                                L_est_reps.append(L_tt)
                            else:
                                raise Exception('unrecognized quantity')

                        # compute the mean queue error
                        L_est_mean = np.sum(L_est_reps) / len(replications)
                        err_estimate[config_id][quan].append(L_est_mean)

        # ==================================
        # visualize the queue
        fig = plt.figure(figsize=(15, 10), dpi=100)
        ax = fig.add_axes([0.1, 0.15, 0.82, 0.75])
        ax2 = ax.twinx()

        if unit == 'imperial':
            # I57_queue_const = self.__metric2imperial(I57_queue_const, 'distance')
            # I57_queue_linear = self.__metric2imperial(I57_queue_linear, 'distance')
            # I57_queue_enkf = self.__metric2imperial(I57_queue_enkf, 'distance')
            unit_str = 'mile'
            for key in err_estimate.keys():
                for quan in err_estimate[key].keys():
                    if quan == 'queue':
                        err_estimate[key][quan] = self.__metric2imperial(np.array(
                            err_estimate[key][quan]), 'distance')
                    elif quan == 'tt':
                        # change time to mins
                        err_estimate[key][quan] = np.array(err_estimate[key][quan]) / 60.0

        elif unit == 'metric':
            unit_str = 'm'

        # print('VerMac queue: const:{0}; linear: {1}; enkf: {2}'.format(I57_queue_const, I57_queue_linear, I57_queue_enkf))

        # ax.scatter([800, 800], [I57_queue_linear, I57_queue_enkf], marker='o' )

        # setup bar plot
        if plot_style == 'bar':
            loc = np.arange(len(xticklabels)).astype(float)
            width = 0.2
            loc_interpolation = loc - width * 1.5
            loc_filtering = loc - width * 0.5
            loc_enkf = loc + width * 0.5

        # plot
        for key in err_estimate.keys():

            for quan in err_estimate[key].keys():

                # set labels, linestyle, marker, color
                if x_axis == 'config':
                    if 'linear' in key:
                        col = 'r'
                        linestyle = ':'
                        if quan == 'queue':
                            label = 'Interpolation - queue length'
                            marker = '^'
                        else:
                            label = 'Interpolation - travel time'
                            marker = 'o'

                        if plot_style == 'bar':
                            loc = loc_interpolation
                            rect_interpolation = ax.bar(loc, err_estimate[key][quan], width, color=col)

                    elif 'fisb' in key:
                        col = 'b'
                        linestyle = '--'
                        if quan == 'queue':
                            label = 'Kernel filtering - queue length'
                            marker = '^'
                        else:
                            label = 'Kernel filtering - travel time'
                            marker = 'o'
                        if plot_style == 'bar':
                            loc = loc_filtering
                            rect_filtering = ax.bar(loc, err_estimate[key][quan], width, color=col)

                    elif 'enkf' in key:
                        col = 'g'
                        linestyle = '-'
                        if quan == 'queue':
                            label = 'Kalman filter - queue length'
                            marker = '^'
                        else:
                            label = 'Kalman filter - travel time'
                            marker = 'o'

                        if plot_style == 'bar':
                            loc = loc_enkf
                            rect_enkf = ax.bar(loc, err_estimate[key][quan], width, color=col)

                elif x_axis == 'alg':
                    # label configurations, which mainly differ in sensors
                    label = key.split('_')[-1]
                else:
                    raise Exception('unrecognized x_axis')

                if plot_style == 'line':

                    if quan == 'queue':
                        # plot queue
                        ax.plot(err_estimate[key][quan], label=label,
                                linewidth=2, marker=marker,
                                linestyle=linestyle, color=col,
                                markersize=10, fillstyle='full')
                    else:
                        # print('err_estimate keys:{0}'.format(err_estimate.keys()))

                        # plot travel time
                        ax2.plot(err_estimate[key][quan], label=label,
                                 linewidth=2, marker=marker,
                                 linestyle=linestyle, color=col,
                                 markersize=10, fillstyle='full')

        if norm == 'L1':
            option_str = 'MAE'
        elif norm == 'L2':
            option_str = 'RMSE'

        if title is None:
            ax.set_title('{0} norm of estimates error'.format(option_str), fontsize=32)
        else:
            ax.set_title(title, fontsize=32)

        ax.set_xlabel('{0}'.format(xlabel), fontsize=30)
        if plot_style == 'line':
            ax.set_xlim([0, len(xticklabels) - 1])
        elif plot_style == 'bar':
            ax.set_xlim([-0.5, len(xticklabels) - 0.5])
        x_ticks = np.arange(0, len(xticklabels))
        x_ticklabels = xticklabels

        ax.set_xticks(x_ticks)
        ax.set_xticklabels(x_ticklabels, fontsize=28)
        ax.tick_params(axis='both', which='major', labelsize=28)
        ax2.tick_params(axis='both', which='major', labelsize=28)

        # set queue label
        ax.set_ylabel('Queue error ({0})'.format(unit_str), fontsize=28)
        ax.set_ylim([0, 1.5])
        # set travel time label
        ax2.set_ylabel('Travel time error (min)', fontsize=28)
        ax2.set_ylim([0, 30])

        # ax.grid(True)
        if plot_style == 'line':
            ax.legend(prop={'size': 20}, loc=2)
            ax2.legend(prop={'size': 20}, loc=1)
        elif plot_style == 'bar':
            ax.legend((rect_interpolation[0], rect_filtering[0], rect_enkf[0]),
                      ('Interpolaiton', 'Kernel filtering', 'Kalman filter'), prop={'size': 20})

        if save_fig is True:
            plt.savefig('{0}'.format(save_fig_name), bbox_inches='tight')
            plt.clf()
            plt.close()
        else:
            plt.draw()

        # ==================================
        # plot the cost-effectiveness: how much more improvement in average by adding each sensor
        if plot_cost_eff is True:
            fig = plt.figure(figsize=(15, 6), dpi=100)
            ax = fig.add_axes([0.1, 0.15, 0.82, 0.75])

            # print('VerMac speed: const:{0}; linear: {1}; enkf: {2}'.format(I57_speed_const, I57_speed_linear, I57_speed_enkf))
            # ax.scatter([800.0, 800.0], [I57_speed_linear, I57_speed_enkf], marker='o' )

            if cost_eff_num_sensors is None:
                print('Warning: Please specify cost_eff_num_sensors for cost-effectiveness plot')
                num_sensors = np.arange(0, len(xticklabels))
            else:
                num_sensors = np.array([float(i) for i in cost_eff_num_sensors])

            for key in err_estimate.keys():

                if x_axis == 'config':
                    if 'linear' in key:
                        label = 'Interpolation'
                    elif 'fisb' in key:
                        label = 'Kernel filtering'
                    elif 'enkf' in key:
                        label = 'Kalman filter'
                elif x_axis == 'alg':
                    label = key.split('_')[-1]

                eff_per_sensor = []
                eff_per_sensor.append(100)
                for num in range(1, len(num_sensors)):
                    rel_improvment = -100.0 * (err_estimate[key][num] - err_estimate[key][num - 1]) / err_estimate[key][
                        num - 1]
                    eff_per_sensor.append(rel_improvment / (num_sensors[num] - num_sensors[num - 1]))

                ax.plot(eff_per_sensor, label=label, linewidth=2)

            if norm == 'L1':
                option_str = 'MAE'
            elif norm == 'L2':
                option_str = 'RMSE'

            if cost_eff_title is None:
                ax.set_title('Cost-effectiveness of Queue length {0}'.format(option_str), fontsize=32)
            else:
                ax.set_title(title, fontsize=32)

            ax.set_xlabel('{0}'.format(xlabel), fontsize=30)
            ax.set_xlim([0, len(xticklabels) - 1])

            x_ticks = np.arange(0, len(xticklabels))
            x_ticklabels = xticklabels

            ax.set_xticks(x_ticks)
            ax.set_xticklabels(x_ticklabels, fontsize=28)
            ax.tick_params(axis='both', which='major', labelsize=28)

            ax.set_ylabel('Improvement per sensor (%)', fontsize=30)
            ax.grid(True)

            plt.legend(prop={'size': 20})

            if save_fig is True:
                plt.savefig('cost_eff_{0}'.format(save_fig_name), bbox_inches='tight')
                plt.clf()
                plt.close()
            else:
                plt.draw()

    def cost_compare_speed_error(self, replications, config_ids_prefix, alg_ids, grid, norm='L1',
                                 area='all', sensors=('RTMS', 'RADAR'), unit_costs=None, queue_buffer=800,
                                 x_axis='config', xlabel='Number of sensors', xticklabels=None,
                                 unit='metric', plot_style='line', title=None,
                                 save_fig=False, save_fig_name='speed_num_sensors.pdf',
                                 ylim=(0, 30), fontsize=(38, 36, 34), figsize=(10, 10)):
        """
        This function conducts the cost accuracy evaluation for the speed MAE. It computes and visualizes the cost
        accuracy for all configurations: config_ids_prefix + sensors

        :param replications: a list of replication ids, int
        :param config_ids_prefix: the configuration prefix, e.g. configHOMO_6,
            if added by _RTMS will be a configuration id.
        :param alg_ids: the algorithm ids.
        :param grid: the grid (s, m) of the true states used for computing the error.
        :param norm: 'L1', 'L2', the norm used.
        :param area: 'all', 'freeflow', 'congflow', 'aroundqueue'. Set this option to compute the speed error in
            corresponding areas.
        :param sensors: the list of sensors that will be compared in terms of cost effectiveness.
        :param unit_costs: a dictionary, unit_costs['sensor_type'] = $ (int)
        :param queue_buffer: meter +/- the position of the queue, used to determine the freeflow, congflow,
            and aroundqueue areas.
        :param x_axis: 'config' or 'alg', plot the results against x axis
        :param xlabel: str, the x label
        :param xticklabels: list of str, x tick labels
        :param unit: 'metric' (m/s), or 'imperial' (mph)
        :param plot_style: 'line' or 'bar'
        :param title: the tile of the figure
        :param save_fig: true or false, save the visualized figure
        :param save_fig_name: save the figure in this name
        :param ylim: the y limit in the unit
        :param fontsize: a tuple of three integers, corresponding to (title, label and legend, ticklabels)
        :param figsize: a tuple, the figure size
        :return: A plot of saved figured.
        """

        if x_axis == 'alg' and plot_style == 'bar':
            raise Exception('bar plot across configurations is not supported.')

        # ==================================
        # compute the errors
        err_speed = OrderedDict()
        cost_sys = OrderedDict()

        if x_axis == 'config':
            for alg_id in alg_ids:
                err_speed[alg_id] = OrderedDict()
                cost_sys[alg_id] = OrderedDict()

                for sensor in sensors:
                    err_speed[alg_id][sensor] = []
                    cost_sys[alg_id][sensor] = []

                    for config_id_prefix in config_ids_prefix:
                        # different types of error metric
                        config_id = config_id_prefix + '_' + sensor
                        _num_sensors = int(config_id_prefix.split('_')[1])

                        L_speed_reps = []
                        for rep in replications:
                            L_speed = self.compute_speed_error_for_scenario_on_grid(rep, config_id, alg_id, grid,
                                                                                    norm=norm, area=area,
                                                                                    queue_buffer=queue_buffer)
                            L_speed_reps.append(L_speed)
                        # compute the mean error across replications
                        L_speed_mean = np.sum(L_speed_reps) / len(replications)
                        # print('appended {0}: {1}'.format( config_id ,L_speed_mean))
                        err_speed[alg_id][sensor].append(L_speed_mean)

                        # compute the cost of this configuration
                        if unit_costs is not None:
                            cost_sys[alg_id][sensor].append(_num_sensors * unit_costs[sensor] / 1000.0)
                            # print('{0}+{1}: {2}'.format(alg_id, area, err_speed[alg_id][area]))
                            print('{0} + {1} {2}: cost, {3} k thousand; {4} speed error, {5} (mph)'.format(alg_id, _num_sensors, sensor,
                                                                                      cost_sys[alg_id][sensor][-1], area, 
                                                                                     self.__metric2imperial(err_speed[alg_id][sensor][-1], 'speed')))


                print('Status: finished computing errors for {0}'.format(alg_id))

        else:
            raise Exception('Cost analysis is only supported across configurations.')

        # ==================================
        # visualize the speed
        fig = plt.figure(figsize=figsize, dpi=100)
        ax = fig.add_axes([0.1, 0.15, 0.82, 0.75])

        if unit == 'imperial':
            unit_str = 'mph'
            for key in err_speed.keys():
                for a in err_speed[key].keys():
                    err_speed[key][a] = self.__metric2imperial(np.array(err_speed[key][a]), 'speed')
        elif unit == 'metric':
            unit_str = 'm/s'

        for key in err_speed.keys():
            # for each algorithm
            for sensor in err_speed[key].keys():
                # set labels, linestyle, marker, color
                if x_axis == 'config':
                    if 'linear' in key:
                        col = 'r'
                        linestyle = ':'
                        if sensor == 'RTMS':
                            label = 'Spatial + RTMS'
                            marker = '^'
                        else:
                            label = 'Spatial + RADAR'
                            marker = 'o'
                    elif 'fisb' in key:
                        col = 'deepskyblue'
                        linestyle = '--'
                        if sensor == 'RTMS':
                            label = 'Spatio-temporal + RTMS'
                            marker = '^'
                        else:
                            label = 'Spatio-temporal + RADAR'
                            marker = 'o'
                    elif 'enkf' in key:
                        col = 'g'
                        linestyle = '-'
                        if sensor == 'RTMS':
                            label = 'EnKF + RTMS'
                            marker = '^'
                        else:
                            label = 'EnKF + RADAR'
                            marker = 'o'
                else:
                    raise Exception('unrecognized x_axis')

                # ==============================================
                # plot the line
                if plot_style == 'line':
                    ax.plot(cost_sys[key][sensor], err_speed[key][sensor], label=label,
                            linewidth=2, marker=marker,
                            linestyle=linestyle, color=col,
                            markersize=10, fillstyle='full')
                    # print out the algorithm and sensor combinations, and the cost
        if norm == 'L1':
            option_str = 'MAE'
        elif norm == 'L2':
            option_str = 'RMSE'
        else:
            raise Exception('Unrecognized norm')

        if title is None:
            ax.set_title('{0} norm of speed error', fontsize=fontsize[0])
        else:
            ax.set_title(title, fontsize=fontsize[0])

        ax.set_xlabel('{0}'.format(xlabel), fontsize=fontsize[1])

        if xticklabels is not None:
            ax.set_xlim([0, len(xticklabels) - 1])
            x_ticks = np.arange(0, len(xticklabels))
            x_ticklabels = xticklabels
            ax.set_xticks(x_ticks)
            ax.set_xticklabels(x_ticklabels, fontsize=fontsize[2])
        ax.tick_params(axis='both', which='major', labelsize=fontsize[2])

        ax.set_ylabel('Velocity error ({0})'.format(unit_str), fontsize=fontsize[1])
        ax.set_ylim(ylim)

        if plot_style == 'line':
            ax.legend(prop={'size': fontsize[1]}, loc='best')
        else:
            raise Exception('Use line plot for cost analysis.')

        if save_fig is True:
            plt.savefig('{0}'.format(save_fig_name))
            plt.clf()
            plt.close()
        else:
            plt.draw()

    def cost_compare_queue_error(self, replications, config_ids_prefix, alg_ids, grid, norm='L1',
                                 sensors=('RTMS', 'RADAR'), unit_costs=None,
                                 x_axis='config', xlabel='Number of sensors', xticklabels=None,
                                 unit='metric', plot_style='line', title=None,
                                 save_fig=False, save_fig_name='speed_num_sensors.pdf',
                                 ylim=(0, 30), fontsize=(38, 36, 34), figsize=(10, 10)):
        """
        This function conducts the cost accuracy evaluation for the queue MAE. It computes and visualizes the cost
        accuracy for all configurations: config_ids_prefix + sensors

        :param replications: a list of replication ids, int
        :param config_ids_prefix: the configuration prefix, e.g. configHOMO_6,
            if added by _RTMS will be a configuration id.
        :param alg_ids: the algorithm ids.
        :param grid: the grid (s, m) of the true states used for computing the error.
        :param norm: 'L1', 'L2', the norm used.
        :param sensors: the list of sensors that will be compared in terms of cost effectiveness.
        :param unit_costs: a dictionary, unit_costs['sensor_type'] = $ (int)
        :param x_axis: 'config' or 'alg', plot the results against x axis
        :param xlabel: str, the x label
        :param xticklabels: list of str, x tick labels
        :param unit: 'metric' (m), or 'imperial' (miles)
        :param plot_style: 'line' or 'bar'
        :param title: the tile of the figure
        :param save_fig: true or false, save the visualized figure
        :param save_fig_name: save the figure in this name
        :param ylim: the y limit in the unit
        :param fontsize: a tuple of three integers, corresponding to (title, label and legend, ticklabels)
        :param figsize: a tuple, the figure size
        :return: A plot of saved figured.
        """

        if x_axis == 'alg' and plot_style == 'bar':
            raise Exception('bar plot across configurations is not supported.')

        # ==================================
        # compute the errors
        err_queue = OrderedDict()
        cost_sys = OrderedDict()

        if x_axis == 'config':
            for alg_id in alg_ids:
                err_queue[alg_id] = OrderedDict()
                cost_sys[alg_id] = OrderedDict()

                for sensor in sensors:
                    err_queue[alg_id][sensor] = []
                    cost_sys[alg_id][sensor] = []

                    for config_id_prefix in config_ids_prefix:
                        # different types of error metric
                        config_id = config_id_prefix + '_' + sensor
                        _num_sensors = int(config_id_prefix.split('_')[1])

                        L_queue_reps = []
                        for rep in replications:
                            L_queue = self.compute_queuelength_error_for_scenario_on_grid(rep, config_id, alg_id, grid,
                                                                                          norm=norm)
                            L_queue_reps.append(L_queue)
                        # compute the mean error across replications
                        L_mean = np.sum(L_queue_reps) / len(replications)

                        err_queue[alg_id][sensor].append(L_mean)

                        # compute the cost of this configuration
                        if unit_costs is not None:
                            cost_sys[alg_id][sensor].append(_num_sensors * unit_costs[sensor] / 1000.0)
                            # print('{0}+{1}: {2}'.format(alg_id, area, err_queue[alg_id][area]))
                            # print out the system error and system cost
                            print('{0} + {1} {2}: cost, {3} k dollars; queue error, {4} (mile)'.format(alg_id, _num_sensors, sensor,
                                                                                      cost_sys[alg_id][sensor][-1], 
                                                                        self.__metric2imperial(err_queue[alg_id][sensor][-1], 'distance')))

                print('Status: finished computing errors for {0}'.format(alg_id))

        else:
            raise Exception('Cost analysis is only supported across configurations.')

        # ==================================
        # visualize the queue
        fig = plt.figure(figsize=figsize, dpi=100)
        ax = fig.add_axes([0.1, 0.15, 0.82, 0.75])

        if unit == 'imperial':

            unit_str = 'mile'
            for key in err_queue.keys():
                for a in err_queue[key].keys():
                    err_queue[key][a] = self.__metric2imperial(np.array(err_queue[key][a]), 'distance')
        elif unit == 'metric':
            unit_str = 'm'

        for key in err_queue.keys():
            # for each algorithm
            for sensor in err_queue[key].keys():
                # set labels, linestyle, marker, color
                if x_axis == 'config':
                    if 'linear' in key:
                        col = 'r'
                        linestyle = ':'
                        if sensor == 'RTMS':
                            label = 'Spatial + RTMS'
                            marker = '^'
                        else:
                            label = 'Spatial + RADAR'
                            marker = 'o'
                    elif 'fisb' in key:
                        col = 'deepskyblue'
                        linestyle = '--'
                        if sensor == 'RTMS':
                            label = 'Spatio-temporal + RTMS'
                            marker = '^'
                        else:
                            label = 'Spatio-temporal + RADAR'
                            marker = 'o'
                    elif 'enkf' in key:
                        col = 'g'
                        linestyle = '-'
                        if sensor == 'RTMS':
                            label = 'EnKF + RTMS'
                            marker = '^'
                        else:
                            label = 'EnKF + RADAR'
                            marker = 'o'
                else:
                    raise Exception('unrecognized x_axis')

                # ==============================================
                # plot the line
                if plot_style == 'line':
                    print('ploting:\n{0}\n{1}'.format(cost_sys[key][sensor], err_queue[key][sensor]))
                    ax.plot(cost_sys[key][sensor], err_queue[key][sensor], label=label,
                            linewidth=2, marker=marker,
                            linestyle=linestyle, color=col,
                            markersize=10, fillstyle='full')

        if norm == 'L1':
            option_str = 'MAE'
        elif norm == 'L2':
            option_str = 'RMSE'
        else:
            raise Exception('Unrecognized norm')

        if title is None:
            ax.set_title('{0} norm of speed error', fontsize=fontsize[0])
        else:
            ax.set_title(title, fontsize=fontsize[0])

        ax.set_xlabel('{0}'.format(xlabel), fontsize=fontsize[1])

        if xticklabels is not None:
            ax.set_xlim([0, len(xticklabels) - 1])
            x_ticks = np.arange(0, len(xticklabels))
            x_ticklabels = xticklabels
            ax.set_xticks(x_ticks)
            ax.set_xticklabels(x_ticklabels, fontsize=fontsize[2])
        ax.tick_params(axis='both', which='major', labelsize=fontsize[2])

        ax.set_ylabel('Queue length error ({0})'.format(unit_str), fontsize=fontsize[1])
        ax.set_ylim(ylim)

        if plot_style == 'line':
            ax.legend(prop={'size': fontsize[1]}, loc='best')
        else:
            raise Exception('Use line plot for cost analysis.')

        if save_fig is True:
            plt.savefig('{0}'.format(save_fig_name))
            plt.clf()
            plt.close()
        else:
            plt.draw()

    def cost_compare_traveltime_error(self, replications, config_ids_prefix, alg_ids, grid, norm='L1',
                                      sensors=('RTMS', 'RADAR'), unit_costs=None,
                                      x_axis='config', xlabel='Number of sensors', xticklabels=None,
                                      unit='metric', plot_style='line', title=None,
                                      save_fig=False, save_fig_name='speed_num_sensors.pdf',
                                      ylim=(0, 30), fontsize=(38, 36, 34), figsize=(10, 10)):
        """
        This function conducts the cost accuracy evaluation for the travel time MAE. It computes and visualizes the cost
        accuracy for all configurations: config_ids_prefix + sensors

        :param replications: a list of replication ids, int
        :param config_ids_prefix: the configuration prefix, e.g. configHOMO_6,
            if added by _RTMS will be a configuration id.
        :param alg_ids: the algorithm ids.
        :param grid: the grid (s, m) of the true states used for computing the error.
        :param norm: 'L1', 'L2', the norm used.
        :param sensors: the list of sensors that will be compared in terms of cost effectiveness.
        :param unit_costs: a dictionary, unit_costs['sensor_type'] = $ (int)
        :param x_axis: 'config' or 'alg', plot the results against x axis
        :param xlabel: str, the x label
        :param xticklabels: list of str, x tick labels
        :param unit: 'metric' (s), or 'imperial' (s)
        :param plot_style: 'line' or 'bar'
        :param title: the tile of the figure
        :param save_fig: true or false, save the visualized figure
        :param save_fig_name: save the figure in this name
        :param ylim: the y limit in the unit
        :param fontsize: a tuple of three integers, corresponding to (title, label and legend, ticklabels)
        :param figsize: a tuple, the figure size
        :return: A plot of saved figured.
        """

        if x_axis == 'alg' and plot_style == 'bar':
            raise Exception('bar plot across configurations is not supported.')

        # ==================================
        # compute the errors
        err_tt = OrderedDict()
        cost_sys = OrderedDict()

        if x_axis == 'config':
            for alg_id in alg_ids:
                err_tt[alg_id] = OrderedDict()
                cost_sys[alg_id] = OrderedDict()

                for sensor in sensors:
                    err_tt[alg_id][sensor] = []
                    cost_sys[alg_id][sensor] = []

                    for config_id_prefix in config_ids_prefix:
                        # different types of error metric
                        config_id = config_id_prefix + '_' + sensor
                        _num_sensors = int(config_id_prefix.split('_')[1])

                        L_tt_reps = []
                        for rep in replications:
                            L_tt = self.compute_traveltime_error_for_scenario_on_grid(rep, config_id, alg_id, grid,
                                                                                         norm=norm)
                            L_tt_reps.append(L_tt)
                        # compute the mean error across replications
                        L_mean = np.sum(L_tt_reps) / len(replications)

                        err_tt[alg_id][sensor].append(L_mean)

                        # compute the cost of this configuration
                        if unit_costs is not None:
                            cost_sys[alg_id][sensor].append(_num_sensors * unit_costs[sensor] / 1000.0)
                            # print('{0}+{1}: {2}'.format(alg_id, area, err_queue[alg_id][area]))
                            print('{0} + {1} {2}: cost, {3} k dollars; travel time error, {4} (min)'.format(alg_id, _num_sensors, sensor,
                                                                                      cost_sys[alg_id][sensor][-1], err_tt[alg_id][sensor][-1]/60.0))


                print('Status: finished computing errors for {0}'.format(alg_id))

        else:
            raise Exception('Cost analysis is only supported across configurations.')

        # ==================================
        # convert travel time to mins
        for key in err_tt.keys():
            for sensor in err_tt[key].keys():
                err_tt[key][sensor] = np.array(err_tt[key][sensor]) / 60.0

        # ==================================
        # visualize the travel time
        fig = plt.figure(figsize=figsize, dpi=100)
        ax = fig.add_axes([0.1, 0.15, 0.82, 0.75])

        for key in err_tt.keys():
            # for each algorithm
            for sensor in err_tt[key].keys():
                # set labels, linestyle, marker, color
                if x_axis == 'config':
                    if 'linear' in key:
                        col = 'r'
                        linestyle = ':'
                        if sensor == 'RTMS':
                            label = 'Spatial + RTMS'
                            marker = '^'
                        else:
                            label = 'Spatial + RADAR'
                            marker = 'o'
                    elif 'fisb' in key:
                        col = 'deepskyblue'
                        linestyle = '--'
                        if sensor == 'RTMS':
                            label = 'Spatio-temporal + RTMS'
                            marker = '^'
                        else:
                            label = 'Spatio-temporal + RADAR'
                            marker = 'o'
                    elif 'enkf' in key:
                        col = 'g'
                        linestyle = '-'
                        if sensor == 'RTMS':
                            label = 'EnKF + RTMS'
                            marker = '^'
                        else:
                            label = 'EnKF + RADAR'
                            marker = 'o'
                else:
                    raise Exception('unrecognized x_axis')

                # ==============================================
                # plot the line
                if plot_style == 'line':
                    print('plotting travel time error:\n{0}\n{1}'.format(cost_sys[key][sensor], err_tt[key][sensor]))
                    ax.plot(cost_sys[key][sensor], err_tt[key][sensor], label=label,
                            linewidth=2, marker=marker,
                            linestyle=linestyle, color=col,
                            markersize=10, fillstyle='full')

        if norm == 'L1':
            option_str = 'MAE'
        elif norm == 'L2':
            option_str = 'RMSE'
        else:
            raise Exception('Unrecognized norm')

        if title is None:
            ax.set_title('{0} norm of speed error', fontsize=fontsize[0])
        else:
            ax.set_title(title, fontsize=fontsize[0])

        ax.set_xlabel('{0}'.format(xlabel), fontsize=fontsize[1])

        if xticklabels is not None:
            ax.set_xlim([0, len(xticklabels) - 1])
            x_ticks = np.arange(0, len(xticklabels))
            x_ticklabels = xticklabels
            ax.set_xticks(x_ticks)
            ax.set_xticklabels(x_ticklabels, fontsize=fontsize[2])
        ax.tick_params(axis='both', which='major', labelsize=fontsize[2])

        ax.set_ylabel('Travel time error (min)', fontsize=fontsize[1])
        ax.set_ylim(ylim)

        if plot_style == 'line':
            ax.legend(prop={'size': fontsize[1]}, loc='best')
        else:
            raise Exception('Use line plot for cost analysis.')

        if save_fig is True:
            plt.savefig('{0}'.format(save_fig_name))
            plt.clf()
            plt.close()
        else:
            plt.draw()

    def compare_speed_error(self, replications, config_ids, alg_ids, grid, norm='L1', area='all', queue_buffer=800,
                            x_axis='config', xlabel='Number of sensors', xticklabels=(), ff_speed=27.2,
                            unit='metric', plot_style='line', title=None, ylim=(0, 30), fontsize=(38, 36, 34),
                            save_fig=False, save_fig_name='speed_num_sensors.pdf', figsize=(10, 10),
                            plot_cost_eff=False, cost_eff_title=None, cost_eff_num_sensors=None,
                            cost_eff_fontsize=(38, 36, 34)):
        """
        This function compares and visualizes the speed MAE in corresponding area

        :param replications: a list of replication ids, int
        :param config_ids: the configuration ids.
        :param alg_ids: the algorithm ids.
        :param grid: the grid (s, m) of the true states used for computing the error.
        :param norm: 'L1', 'L2', the norm used.
        :param area: 'all', 'freeflow', 'congflow', 'aroundqueue'. Set this option to compute the speed error in
            corresponding areas.
        :param queue_buffer: meter +/- the position of the queue, used to determine the freeflow, congflow,
            and aroundqueue areas.
        :param x_axis: 'config' or 'alg', plot the results against x axis
        :param xlabel: str, the x label
        :param xticklabels: list of str, x tick labels
        :param ff_speed: m/s, the maximum free flow speed.
        :param unit: 'metric' (m/s), or 'imperial' (mph)
        :param plot_style: 'line' or 'bar'
        :param title: the tile of the figure
        :param save_fig: true or false, save the visualized figure
        :param save_fig_name: save the figure in this name
        :param ylim: the y limit in the unit
        :param fontsize: a tuple of three integers, corresponding to (title, label and legend, ticklabels)
        :param figsize: a tuple, the figure size
        :param plot_cost_eff: True or False, whether to plot the derivative of the MAE curves.
        :param cost_eff_title: str, the title for the cost_eff figure
        :param cost_eff_num_sensors: the numebr of sensors in each configuration for computing the derivative.
        :param cost_eff_fontsize: corresponds to (title, label and legend, ticklabels)
        :return: A plot of saved figured.
        """

        if x_axis == 'alg' and plot_style == 'bar':
            raise Exception('bar plot across configurations is not supported.')

        # ==================================
        # compute the errors
        err_speed = OrderedDict()

        if x_axis == 'config':
            for alg_id in alg_ids:
                err_speed[alg_id] = []

                for config_id in config_ids:
                    # different types of error metric
                    L_speed_reps = []
                    for rep in replications:
                        L_speed = self.compute_speed_error_for_scenario_on_grid(rep, config_id, alg_id, grid,
                                                                                norm=norm, area=area,
                                                                                queue_buffer=queue_buffer,
                                                                                max_speed=ff_speed)
                        L_speed_reps.append(L_speed)

                    # compute the mean error across replications
                    L_speed_mean = np.sum(L_speed_reps) / len(replications)
                    err_speed[alg_id].append(L_speed_mean)

        elif x_axis == 'alg':
            for config_id in config_ids:
                err_speed[config_id] = []

                for alg_id in alg_ids:
                    # different types of error metric
                    L_speed_reps = []
                    for rep in replications:
                        L_speed = self.compute_speed_error_for_scenario_on_grid(rep, config_id, alg_id, grid,
                                                                                norm=norm, area=area,
                                                                                queue_buffer=queue_buffer,
                                                                                max_speed=ff_speed)
                        L_speed_reps.append(L_speed)
                    # compute the mean error across replications
                    L_speed_mean = np.sum(L_speed_reps) / len(replications)
                    err_speed[config_id].append(L_speed_mean)

        # ==================================
        # visualize the speed
        fig = plt.figure(figsize=figsize, dpi=100)
        ax = fig.add_axes([0.1, 0.15, 0.82, 0.75])

        if unit == 'imperial':
            unit_str = 'mph'
            for key in err_speed.keys():
                err_speed[key] = self.__metric2imperial(np.array(err_speed[key]), 'speed')
        elif unit == 'metric':
            unit_str = 'm/s'

        # setup bar plot
        if plot_style == 'bar':
            loc = np.arange(len(xticklabels)).astype(float)
            width = 0.2
            loc_interpolation = loc - width * 1.5
            loc_filtering = loc - width * 0.5
            loc_enkf = loc + width * 0.5

        for key in err_speed.keys():

            if x_axis == 'config':
                if 'linear' in key:
                    label = 'Spatial'
                    marker = '^'
                    col = 'r'
                    if plot_style == 'bar':
                        loc = loc_interpolation
                        rect_interpolation = ax.bar(loc, err_speed[key], width, color=col, hatch='/')
                elif 'fisb' in key:
                    label = 'Spatio-temporal'
                    marker = 'o'
                    col = 'deepskyblue'
                    if plot_style == 'bar':
                        loc = loc_filtering
                        rect_filtering = ax.bar(loc, err_speed[key], width, color=col, hatch='*')
                elif 'enkf' in key:
                    label = 'EnKF'
                    col = 'g'
                    marker = '*'
                    if plot_style == 'bar':
                        loc = loc_enkf
                        rect_enkf = ax.bar(loc, err_speed[key], width, color=col, hatch='x')
            elif x_axis == 'alg':
                label = key.split('_')[-1]
            else:
                raise Exception('unrecognized x_axis')

            if plot_style == 'line':
                ax.plot(err_speed[key], label=label, linewidth=2, color=col,
                        marker=marker, markersize=10, fillstyle='full')

        if norm == 'L1':
            option_str = 'MAE'
        elif norm == 'L2':
            option_str = 'RMSE'
        else:
            raise Exception('Unrecognized norm')

        if title is None:
            ax.set_title('{0} norm of {1} speed error'.format(option_str, area), fontsize=fontsize[0])
        else:
            ax.set_title(title, fontsize=fontsize[0])

        ax.set_xlabel('{0}'.format(xlabel), fontsize=fontsize[1])
        if plot_style == 'line':
            ax.set_xlim([0, len(xticklabels) - 1])
        elif plot_style == 'bar':
            ax.set_xlim([-0.5, len(xticklabels) - 0.5])
        x_ticks = np.arange(0, len(xticklabels))
        x_ticklabels = xticklabels

        ax.set_xticks(x_ticks)
        ax.set_xticklabels(x_ticklabels, fontsize=fontsize[2])
        ax.tick_params(axis='both', which='major', labelsize=fontsize[2])

        ax.set_ylabel('Velocity error ({0})'.format(unit_str), fontsize=fontsize[1])
        ax.set_ylim(ylim)
        # ax.grid(True)

        if plot_style == 'line':
            plt.legend(prop={'size': fontsize[1]})
        elif plot_style == 'bar':
            ax.legend((rect_interpolation[0], rect_filtering[0], rect_enkf[0]),
                      ('Spatial', 'Spatio-temporal', 'EnKF'), prop={'size': fontsize[1]}, loc='best')

        if save_fig is True:
            plt.savefig('{0}'.format(save_fig_name), bbox_inches='tight')
            plt.clf()
            plt.close()
        else:
            plt.draw()

        # ================================================================================
        # plot the cost-effectiveness: how much more improvement in average by adding each sensor
        # ================================================================================
        if plot_cost_eff is True:
            fig = plt.figure(figsize=figsize, dpi=100)
            ax = fig.add_axes([0.1, 0.15, 0.82, 0.75])

            if cost_eff_num_sensors is None:
                print('Warning: Please specify cost_eff_num_sensors for cost-effectiveness plot')
                num_sensors = np.arange(0, len(xticklabels))
            else:
                num_sensors = np.array([float(i) for i in cost_eff_num_sensors])

            for key in err_speed.keys():

                if x_axis == 'config':
                    if 'linear' in key:
                        label = 'Spatial'
                        marker = '^'
                        col = 'r'
                    elif 'fisb' in key:
                        label = 'Spatio-temporal'
                        marker = 'o'
                        col = 'deepskyblue'
                    elif 'enkf' in key:
                        label = 'EnKF'
                        col = 'g'
                        marker = '*'
                elif x_axis == 'alg':
                    label = key.split('_')[-1]

                eff_per_sensor = []
                eff_per_sensor.append(100)
                for num in range(1, len(num_sensors)):
                    rel_improvment = -100.0 * (err_speed[key][num] -
                                               err_speed[key][num - 1]) / err_speed[key][num - 1]
                    eff_per_sensor.append(rel_improvment / (num_sensors[num] - num_sensors[num - 1]))

                ax.plot(eff_per_sensor, label=label, linewidth=2,
                        marker=marker, color=col, markersize=10)

            if norm == 'L1':
                option_str = 'MAE'
            elif norm == 'L2':
                option_str = 'RMSE'

            if cost_eff_title is None:
                ax.set_title('Improvement of velocity {0}'.format(option_str),
                             fontsize=cost_eff_fontsize[0])
            else:
                ax.set_title(cost_eff_title, fontsize=cost_eff_fontsize[0])

            ax.set_xlabel('{0}'.format(xlabel), fontsize=cost_eff_fontsize[1])
            ax.set_xlim([0, len(xticklabels) - 1])

            x_ticks = np.arange(0, len(xticklabels))
            x_ticklabels = xticklabels

            ax.set_xticks(x_ticks)
            ax.set_xticklabels(x_ticklabels, fontsize=cost_eff_fontsize[2])
            ax.tick_params(axis='both', which='major', labelsize=cost_eff_fontsize[2])

            ax.set_ylabel('Improvement per sensor (%)', fontsize=cost_eff_fontsize[1])
            # ax.grid(True)

            plt.legend(prop={'size': cost_eff_fontsize[1]})

            if save_fig is True:
                plt.savefig('cost_eff_{0}'.format(save_fig_name), bbox_inches='tight')
                plt.clf()
                plt.close()
            else:
                plt.draw()

    def compare_queuelength_error(self, replications, config_ids, alg_ids, grid, norm='L1',
                                  x_axis='config', xlabel='Number of sensors', xticklabels=(),
                                  unit='metric', plot_style='line', title=None, figsize=(10, 10),
                                  save_fig=False, save_fig_name='queue_num_sensors.pdf', ylim=(0, 2),
                                  fontsize=(38, 36, 34),
                                  plot_cost_eff=False, cost_eff_title=None, cost_eff_num_sensors=None,
                                  cost_eff_fontsize=(38, 36, 34)):
        """
        This function compares and visualizes the queue length MAE

        :param replications: a list of replication ids, int
        :param config_ids: the configuration ids.
        :param alg_ids: the algorithm ids.
        :param grid: the grid (s, m) of the true states used for computing the error.
        :param norm: 'L1', 'L2', the norm used.
        :param x_axis: 'config' or 'alg', plot the results against x axis
        :param xlabel: str, the x label
        :param xticklabels: list of str, x tick labels
        :param unit: 'metric' (m), or 'imperial' (miles)
        :param plot_style: 'line' or 'bar'
        :param title: the tile of the figure
        :param save_fig: true or false, save the visualized figure
        :param save_fig_name: save the figure in this name
        :param ylim: the y limit in the unit
        :param fontsize: a tuple of three integers, corresponding to (title, label and legend, ticklabels)
        :param figsize: a tuple, the figure size
        :param plot_cost_eff: True or False, whether to plot the derivative of the MAE curves.
        :param cost_eff_title: str, the title for the cost_eff figure
        :param cost_eff_num_sensors: the numebr of sensors in each configuration for computing the derivative.
        :param cost_eff_fontsize: corresponds to (title, label and legend, ticklabels)
        :return: A plot of saved figured.
        """

        if x_axis == 'alg' and plot_style == 'bar':
            raise Exception('bar plot across configurations is not supported.')

        # ==================================
        # compute the errors
        err_queue = OrderedDict()

        if x_axis == 'config':
            for alg_id in alg_ids:
                err_queue[alg_id] = []
                for config_id in config_ids:

                    L_queue_reps = []
                    for rep in replications:
                        L_queue = \
                            self.compute_queuelength_error_for_scenario_on_grid(rep, config_id, alg_id, grid, norm)
                        L_queue_reps.append(L_queue)

                    # compute the mean queue error
                    L_queue_mean = np.sum(L_queue_reps) / len(replications)
                    err_queue[alg_id].append(L_queue_mean)

        elif x_axis == 'alg':
            for config_id in config_ids:
                err_queue[config_id] = []
                for alg_id in alg_ids:
                    L_queue_reps = []
                    for rep in replications:
                        L_queue = \
                            self.compute_queuelength_error_for_scenario_on_grid(rep, config_id, alg_id, grid, norm)
                        L_queue_reps.append(L_queue)

                    # compute the mean queue error
                    L_queue_mean = np.sum(L_queue_reps) / len(replications)
                    err_queue[config_id].append(L_queue_mean)

        # ==================================
        # visualize the queue
        fig = plt.figure(figsize=figsize, dpi=100)
        ax = fig.add_axes([0.1, 0.15, 0.82, 0.75])

        if unit == 'imperial':

            unit_str = 'mile'
            for key in err_queue.keys():
                err_queue[key] = self.__metric2imperial(np.array(err_queue[key]), 'distance')
        elif unit == 'metric':
            unit_str = 'm'

        # setup bar plot
        if plot_style == 'bar':
            loc = np.arange(len(xticklabels)).astype(float)
            width = 0.2
            loc_interpolation = loc - width * 1.5
            loc_filtering = loc - width * 0.5
            loc_enkf = loc + width * 0.5

        # set line labels
        for key in err_queue.keys():
            if x_axis == 'config':
                if 'linear' in key:
                    label = 'Spatial'
                    marker = '^'
                    col = 'r'
                    if plot_style == 'bar':
                        loc = loc_interpolation
                        col = 'r'
                        print('current color interpolation: {0}'.format(col))
                        rect_interpolation = ax.bar(loc, err_queue[key], width, color=col, hatch='/')
                elif 'fisb' in key:
                    label = 'Spatio-temporal'
                    marker = 'o'
                    col = 'deepskyblue'
                    if plot_style == 'bar':
                        loc = loc_filtering
                        rect_filtering = ax.bar(loc, err_queue[key], width, color=col, hatch='*')
                elif 'enkf' in key:
                    label = 'EnKF'
                    col = 'g'
                    marker = '*'
                    if plot_style == 'bar':
                        loc = loc_enkf
                        rect_enkf = ax.bar(loc, err_queue[key], width, color=col, hatch='x')
            elif x_axis == 'alg':
                # label configurations, which mainly differ in sensors
                label = key.split('_')[-1]
            else:
                raise Exception('unrecognized x_axis')

            if plot_style == 'line':
                ax.plot(err_queue[key], label=label, linewidth=2,
                        marker=marker, markersize=10, fillstyle='full', color=col)

        if norm == 'L1':
            option_str = 'MAE'
        elif norm == 'L2':
            option_str = 'RMSE'

        if title is None:
            ax.set_title('{0} norm of queue error'.format(option_str), fontsize=fontsize[0])
        else:
            ax.set_title(title, fontsize=fontsize[0])

        ax.set_xlabel('{0}'.format(xlabel), fontsize=fontsize[1])
        if plot_style == 'line':
            ax.set_xlim([0, len(xticklabels) - 1])
        elif plot_style == 'bar':
            ax.set_xlim([-0.5, len(xticklabels) - 0.5])

        x_ticks = np.arange(0, len(xticklabels))
        x_ticklabels = xticklabels

        ax.set_xticks(x_ticks)
        ax.set_xticklabels(x_ticklabels, fontsize=fontsize[2])
        ax.tick_params(axis='both', which='major', labelsize=fontsize[2])

        ax.set_ylabel('queue error ({0})'.format(unit_str), fontsize=fontsize[1])
        ax.set_ylim(ylim)

        # ax.grid(True)
        if plot_style == 'line':
            plt.legend(prop={'size': fontsize[1]})
        elif plot_style == 'bar':
            ax.legend((rect_interpolation[0], rect_filtering[0], rect_enkf[0]),
                      ('Spatial', 'Spatio-temporal', 'EnKF'), prop={'size': fontsize[1]})

        if save_fig is True:
            plt.savefig('{0}'.format(save_fig_name), bbox_inches='tight')
            plt.clf()
            plt.close()
        else:
            plt.draw()

        # ================================================================================
        # plot the cost-effectiveness: how much more improvement in average by adding each sensor
        # ================================================================================
        if plot_cost_eff is True:
            fig = plt.figure(figsize=figsize, dpi=100)
            ax = fig.add_axes([0.1, 0.15, 0.82, 0.75])

            if cost_eff_num_sensors is None:
                print('Warning: Please specify cost_eff_num_sensors for cost-effectiveness plot')
                num_sensors = np.arange(0, len(xticklabels))
            else:
                num_sensors = np.array([float(i) for i in cost_eff_num_sensors])

            for key in err_queue.keys():

                if x_axis == 'config':
                    if 'linear' in key:
                        label = 'Spatial'
                        marker = '^'
                        col = 'r'
                    elif 'fisb' in key:
                        label = 'Spatio-temporal'
                        marker = 'o'
                        col = 'deepskyblue'
                    elif 'enkf' in key:
                        label = 'EnKF'
                        col = 'g'
                        marker = '*'
                elif x_axis == 'alg':
                    label = key.split('_')[-1]

                eff_per_sensor = []
                eff_per_sensor.append(100)
                for num in range(1, len(num_sensors)):
                    rel_improvment = -100.0 * (err_queue[key][num] -
                                               err_queue[key][num - 1]) / err_queue[key][num - 1]
                    eff_per_sensor.append(rel_improvment / (num_sensors[num] - num_sensors[num - 1]))

                ax.plot(eff_per_sensor, label=label, linewidth=2,
                        marker=marker, color=col, markersize=10)

            if norm == 'L1':
                option_str = 'MAE'
            elif norm == 'L2':
                option_str = 'RMSE'

            if cost_eff_title is None:
                ax.set_title('Cost-effectiveness of Queue length {0}'.format(option_str),
                             fontsize=cost_eff_fontsize[0])
            else:
                ax.set_title(cost_eff_title, fontsize=cost_eff_fontsize[0])

            ax.set_xlabel('{0}'.format(xlabel), fontsize=cost_eff_fontsize[1])
            ax.set_xlim([0, len(xticklabels) - 1])

            x_ticks = np.arange(0, len(xticklabels))
            x_ticklabels = xticklabels

            ax.set_xticks(x_ticks)
            ax.set_xticklabels(x_ticklabels, fontsize=cost_eff_fontsize[2])
            ax.tick_params(axis='both', which='major', labelsize=cost_eff_fontsize[2])

            ax.set_ylabel('Improvement per sensor (%)', fontsize=cost_eff_fontsize[1])
            # ax.grid(True)

            plt.legend(prop={'size': cost_eff_fontsize[1]})

            # plt.gca().tight_layout()
            if save_fig is True:
                plt.savefig('cost_eff_{0}'.format(save_fig_name), bbox_inches='tight')
                plt.clf()
                plt.close()
            else:
                plt.draw()

    def compare_traveltime_error(self, replications, config_ids, alg_ids, grid, norm='L2',
                                 x_axis='config', xlabel='Number of sensors', xticklabels=None,
                                 unit='metric', plot_style='line', figsize=(10, 10),
                                 save_fig=False, save_fig_name='traveltime_num_sensors.pdf',
                                 ylim=(0, 30), fontsize=(38, 36, 34), title=None,
                                 plot_cost_eff=False, cost_eff_title=None, cost_eff_fontsize=(38, 36, 34),
                                 cost_eff_num_sensors=None):
        """
        This function compares and visualizes the travel time MAE

        :param replications: a list of replication ids, int
        :param config_ids: the configuration ids.
        :param alg_ids: the algorithm ids.
        :param grid: the grid (s, m) of the true states used for computing the error.
        :param norm: 'L1', 'L2', the norm used.
        :param x_axis: 'config' or 'alg', plot the results against x axis
        :param xlabel: str, the x label
        :param xticklabels: list of str, x tick labels
        :param unit: 'metric' (s), or 'imperial' (s)
        :param plot_style: 'line' or 'bar'
        :param title: the tile of the figure
        :param save_fig: true or false, save the visualized figure
        :param save_fig_name: save the figure in this name
        :param ylim: the y limit in the unit
        :param fontsize: a tuple of three integers, corresponding to (title, label and legend, ticklabels)
        :param figsize: a tuple, the figure size
        :param plot_cost_eff: True or False, whether to plot the derivative of the MAE curves.
        :param cost_eff_title: str, the title for the cost_eff figure
        :param cost_eff_num_sensors: the numebr of sensors in each configuration for computing the derivative.
        :param cost_eff_fontsize: corresponds to (title, label and legend, ticklabels)
        :return: A plot of saved figured.
        """

        if x_axis == 'alg' and plot_style == 'bar':
            raise Exception('bar plot across configurations is not supported.')

        # compute the errors
        err_traveltime = OrderedDict()

        # the x-axis is in default config
        if xticklabels is None:
            xticklabels = config_ids

        if x_axis == 'config':
            for alg_id in alg_ids:
                err_traveltime[alg_id] = []
                for config_id in config_ids:

                    L_tt_reps = []
                    for rep in replications:
                        L_traveltime = \
                            self.compute_traveltime_error_for_scenario_on_grid(rep, config_id,
                                                                               alg_id, grid, norm=norm)
                        L_tt_reps.append(L_traveltime)

                    # compute the mean
                    L_tt_mean = np.sum(L_tt_reps) / len(replications)
                    err_traveltime[alg_id].append(L_tt_mean)

        elif x_axis == 'alg':
            for config_id in config_ids:
                err_traveltime[config_id] = []
                for alg_id in alg_ids:

                    L_tt_reps = []
                    for rep in replications:
                        L_traveltime = \
                            self.compute_traveltime_error_for_scenario_on_grid(rep, config_id,
                                                                               alg_id, grid, norm=norm)
                        L_tt_reps.append(L_traveltime)
                    # compute the mean
                    L_tt_mean = np.sum(L_tt_reps) / len(replications)
                    err_traveltime[config_id].append(L_tt_mean)

        # ==================================
        # convert travel time to mins
        for key in err_traveltime.keys():
            err_traveltime[key] = np.array(err_traveltime[key]) / 60.0

        # ==================================
        # visualize the travel time
        fig = plt.figure(figsize=figsize, dpi=100)
        ax = fig.add_axes([0.1, 0.15, 0.82, 0.75])

        # setup bar plot
        if plot_style == 'bar':
            loc = np.arange(len(xticklabels)).astype(float)
            width = 0.2
            loc_interpolation = loc - width * 1.5
            loc_filtering = loc - width * 0.5
            loc_enkf = loc + width * 0.5

        for key in err_traveltime.keys():

            if x_axis == 'config':
                if 'linear' in key:
                    label = 'Spatial'
                    marker = '^'
                    col = 'r'
                    if plot_style == 'bar':
                        loc = loc_interpolation
                        rect_interpolation = ax.bar(loc, err_traveltime[key], width, color=col, hatch='/')
                elif 'fisb' in key:
                    label = 'Spatio-temporal'
                    marker = 'o'
                    col = 'deepskyblue'
                    if plot_style == 'bar':
                        loc = loc_filtering
                        rect_filtering = ax.bar(loc, err_traveltime[key], width, color=col, hatch='*')
                elif 'enkf' in key:
                    label = 'EnKF'
                    col = 'g'
                    marker = '*'
                    if plot_style == 'bar':
                        loc = loc_enkf
                        rect_enkf = ax.bar(loc, err_traveltime[key], width, color=col, hatch='x')
            elif x_axis == 'alg':
                # label configurations, which mainly differ in sensors
                label = key.split('_')[-1]
            else:
                raise Exception('unrecognized x_axis')

            if plot_style == 'line':
                ax.plot(np.array(err_traveltime[key]), label=label, linewidth=2,
                        marker=marker, markersize=10, fillstyle='full', color=col)

        if norm == 'L1':
            option_str = 'MAE'
        elif norm == 'L2':
            option_str = 'RMSE'

        if title is None:
            ax.set_title('{0} norm of travel time error'.format(option_str), fontsize=fontsize[0])
        else:
            ax.set_title(title, fontsize=fontsize[0])

        ax.set_xlabel('{0}'.format(xlabel), fontsize=fontsize[1])
        if plot_style == 'line':
            ax.set_xlim([0, len(xticklabels) - 1])
        elif plot_style == 'bar':
            ax.set_xlim([-0.5, len(xticklabels) - 0.5])
        ax.set_ylim(ylim)
        x_ticks = np.arange(0, len(xticklabels))
        x_ticklabels = xticklabels

        ax.set_xticks(x_ticks)
        ax.set_xticklabels(x_ticklabels, fontsize=fontsize[2])
        ax.tick_params(axis='both', which='major', labelsize=fontsize[2])

        ax.set_ylabel('travel time error (min)', fontsize=fontsize[1])
        # ax.grid(True)
        if plot_style == 'line':
            plt.legend(prop={'size': fontsize[1]}, loc='best')
        elif plot_style == 'bar':
            ax.legend((rect_interpolation[0], rect_filtering[0], rect_enkf[0]),
                      ('Spatial', 'Spatio-temporal', 'EnKF'), prop={'size': fontsize[1]})

        if save_fig is True:
            plt.savefig('{0}'.format(save_fig_name), bbox_inches='tight')
            plt.clf()
            plt.close()
        else:
            plt.draw()

        # ================================================================================
        # plot the cost-effectiveness: how much more improvement in average by adding each sensor
        # ================================================================================
        if plot_cost_eff is True:
            fig = plt.figure(figsize=figsize, dpi=100)
            ax = fig.add_axes([0.1, 0.15, 0.82, 0.75])

            if cost_eff_num_sensors is None:
                print('Warning: Please specify cost_eff_num_sensors for cost-effectiveness plot')
                num_sensors = np.arange(0, len(xticklabels))
            else:
                num_sensors = np.array([float(i) for i in cost_eff_num_sensors])

            for key in err_traveltime.keys():

                if x_axis == 'config':
                    if 'linear' in key:
                        label = 'Spatial'
                        marker = '^'
                        col = 'r'
                    elif 'fisb' in key:
                        label = 'Spatio-temporal'
                        marker = 'o'
                        col = 'deepskyblue'
                    elif 'enkf' in key:
                        label = 'EnKF'
                        col = 'g'
                        marker = '*'
                else:
                    label = key.split('_')[-1]

                eff_per_sensor = []
                eff_per_sensor.append(100)
                for num in range(1, len(num_sensors)):
                    rel_improvment = -100.0 * (err_traveltime[key][num] -
                                               err_traveltime[key][num - 1]) / err_traveltime[key][num - 1]
                    eff_per_sensor.append(rel_improvment / (num_sensors[num] - num_sensors[num - 1]))

                ax.plot(eff_per_sensor, label=label, linewidth=2,
                        marker=marker, color=col, markersize=10)

            if norm == 'L1':
                option_str = 'MAE'
            elif norm == 'L2':
                option_str = 'RMSE'

            if cost_eff_title is None:
                ax.set_title('Improvement of Travel time {0}'.format(option_str),
                             fontsize=cost_eff_fontsize[0])
            else:
                ax.set_title(cost_eff_title, fontsize=cost_eff_fontsize[0])

            ax.set_xlabel('{0}'.format(xlabel), fontsize=cost_eff_fontsize[1])
            ax.set_xlim([0, len(xticklabels) - 1])
            x_ticks = np.arange(0, len(xticklabels))
            x_ticklabels = xticklabels

            ax.set_xticks(x_ticks)
            ax.set_xticklabels(x_ticklabels, fontsize=cost_eff_fontsize[2])
            ax.tick_params(axis='both', which='major', labelsize=cost_eff_fontsize[2])

            ax.set_ylabel('Improvement per sensor (%)', fontsize=cost_eff_fontsize[1])
            # ax.grid(True)

            plt.legend(prop={'size': cost_eff_fontsize[1]})

            if save_fig is True:
                plt.savefig('cost_eff_{0}'.format(save_fig_name), bbox_inches='tight')
                plt.clf()
                plt.close()
            else:
                plt.draw()

    # deprecated
    def sorted_bar_speed(self, rep, scenarios, norm='L1', area='all', queue_buffer=800,
                         grid=(5, 50),
                         xlabel='Configurations', xticklabels=None,
                         title='MAE error for speed',
                         unit='imperial', save_fig=False, save_fig_name='bar_speed.pdf'):
        """
        This function plots the sorted bar plot to compare across the algorithms and configurations (num_sensors).
        :param rep: int, replication id
        :param scenarios: a list of tuples (config_id, alg_id)
        :param norm: 'L1', 'L2'
        :param area: 'all', 'freeflow', 'congflow', 'aroundqueue', for speed MAE in corresponding area
        :param queue_buffer: meters, +/- the queue to define the areas
        :param grid: (s, m) the grid to compute against.
        :param xlabel: str, the x label
        :param xticklabels: list of str, the x tick labels
        :param title: str, the title
        :param unit: 'metric'(m/s) or 'imperial' (mph)
        :param save_fig: True of False
        :param save_fig_name: name of the saved figure
        :return: a plotted or saved figure
        """

        # ======================================================
        # compute and get the error for all configurations.
        err_speed = []
        for tup in scenarios:
            config_id, alg_id = tup

            L_speed = \
                self.compute_speed_error_for_scenario_on_grid(rep, config_id,
                                                              alg_id, grid,
                                                              norm=norm, area=area,
                                                              queue_buffer=queue_buffer)
            err_speed.append(L_speed)

        # ======================================================
        # sort the results
        zipped_sces = zip(err_speed, scenarios)
        zipped_sces.sort()
        sorted_err, sorted_sce = zip(*zipped_sces)
        sorted_err = np.array(list(sorted_err))

        if unit == 'imperial':
            unit_str = 'mph'
            sorted_err = self.__metric2imperial(sorted_err, 'speed')
        elif unit == 'metric':
            unit_str = 'm/s'
        else:
            raise Exception('Unrecognized unit')

        # ======================================================
        # visualize the speed
        width = 0.45
        fig = plt.figure(figsize=(15, 6), dpi=100)
        ax = fig.add_axes([0.1, 0.15, 0.82, 0.75])
        barlist = ax.bar(np.arange(len(sorted_err)), sorted_err, width=width)

        # update the color
        for i, tup in enumerate(sorted_sce):
            alg = tup[1]
            if 'linear' in alg:
                barlist[i].set_color('b')
            elif 'fisb' in alg:
                barlist[i].set_color('g')
            elif 'enkf' in alg:
                barlist[i].set_color('r')

        # add legend
        linear_alg = mpatches.Patch(color='blue', label='Interpolation')
        fisb_alg = mpatches.Patch(color='green', label='Kernel filtering')
        enkf_alg = mpatches.Patch(color='red', label='Kalman filter')
        ax.legend(handles=[linear_alg, fisb_alg, enkf_alg], loc=2, prop={'size': 20})

        # set x ticks
        ax.set_title(title, fontsize=32)

        if xticklabels is None:
            xticklabels = ['{0}'.format(i[0]).strip('configHOMO_').strip('_RTMS') for i in sorted_sce]
        ax.set_xticks(np.arange(0, len(sorted_err)) + width / 2)
        ax.set_xticklabels(xticklabels, fontsize=28)
        ax.set_xlabel('{0}'.format(xlabel), fontsize=30)
        ax.set_xlim([-0.5, len(sorted_err) + 0.5])

        # set y ticks
        ax.set_ylabel('Speed error ({0})'.format(unit_str), fontsize=28)
        ax.grid(True)

        if save_fig is True:
            plt.savefig('{0}'.format(save_fig_name))
            plt.clf()
            plt.close()
        else:
            plt.draw()

    # deprecated
    def sorted_bar_queue(self, rep, scenarios, norm='L1',
                         grid=(5, 50),
                         xlabel='Configurations', xticklabels=None,
                         title='MAE error for queue length',
                         unit='imperial', save_fig=False, save_fig_name='bar_queue.pdf'):
        """
        This function plots the sorted bar plot to compare across the algorithms and configurations (num_sensors).
        :param rep: int replication id
        :param scenarios: a list of tuples (config_id, alg_id)
        :param norm: 'L1', 'L2'
        :param grid: (s, m), the grid to compute against
        :param xlabel: x label
        :param xticklabels: list of str, the x tick labels
        :param title: str, the title
        :param unit: 'metric'(m) or 'imperial' (mile)
        :param save_fig: True of False
        :param save_fig_name: name of the saved figure
        :return: a plotted or saved figure
        """

        # ======================================================
        # compute and get the error for all configurations.
        err_queue = []
        for tup in scenarios:
            config_id, alg_id = tup

            L_queue = \
                self.compute_queuelength_error_for_scenario_on_grid(rep, config_id,
                                                                    alg_id, grid,
                                                                    norm=norm)
            err_queue.append(L_queue)

        # ======================================================
        # sort the results
        zipped_sces = zip(err_queue, scenarios)
        zipped_sces.sort()
        sorted_err, sorted_sce = zip(*zipped_sces)
        sorted_err = np.array(list(sorted_err))

        if unit == 'imperial':
            unit_str = 'mile'
            sorted_err = self.__metric2imperial(sorted_err, 'distance')
        elif unit == 'metric':
            unit_str = 'm'
        else:
            raise Exception('Unrecognized unit')

        # ======================================================
        # visualize the speed
        width = 0.45
        fig = plt.figure(figsize=(15, 6), dpi=100)
        ax = fig.add_axes([0.1, 0.15, 0.82, 0.75])
        barlist = ax.bar(np.arange(len(sorted_err)), sorted_err, width=width)

        # update the color
        for i, tup in enumerate(sorted_sce):
            alg = tup[1]
            if 'linear' in alg:
                barlist[i].set_color('b')
            elif 'fisb' in alg:
                barlist[i].set_color('g')
            elif 'enkf' in alg:
                barlist[i].set_color('r')

        # add legend
        linear_alg = mpatches.Patch(color='blue', label='Interpolation')
        fisb_alg = mpatches.Patch(color='green', label='Kernel filtering')
        enkf_alg = mpatches.Patch(color='red', label='Kalman filter')
        ax.legend(handles=[linear_alg, fisb_alg, enkf_alg], loc=2, prop={'size': 20})

        # set x ticks
        ax.set_title(title, fontsize=32)

        if xticklabels is None:
            xticklabels = ['{0}'.format(i[0]).strip('configHOMO_').strip('_RTMS') for i in sorted_sce]
        ax.set_xticks(np.arange(0, len(sorted_err)) + width / 2)
        ax.set_xticklabels(xticklabels, fontsize=28)
        ax.set_xlabel('{0}'.format(xlabel), fontsize=30)
        ax.set_xlim([-0.5, len(sorted_err) + 0.5])

        # set y ticks
        ax.set_ylabel('Queue length error ({0})'.format(unit_str), fontsize=28)
        ax.grid(True)

        if save_fig is True:
            plt.savefig('{0}'.format(save_fig_name))
            plt.clf()
            plt.close()
        else:
            plt.draw()

    # deprecated
    def sorted_bar_traveltime(self, rep, scenarios, norm='L1',
                              grid=(5, 50),
                              xlabel='Configurations', xticklabels=None,
                              title='MAE error for travel time',
                              unit='metric', save_fig=False, save_fig_name='bar_tt.pdf'):
        """
        This function plots the sorted bar plot to compare across the algorithms and configurations (num_sensors).
        :param rep: the replication id
        :param scenarios: a list of tuples (config_id, alg_id)
        :param norm: 'L1', 'L2'
        :param grid: (s, m), the grid to compute against
        :param xlabel: the x label
        :param xticklabels: list of str, the x tick labels
        :param title: str, the title
        :param unit: 'metric'(s) or 'imperial' (s)
        :param save_fig: True of False
        :param save_fig_name: name of the saved figure
        :return: a plotted or saved figure
        """

        # ======================================================
        # compute and get the error for all configurations.
        err_tt = []
        for tup in scenarios:
            config_id, alg_id = tup

            L_tt = \
                self.compute_traveltime_error_for_scenario_on_grid(rep, config_id,
                                                                   alg_id, grid,
                                                                   norm=norm)
            err_tt.append(L_tt)

        # ======================================================
        # sort the results
        zipped_sces = zip(err_tt, scenarios)
        zipped_sces.sort()
        sorted_err, sorted_sce = zip(*zipped_sces)
        sorted_err = np.array(list(sorted_err))

        unit_str = 's'

        # ======================================================
        # visualize the speed
        width = 0.45
        fig = plt.figure(figsize=(15, 6), dpi=100)
        ax = fig.add_axes([0.1, 0.15, 0.82, 0.75])
        barlist = ax.bar(np.arange(len(sorted_err)), sorted_err, width=width)

        # update the color
        for i, tup in enumerate(sorted_sce):
            alg = tup[1]
            if 'linear' in alg:
                barlist[i].set_color('b')
            elif 'fisb' in alg:
                barlist[i].set_color('g')
            elif 'enkf' in alg:
                barlist[i].set_color('r')

        # add legend
        linear_alg = mpatches.Patch(color='blue', label='Interpolation')
        fisb_alg = mpatches.Patch(color='green', label='Kernel filtering')
        enkf_alg = mpatches.Patch(color='red', label='Kalman filter')
        ax.legend(handles=[linear_alg, fisb_alg, enkf_alg], loc=2, prop={'size': 20})

        # set x ticks
        ax.set_title(title, fontsize=32)

        if xticklabels is None:
            xticklabels = ['{0}'.format(i[0]).strip('configHOMO_').strip('_RTMS') for i in sorted_sce]
        ax.set_xticks(np.arange(0, len(sorted_err)) + width / 2)
        ax.set_xticklabels(xticklabels, fontsize=28)
        ax.set_xlabel('{0}'.format(xlabel), fontsize=30)
        ax.set_xlim([-0.5, len(sorted_err) + 0.5])

        # set y ticks
        ax.set_ylabel('Travel time error ({0})'.format(unit_str), fontsize=28)
        ax.grid(True)

        if save_fig is True:
            plt.savefig('{0}'.format(save_fig_name))
            plt.clf()
            plt.close()
        else:
            plt.draw()

    def compute_speed_error_for_scenario_on_grid(self, rep, config_id, alg_id, grid,
                                                 norm='L1', area='all', queue_buffer=800,
                                                 max_speed=27.2):
        """
        This function returns the L1 norm error of the speed, queue, and travel time
        The scenario is defined by the combination of the config_id and the alg_id
        :param rep: the replication data
        :param config_id: the configuration id
        :param alg_id: the algorithm used
        :param grid: tuple (s, m), compute against the true state grid.
        :param norm: 'L1' or 'L2'
        :param area: 'all','freeflow','congflow', 'aroundqueue'
        :param queue_buffer: m, if queue is at 6000 m,
                then [6000-queue_buffer, 6000+queue_buffer] is regarded as queue area.
        :param max_speed: m/s, for cleaning up the speed data
        :return: float, average speed error in m/s
        """
        # get the true speed, queue, and travel time.
        true_speed_file = self.__truestate_result_dir[rep] + 'truestate_{0}s{1}m_speed.txt'.format(grid[0], grid[1])
        true_speed = np.genfromtxt(true_speed_file, delimiter=',')

        # get the estimated result from this scenario
        est_speed_file = self.__est_result_dir[rep] + '{0}_{1}_speed.txt'.format(config_id, alg_id)
        est_speed = np.genfromtxt(est_speed_file, delimiter=',')

        # --------------------------------------------
        # Cleaning the data:
        #   - nan of inf changed to max speed
        #   - 0<0.45 m/s (1 mph), change to 0.45 to prevent outlier caused by completely stopped traffic
        true_speed[np.isinf(true_speed) | np.isnan(true_speed)] = max_speed  # the vm in calibrated FD
        est_speed[np.isinf(est_speed) | np.isnan(est_speed)] = max_speed

        # ======================================================
        # If the dimensions are not the same, then need to expand to the same dimension
        if true_speed.shape != est_speed.shape:

            num_rows = np.max([true_speed.shape[0], est_speed.shape[0]])
            num_cols = np.max([true_speed.shape[1], est_speed.shape[1]])

            if true_speed.shape[0] != num_rows or true_speed.shape[1] != num_cols:
                true_speed = self.__expand_matrix(true_speed, (num_rows, num_cols))

            if est_speed.shape[0] != num_rows or est_speed.shape[1] != num_cols:
                est_speed = self.__expand_matrix(est_speed, (num_rows, num_cols))

        else:
            # rows: time, cols: space
            num_rows = true_speed.shape[0]
            num_cols = true_speed.shape[1]

        # ======================================================
        # now compute the L1 norm error. Normalized to the number of available cells.
        # compute the L1 norm of the speed
        err_speed_matrix = np.abs(est_speed - true_speed)
        notNaN = ~np.isnan(err_speed_matrix) & ~np.isinf(err_speed_matrix)
        err_speed_matrix[~notNaN] = 0

        start_time = 60  # 60*5 = 5 min

        # normalize to the number of cells in respective areas
        if area == 'all':
            # only use the states after 50 steps (60*5 = 300 s)
            active_err_speed_matrix = err_speed_matrix[start_time:, :]

            if norm == 'L1':
                L_speed_err = active_err_speed_matrix.sum() \
                              / notNaN[start_time:, :].sum()
            elif norm == 'L2':
                L_speed_err = np.sqrt(np.power(active_err_speed_matrix, 2).sum()
                                      / notNaN[start_time:, :].sum())
            else:
                raise Exception('Unrecognized norm')

        else:
            # compute the speed error in corresponding areas
            true_queue_file = self.__truestate_result_dir[rep] + \
                              'truestate_{0}s{1}m_queue.txt'.format(grid[0], grid[1])
            queue = np.array(np.genfromtxt(true_queue_file, delimiter=','))

            # convert the queue length into cells [0, queue_cell[t]] is freeflow
            queue_cell = num_cols - queue / grid[1]
            queue_cell = np.array([int(i) for i in queue_cell])
            buffer_cells = int(queue_buffer / grid[1])

            # append the averaged error across cell at each time here.
            avg_speed_err_each_step = []

            if area == 'freeflow':
                # upstream queue buffer; [0:us_queue_cell[t] are free flow cells]
                us_queue_cell = queue_cell - buffer_cells
                us_queue_cell[us_queue_cell < 0] = 0
                us_queue_cell[us_queue_cell >= num_cols] = num_cols - 1

                # skip first 5 min
                for t in range(start_time, num_rows):
                    active_cells = err_speed_matrix[t, 0:us_queue_cell[t]]
                    active_notNaN = notNaN[t, 0:us_queue_cell[t]]
                    num_active_meas = active_notNaN.sum()
                    if num_active_meas != 0:
                        if norm == 'L1':
                            avg_speed_err_each_step.append(active_cells.sum() / num_active_meas)
                        elif norm == 'L2':
                            avg_speed_err_each_step.append(np.sqrt(np.power(active_cells, 2).sum() /
                                                                   num_active_meas))
                        else:
                            raise Exception('Unrecognized norm')

            elif area == 'congflow':
                # downstream queue buffer; [ds_queue_cell[t], num_cols] are congested cells
                ds_queue_cell = queue_cell + buffer_cells
                ds_queue_cell[ds_queue_cell < 0] = 0
                ds_queue_cell[ds_queue_cell >= num_cols] = num_cols - 1

                # skip first 5 min
                for t in range(start_time, num_rows):
                    active_cells = err_speed_matrix[t, ds_queue_cell[t]:num_cols]
                    active_notNaN = notNaN[t, ds_queue_cell[t]:num_cols]
                    num_active_meas = active_notNaN.sum()

                    if num_active_meas != 0:
                        if norm == 'L1':
                            avg_speed_err_each_step.append(active_cells.sum() / num_active_meas)
                        elif norm == 'L2':
                            avg_speed_err_each_step.append(np.sqrt(np.power(active_cells, 2).sum() /
                                                                   num_active_meas))
                        else:
                            raise Exception('Unrecognized norm')

            elif area == 'aroundqueue':
                # [us_queue_cell[t], ds_queue_cell[t]] are around queue cells
                # upstream queue buffer;
                us_queue_cell = queue_cell - buffer_cells
                us_queue_cell[us_queue_cell < 0] = 0
                us_queue_cell[us_queue_cell >= num_cols] = num_cols - 1
                # downstream queue buffer; []
                ds_queue_cell = queue_cell + buffer_cells
                ds_queue_cell[ds_queue_cell < 0] = 0
                ds_queue_cell[ds_queue_cell >= num_cols] = num_cols - 1

                # skip first 5 min
                for t in range(start_time, num_rows):
                    active_cells = err_speed_matrix[t, us_queue_cell[t]:ds_queue_cell[t]]
                    active_notNaN = notNaN[t, us_queue_cell[t]:ds_queue_cell[t]]
                    num_active_meas = active_notNaN.sum()

                    if num_active_meas != active_cells.size:
                        print('num_active_meas is: {0}'.format(num_active_meas))
                        print('active_cells: {0}'.format(active_cells))

                    if num_active_meas != 0:
                        if norm == 'L1':
                            avg_speed_err_each_step.append(active_cells.sum()
                                                           / num_active_meas)
                        elif norm == 'L2':
                            avg_speed_err_each_step.append(np.sqrt(np.power(active_cells, 2).sum() /
                                                                   num_active_meas))
                        else:
                            raise Exception('Unrecognized norm')

            avg_speed_err_each_step = np.array(avg_speed_err_each_step)

            if norm == 'L1':
                L_speed_err = avg_speed_err_each_step.sum() / avg_speed_err_each_step.size
                # print('avg_speed_error_each_step: {0}'.format(avg_speed_err_each_step))
            elif norm == 'L2':
                L_speed_err = np.sqrt(np.power(avg_speed_err_each_step, 2).sum() /
                                      avg_speed_err_each_step.size)
            else:
                raise Exception('Unrecognized norm')

        return L_speed_err

    def compute_queuelength_error_for_scenario_on_grid(self, rep, config_id, alg_id, grid, norm='L1'):
        """
        This function returns the L1 norm error of the speed, queue, and travel time
        The scenario is defined by the combination of the config_id and the alg_id
        :param rep: the replication data
        :param config_id: the configuration id
        :param alg_id: the algorithm used
        :param grid: tuple (s, m), compute against the true state grid.
        :param norm: 'L1' or 'L2'
        :return: float, average queue length error in meters.
        """
        # get the queue length
        true_queue_file = self.__truestate_result_dir[rep] + 'truestate_{0}s{1}m_queue.txt'.format(grid[0], grid[1])
        true_queue = np.genfromtxt(true_queue_file, delimiter=',')

        # get the estimated result from this scenario
        est_queue_file = self.__est_result_dir[rep] + '{0}_{1}_queue.txt'.format(config_id, alg_id)
        est_queue = np.genfromtxt(est_queue_file, delimiter=',')

        # ======================================================
        # If the dimensions are not the same, then need to expand to the same dimension
        # expand the array to the same size
        if len(true_queue) > len(est_queue):
            est_queue = self.__expand_array(est_queue, len(true_queue))
        elif len(true_queue) < len(est_queue):
            true_queue = self.__expand_array(true_queue, len(est_queue))

        # ======================================================
        # now compute the L1 norm error. Normalized to the number of available cells.

        # compute the L1 norm error of the queue length
        err_queue = np.abs(est_queue - true_queue)
        notNaN = ~np.isnan(err_queue) & ~np.isinf(err_queue)
        err_queue[~notNaN] = 0

        # only use the states after 60 steps (60*5 = 300 s)
        start_time = 60
        active_err_queue_matrix = err_queue[start_time:]

        # normallize to the number of cells
        if norm == 'L1':
            L_queue_err = active_err_queue_matrix.sum() / notNaN[start_time:].sum()
        elif norm == 'L2':
            L_queue_err = np.sqrt(np.power(active_err_queue_matrix, 2).sum() /
                                  (notNaN[start_time:].sum()))
        else:
            raise Exception('unrecognized norm')

        return L_queue_err

    def compute_traveltime_error_for_scenario_on_grid(self, rep, config_id, alg_id, grid,
                                                      norm='L1'):
        """
        This function returns the L1 norm error of the speed, queue, and travel time
        The scenario is defined by the combination of the config_id and the alg_id
        :param rep: the replication data
        :param config_id: the configuration id
        :param alg_id: the algorithm used
        :param grid: tuple (s, m), compute against the true state grid.
        :param norm: 'L1' of 'L2'
        :return: float, average travel time error in seconds
        """
        # get the true speed, queue, and travel time.
        true_traveltime_file = self.__truestate_result_dir[rep] + 'truestate_{0}s{1}m_true_traveltime.txt'.format(
            grid[0], grid[1])

        true_traveltime = np.genfromtxt(true_traveltime_file, delimiter=',')

        # get the estimated result from this scenario
        est_traveltime_file = self.__est_result_dir[rep] + '{0}_{1}_traveltime.txt'.format(config_id, alg_id)

        est_traveltime = np.genfromtxt(est_traveltime_file, delimiter=',')

        # ======================================================
        # If the dimensions are not the same, then need to expand to the same dimension

        # expand the array to the same size
        if len(true_traveltime) > len(est_traveltime):
            est_traveltime = self.__expand_array(est_traveltime, len(true_traveltime))
        elif len(true_traveltime) < len(est_traveltime):
            true_traveltime = self.__expand_array(true_traveltime, len(est_traveltime))

        # ======================================================
        # now compute the L1 norm error. Normalized to the number of available cells.
        # compute the L1 norm error of the travel time
        err_traveltime = np.abs(est_traveltime - true_traveltime)

        notNaN = ~np.isnan(err_traveltime) & ~np.isinf(err_traveltime)
        err_traveltime[~notNaN] = 0

        # only use the states after 50 steps (50*5 = 250 s)
        start_time = 60
        active_err_traveltime = err_traveltime[start_time:]

        # normallize to the number of cells
        if norm == 'L1':
            L_traveltime_error = (active_err_traveltime.sum()) / (notNaN[start_time:].sum())
        elif norm == 'L2':
            L_traveltime_error = np.sqrt(np.power(active_err_traveltime, 2).sum() /
                                         notNaN[start_time].sum())
        else:
            raise Exception('Unrecognized norm')

        return L_traveltime_error

    def __expand_matrix(self, mat, shape):
        """
        This function expand the matrix mat to the size defined by shape
        :param mat: a 2d array/matrix
        :param shape: (num_row, num_col)
        :return: mat_expanded
        """

        mat_expanded = None
        # If only need to expand rows
        if mat.shape[0] != shape[0] and mat.shape[1] == shape[1]:
            mat_expanded = np.array([]).reshape(0, mat.shape[1])
            # only duplicate rows
            num_dup = shape[0] / mat.shape[0]
            for row in range(0, mat.shape[0]):
                for i in range(0, num_dup):
                    mat_expanded = np.vstack([mat_expanded, mat[row, :].reshape(1, mat.shape[1])])

        # if only need to expand columns
        elif mat.shape[0] == shape[0] and mat.shape[1] != shape[1]:
            mat_expanded = np.array([]).reshape(mat.shape[0], 0)

            # only duplicate cols
            num_dup = shape[1] / mat.shape[1]
            for col in range(0, mat.shape[1]):
                for i in range(0, num_dup):
                    # print('mat_expanded:{0}; mat[:,col].shape {1}'.format(mat_expanded.shape,
                    #                                                       mat[:, col].shape))
                    mat_expanded = np.hstack([mat_expanded, mat[:, col].reshape(mat.shape[0], 1)])

        # if both axis need to be expanded
        elif mat.shape[0] != shape[0] and mat.shape[1] == shape[1]:
            mat_tmp = np.array([]).reshape(0, mat.shape[1])
            # first duplicate rows
            num_dup = shape[0] / mat.shape[0]
            for row in range(0, mat.shape[0]):
                for i in range(0, num_dup):
                    mat_tmp = np.vstack([mat_tmp, mat[row, :].reshape(1, mat.shape[1])])

            mat_expanded = np.array([]).reshape(mat_tmp.shape[0], 0)
            # then duplicate cols
            num_dup = shape[1] / mat.shape[1]
            for col in range(0, mat.shape[1]):
                for i in range(0, num_dup):
                    mat_expanded = np.hstack([mat_expanded, mat_tmp[:, col].reshape(mat.shape[0], 1)])

        if mat_expanded is None:
            raise Exception('Unable to expand matrix')

        return mat_expanded

    def __expand_array(self, array, length):
        """
        This function expand the array to the size defined by length.
        :param array: a 1d array/matrix
        :param length: int
        :return: array_expanded
        """

        array_expanded = None
        if len(array) != length:

            num_dup = length / len(array)
            array_expanded = []

            for val in array:
                for i in range(0, num_dup):
                    array_expanded.append(val)

            array_expanded = np.array(array_expanded)

            return array_expanded
        else:
            return array

    # ================================================================================ #
    # Compute the queue length and travel time from estimated speed
    # ================================================================================ #
    def compute_queue_for_scenario(self, rep, config_id, alg_id, v_threshold,
                                   estimate_grid=(5, 50)):
        """
        This function is the standard methods for computing the queue length based on the speed field for all algorithms
        :param rep: int, the replication id
        :param config_id: str,
        :param alg_id:  str,
        :param v_threshold: m/s, the speed threshold under which will be regarded as congested
        :param estimate_grid: (s, m), the estimation grid for the speed for this configuration and algorithm
        :return: computed queue will be saved in the result folder.
        """

        # read speed estimates from the file
        # read the result from file
        est_speed_file = self.__est_result_dir[rep] + '{0}_{1}_speed.txt'.format(config_id, alg_id)

        speed_data = np.genfromtxt(est_speed_file, delimiter=',')
        # each row is a measurement at each cell
        speed_data = np.matrix(speed_data)

        # find the first index that is congested at each step
        len_cell = estimate_grid[1]
        num_cells = speed_data.shape[1]
        queue_length = []
        for step in range(0, speed_data.shape[0]):

            # if there are np.nan values
            try:
                congested = np.array(speed_data[step, :] <= v_threshold)[0]
            except RuntimeWarning:
                queue_length.append(np.nan)
                continue

            # starting form the last cell, stop at the cell that is not congested
            queue_cell = num_cells
            found_queue_end = False
            for cell in range(num_cells - 1, -1, -1):
                if congested[cell]:
                    queue_cell = cell
                else:
                    # found a uncongested cell, check if the next cell is also un congested
                    if ~congested[cell - 1]:
                        # if not congested, then found the location of the queue
                        queue_length.append(len_cell * (num_cells - queue_cell))
                        found_queue_end = True
                        break
                    else:
                        # otherwise, just noise in one cell, continue
                        queue_cell = cell

            # if queue end not found, then the entire road is congested
            if found_queue_end is False:
                queue_length.append(len_cell * num_cells)

        # ========================================================
        # a smoothing step is performed.
        # Assume queue can not propogate faster than 3 cells per step
        smooth_win = 10
        for i in range(1, len(queue_length)):
            if queue_length[i] - queue_length[i - 1] > 3 * len_cell:
                # grows faster than 3 cells per step, smooth based on the last 3 measure
                queue_length[i] = np.mean(queue_length[np.max([0, i - smooth_win]):i + 1])
            elif queue_length[i] - queue_length[i - 1] < -3 * len_cell:
                # dissipate faster then 3 cells per step, change to one cell
                queue_length[i] = np.mean(queue_length[np.max([0, i - smooth_win]):i + 1])
        # ========================================================

        # save the queue length in file.
        est_queue_file = self.__est_result_dir[rep] + '{0}_{1}_queue.txt'.format(config_id, alg_id)
        queue_length = np.array(queue_length)
        queue_length = queue_length.reshape((speed_data.shape[0], 1))
        np.savetxt(est_queue_file, queue_length, delimiter=",")

    def compute_traveltime_for_scenario(self, rep, config_id, alg_id, grid_res=(5, 50),
                                        max_speed=27.2, min_speed=1.34):
        """
        This function is the standard method for computing the instantaneous travel time for all algorithms
        :param rep: int, the replication id
        :param config_id: str, the configuration id
        :param alg_id: str, the algorithm id
        :param max_speed: (m/s) all speed will be truncated to [min_speed, max_speed]
        :param min_speed: (m/s)
        :return: will be saved in _traveltime.txt in the corresponding folder
        """
        # read speed estimates from the file
        # read the result from file
        est_speed_file = self.__est_result_dir[rep] + '{0}_{1}_speed.txt'.format(config_id, alg_id)

        speed_data = np.genfromtxt(est_speed_file, delimiter=',')
        # each row is a measurement at each cell
        speed_data = np.matrix(speed_data)

        # --------------------------------------------
        # Cleaning the data:
        #   - inf, or nan, all changed to max_speed
        #   - 0<min_speed m/s (1 mph), change to min_speed to prevent outlier caused by completely stopped traffic
        speed_data[np.isnan(speed_data) | np.isinf(speed_data)] = max_speed

        print('number of substituted stop speed cells: {0}'.format(sum(sum(np.array(speed_data < min_speed)))))
        speed_data[speed_data < min_speed] = min_speed
        # --------------------------------------------

        traveltime = []
        len_cell = grid_res[1]
        for step in range(0, speed_data.shape[0]):
            traveltime.append(np.sum(len_cell / speed_data[step, :]))

        # save the travel time in file
        est_traveltime_file = self.__est_result_dir[rep] + '{0}_{1}_traveltime.txt'.format(config_id, alg_id)
        traveltime = np.array(traveltime).reshape((speed_data.shape[0], 1))
        np.savetxt(est_traveltime_file, traveltime, delimiter=",")

    def compute_true_queue_for_grid(self, rep, grid_res, v_threshold):
        """
        This function computes the true queue based on the true speed field.
        :param rep: int, the replication id
        :param grid_res: (s, m) grid for the true speed field
        :param v_threshold: m/s, the speed threshold under which will be regarded as congested
        :return: computed queue will be saved in the result folder.
        """

        # read the result from file
        true_speed_file = self.__truestate_result_dir[rep] + \
                          'truestate_{0}s{1}m_speed.txt'.format(grid_res[0], grid_res[1])

        speed_data = np.genfromtxt(true_speed_file, delimiter=',')
        speed_data = np.array(speed_data)

        # find the first index that is congested at each step
        len_cell = grid_res[1]
        num_cells = speed_data.shape[1]
        queue_length = []
        for step in range(0, speed_data.shape[0]):

            # if there are np.nan values
            try:
                congested = np.array(speed_data[step, :] <= v_threshold)
            except RuntimeWarning:
                queue_length.append(np.nan)
                continue
            # print('congested: {0}'.format(congested))
            # starting form the last cell, stop at the cell that is not congested
            queue_cell = num_cells
            found_queue_end = False
            # print('num_cells: {0}'.format(num_cells))
            for cell in range(num_cells - 1, -1, -1):

                if congested[cell]:
                    queue_cell = cell
                else:

                    # found an un congested cell, check if the next cell is also un congested
                    if ~congested[cell - 1]:
                        # if not congested, then found the location of the queue
                        queue_length.append(len_cell * (num_cells - queue_cell))
                        found_queue_end = True
                        break
                    else:
                        # otherwise, just noise in one cell, continue
                        queue_cell = cell

            # if queue end not found, then the entire road is congested
            if found_queue_end is False:
                queue_length.append(len_cell * num_cells)

        # save the queue length in file.
        true_queue_file = self.__truestate_result_dir[rep] + \
                          'truestate_{0}s{1}m_queue.txt'.format(grid_res[0], grid_res[1])
        queue_length = np.array(queue_length).reshape((speed_data.shape[0], 1))
        np.savetxt(true_queue_file, queue_length, delimiter=",")

    def compute_trueinst_traveltime_for_grid(self, rep, grid_res, max_speed=27.2, min_speed=1.34):
        """
        This function computes the instantaneous travel time from the true speed field.
        :param rep: int, the replication id
        :param grid_res: (s, m), the grid resolution for the true speed field.
        :return: will be saved in _trueinst_traveltime.txt
        """
        # read the result from file
        true_speed_file = self.__truestate_result_dir[rep] + \
                          'truestate_{0}s{1}m_speed.txt'.format(grid_res[0], grid_res[1])

        speed_data = np.genfromtxt(true_speed_file, delimiter=',')
        speed_data = np.array(speed_data)

        # --------------------------------------------
        # Cleaning the data:
        #   - nan, inf, changed to max speed
        #   - 0<min_speed m/s (1 mph), change to min_speed to prevent outlier caused by completely stopped traffic
        speed_data[np.isnan(speed_data) | np.isinf(speed_data)] = max_speed  # the vm in calibrated FD

        print('number of substituted stop speed cells: {0}'.format(sum(sum(np.array(speed_data < min_speed)))))
        speed_data[speed_data < min_speed] = min_speed
        # --------------------------------------------

        traveltime = []
        len_cell = grid_res[1]
        for step in range(0, speed_data.shape[0]):
            traveltime.append(np.sum(len_cell / speed_data[step, :]))

        # save the travel time in file
        est_traveltime_file = self.__truestate_result_dir[rep] + \
                              'truestate_{0}s{1}m_trueinst_traveltime.txt'.format(grid_res[0], grid_res[1])
        traveltime = np.array(traveltime).reshape((speed_data.shape[0], 1))
        np.savetxt(est_traveltime_file, traveltime, delimiter=",")

    # ================================================================================ #
    # visualization function for the estimation results
    # ================================================================================ #

    def plot_speed_for_scenario(self, rep, config_id, alg_id, unit='metric', limit=(0, 40),
                                fig_size=(16, 8), fontsize=(40, 36, 34), grid_res=(5, 50),
                                save_fig=False, title=None):
        """
        This function plots the estimated speed profile over the entire time space domain in the specified unit. It will
        first expand the estimated speed field to the true states grid, and then visualize.
        :param rep: int, the replication ide
        :param config_id: string, the configuration id
        :param alg_id: string, the algorithm id
        :param unit: 'metric', 'imperial'; respectively 'm, s, m/s', and 'mile, hour, mph'
        :param limit: The limit of the colorbar in above units
        :param fig_size: tuple, the figure size
        :param fontsize: tuple of three: title, label and legend, tick label size
        :param grid_res: (s, m), the grid for the true states.
        :param save_fig: True or False, save figure or not
        :param title: the title of the figure
        :return: A figure profile with x-axis being the time, and y-axis being the space. (flow direction upwards)
        """
        # read the result from file
        est_speed_file = self.__est_result_dir[rep] + '{0}_{1}_speed.txt'.format(config_id, alg_id)
        true_speed_file = self.__truestate_result_dir[rep] + \
                          'truestate_{0}s{1}m_speed.txt'.format(self.grid_res[0], self.grid_res[1])

        est_speed = np.genfromtxt(est_speed_file, delimiter=',')
        est_speed = np.flipud(np.matrix(est_speed).T)

        true_speed = np.genfromtxt(true_speed_file, delimiter=',')
        true_speed = np.flipud(np.matrix(true_speed).T)

        # ======================================================
        # If the dimensions are not the same, then need to expand to the same dimension
        if true_speed.shape != est_speed.shape:

            num_rows = np.max([true_speed.shape[0], est_speed.shape[0]])
            num_cols = np.max([true_speed.shape[1], est_speed.shape[1]])

            if true_speed.shape[0] != num_rows or true_speed.shape[1] != num_cols:
                true_speed = self.__expand_matrix(true_speed, (num_rows, num_cols))

            if est_speed.shape[0] != num_rows or est_speed.shape[1] != num_cols:
                est_speed = self.__expand_matrix(est_speed, (num_rows, num_cols))

        speed_data = est_speed
        # ======================================================

        if unit == 'metric':
            # all internal values are in metric, so plot directly
            speed = np.flipud(speed_data)
            unit_str = 'm/s'
        elif unit == 'imperial':
            speed = self.__metric2imperial(speed_data, 'speed')
            unit_str = 'mph'
        else:
            raise Exception('Error: Unrecognized unit for plotting speed.')

        # truncate the first 5 mins
        start_time = int(5 * 60 / grid_res[0])

        fig = plt.figure(figsize=fig_size, dpi=100)
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        im = ax.imshow(speed[:, start_time:], cmap=plt.get_cmap('jet_r'),
                       interpolation='nearest',
                       aspect='auto',
                       vmin=limit[0], vmax=limit[1])
        ax.autoscale(False)

        if title is None:
            ax.set_title('Velocity estimate for {0} {1} ({2})'.format(alg_id,
                                                                      config_id,
                                                                      unit_str), fontsize=fontsize[0])
        else:
            ax.set_title(title, fontsize=fontsize[0])

        plt.xlabel('Time (min)', fontsize=fontsize[1])
        plt.ylabel('Space (mile)', fontsize=fontsize[1])

        # set x and y ticks
        x_ticks = np.array([start_time, 360, 720, 1080, 1440, 1800])
        x_ticklabels = 5 * x_ticks / 60  # in min
        x_ticks -= start_time

        if self.workzone == 'I80':
            y_ticks = np.linspace(0, 160, 6).astype(int)
        elif self.workzone == 'I57':
            y_ticks = np.linspace(0, 128, 5).astype(int)
        # y_ticks = y_ticks[:-1]      # remove the last one
        y_ticklabels = y_ticks[::-1] / 32

        ax.set_xticks(x_ticks)
        ax.set_xticklabels(x_ticklabels, fontsize=fontsize[2])
        ax.set_yticks(y_ticks)
        ax.set_yticklabels(y_ticklabels, fontsize=fontsize[2])
        ax.tick_params(axis='both', which='major', labelsize=fontsize[2])

        cax = fig.add_axes([0.95, 0.15, 0.01, 0.7])
        cbar = fig.colorbar(im, cax=cax, orientation='vertical')
        cbar.set_ticks([0, 20, 40, 60, 80])
        cbar.set_ticklabels([0, 20, 40, 60, 80])
        cax.tick_params(labelsize=30)

        if save_fig is True:
            plt.savefig('{0}_{1}_speed_est.pdf'.format(config_id, alg_id), bbox_inches='tight')
            # clear the memory once generated
            plt.clf()
            plt.close()
        else:
            plt.draw()

    # deprecated
    def plot_speed_est_error_for_scenario(self, rep, config_id, alg_id, unit='metric', limit=(0, 40),
                                          save_fig=False, title=None):
        """
        This function plots the estimated speed filed error.
        :param rep: int, the replication id
        :param config_id: str, the configuraiton id
        :param alg_id: the algorithm id
        :param unit: 'metric', 'imperial'; respectively 'm, s, m/s', and 'mile, hour, mph'
        :param limit: The limit of the color bar in above units
        :param save_fig: True or False
        :param title: figure title
        :return: A figure profile with x-axis being the time, and y-axis being the space. (flow direction upwards)
        """
        # read the estimation result from file
        est_speed_file = self.__est_result_dir[rep] + '{0}_{1}_speed.txt'.format(config_id, alg_id)
        true_speed_file = self.__truestate_result_dir[rep] + \
                          'truestate_{0}s{1}m_speed.txt'.format(self.grid_res[0], self.grid_res[1])

        est_speed = np.genfromtxt(est_speed_file, delimiter=',')
        est_speed = np.flipud(np.matrix(est_speed).T)

        true_speed = np.genfromtxt(true_speed_file, delimiter=',')
        true_speed = np.flipud(np.matrix(true_speed).T)

        # ======================================================
        # If the dimensions are not the same, then need to expand to the same dimension
        if true_speed.shape != est_speed.shape:

            num_rows = np.max([true_speed.shape[0], est_speed.shape[0]])
            num_cols = np.max([true_speed.shape[1], est_speed.shape[1]])

            if true_speed.shape[0] != num_rows or true_speed.shape[1] != num_cols:
                true_speed = self.__expand_matrix(true_speed, (num_rows, num_cols))

            if est_speed.shape[0] != num_rows or est_speed.shape[1] != num_cols:
                est_speed = self.__expand_matrix(est_speed, (num_rows, num_cols))

        # ======================================================
        if unit == 'metric':
            # all internal values are in metric, so plot directly
            unit_str = 'm/s'
        elif unit == 'imperial':
            est_speed = self.__metric2imperial(est_speed, 'speed')
            true_speed = self.__metric2imperial(true_speed, 'speed')
            unit_str = 'mph'
        else:
            raise Exception('Error: Unrecognized unit for plotting speed.')

        # absolute speed difference
        diff_speed = np.abs(est_speed - true_speed)

        fig = plt.figure(figsize=(15, 8), dpi=100)
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        im = ax.imshow(diff_speed, cmap=plt.get_cmap('jet'),
                       interpolation='nearest',
                       aspect='auto',
                       vmin=limit[0], vmax=limit[1])
        ax.autoscale(False)

        # plot flow sensors
        flow_sensor_ytick = {}
        flow_sensor_ytick['tick'] = []
        flow_sensor_ytick['label'] = []

        if config_id in self.all_config.keys():

            for s_id in self.all_config[config_id]['sensors'].keys():
                # the flow sensor will be drawn as a line in the beginning of the cells.
                cell_id = self.all_config[config_id]['sensors'][s_id]['cell']

                val = (len(self.space_grid) - 1) - 0.5 - cell_id
                flow_sensor_ytick['tick'].append(val)
                flow_sensor_ytick['label'].append('{0}'.format(s_id))
                ax.plot([-0.5, (len(self.time_grid) - 1) - 0.5], [val, val], color='k', linewidth=1)

            # set y ticks
            ax.set_yticks(flow_sensor_ytick['tick'])
            ax.set_yticklabels(flow_sensor_ytick['label'])
            # ax.set_yticklabels(['SB1', 'SB2', 'SB3', 'SB4', 'SB5', 'SB6', 'SB7', 'SB8', 'SB9'], fontsize=28)
        else:
            print('Warning: did not find configuration {0} in self.all_config'.format(config_id))

        if title is None:
            ax.set_title('Speed estimates error {0} {1} ({2})'.format(alg_id,
                                                                      config_id,
                                                                      unit_str), fontsize=24)
        else:
            ax.set_title(title, fontsize=24)

        plt.xlabel('Time', fontsize=30)
        plt.ylabel('Traffic direction $\mathbf{\Rightarrow}$', fontsize=20)
        cax = fig.add_axes([0.95, 0.25, 0.01, 0.5])
        cbar = fig.colorbar(im, cax=cax, orientation='vertical')

        cbar.ax.tick_params(labelsize=20)

        if save_fig is True:
            plt.savefig('{0}_{1}_speed_error.pdf'.format(config_id, alg_id))
            # clear the memory once generated
            plt.clf()
            plt.close()
        else:
            plt.draw()

    # deprecated
    def plot_density_for_scenario(self, rep, config_id, alg_id, unit='metric', limit=(0, 5),
                                  save_fig=False, title=None):
        """
        This function plots the estimated speed profile over the entire time space domain in the specified unit
        :param rep: int, the replication id
        :param config_id: str, the configuraiton id
        :param alg_id: the algorithm id
        :param unit: 'metric', 'imperial'; respectively 'm, s, m/s', and 'mile, hour, mph'
        :param limit: The limit of the color bar in above units
        :param save_fig: True or False
        :param title: figure title
        :return: A figure profile with x-axis being the time, and y-axis being the space. (flow direction upwards)
        """
        # read the result from file
        est_density_file = self.__est_result_dir[rep] + '{0}_{1}_density.txt'.format(config_id, alg_id)

        density_data = np.genfromtxt(est_density_file, delimiter=',')
        density_data = np.matrix(density_data).T

        if unit == 'metric':
            # all internal values are in metric, so plot directly
            density = np.flipud(density_data)
            unit_str = 'veh/m'
        elif unit == 'imperial':
            density = np.flipud(self.__metric2imperial(density_data, 'density'))
            unit_str = 'veh/mile'
        else:
            raise Exception('Error: Unrecognized unit for plotting speed.')

        fig = plt.figure(figsize=(15, 8), dpi=100)
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        im = ax.imshow(density, cmap=plt.get_cmap('jet'),
                       interpolation='nearest',
                       aspect='auto',
                       vmin=limit[0], vmax=limit[1])
        ax.autoscale(False)
        if title is None:
            ax.set_title('Density estimates ({0}_{1}) ({2})'.format(alg_id, config_id,
                                                                    unit_str), fontsize=24)
        else:
            ax.set_title(title, fontsize=24)

        plt.xlabel('Time')
        plt.ylabel('Space, traffic direction $\Rightarrow$')
        cax = fig.add_axes([0.95, 0.25, 0.01, 0.5])
        fig.colorbar(im, cax=cax, orientation='vertical')

        if save_fig is True:
            plt.savefig('{0}_{1}_density_est.pdf'.format(config_id, alg_id))
            # clear the memory once generated
            plt.clf()
            plt.close()
        else:
            plt.draw()

    def plot_queue_for_scenario(self, rep, config_id, alg_id, unit='metric', ylim=(0, 6600),
                                save_fig=False, title=None, true_grid=None,
                                fontsize=(40, 36, 34), fig_size=(10, 10)):
        """
        This function plots the length of the queue
        :param rep: int, the replication ide
        :param config_id: string, the configuration id
        :param alg_id: string, the algorithm id
        :param unit: 'metric', 'imperial'; respectively 'm, s, m/s', and 'mile, hour, mph'
        :param ylim: The limit of the color bar in above units
        :param fig_size: tuple, the figure size
        :param fontsize: tuple of three: title, label and legend, ticklabel size
        :param true_grid: (s, m), the grid for the true states.
        :param save_fig: True or False, save figure or not
        :param title: the title of the figure
        :return: a plotted or saved figure
        """
        # read the result from file
        est_queue_file = self.__est_result_dir[rep] + '{0}_{1}_queue.txt'.format(config_id, alg_id)
        est_queue = np.genfromtxt(est_queue_file, delimiter=',')

        if true_grid is None:
            true_queue_file = self.__truestate_result_dir[rep] + 'truestate_{0}s{1}m_queue.txt'.format(self.grid_res[0],
                                                                                                       self.grid_res[1])
        else:
            true_queue_file = self.__truestate_result_dir[rep] + 'truestate_{0}s{1}m_queue.txt'.format(true_grid[0],
                                                                                                       true_grid[1])

        true_queue = np.genfromtxt(true_queue_file, delimiter=',')

        if unit == 'imperial':
            est_queue = self.__metric2imperial(est_queue, 'distance')
            true_queue = self.__metric2imperial(true_queue, 'distance')
            unit_str = 'mile'
        else:
            unit_str = 'm'

        fig = plt.figure(figsize=fig_size, dpi=100)
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

        start_time = 60
        ax.plot(est_queue[start_time:], label='Estimated', linewidth=2)
        ax.plot(true_queue[start_time:], label='True', linewidth=2)

        ax.set_ylim(ylim)
        ax.set_xlabel('Time (min)', fontsize=fontsize[1])
        ax.set_ylabel('Queue length (miles)', fontsize=fontsize[1])

        # #==============================================================
        # set the x and y ticks
        # #==============================================================
        x_ticks = np.array([start_time, 360, 720, 1080, 1440, 1800])
        x_ticklabels = 5 * x_ticks / 60  # in min
        x_ticks = x_ticks - start_time
        ax.set_xlim([x_ticks[0], x_ticks[-1]])
        ax.set_xticks(x_ticks)
        ax.set_xticklabels(x_ticklabels, fontsize=fontsize[2])
        ax.tick_params(axis='both', which='major', labelsize=fontsize[2])

        if title is None:
            ax.set_title('Estimated queue length {0} {1} ({2})'.format(alg_id,
                                                                       config_id,
                                                                       unit_str),
                         fontsize=fontsize[0])
        else:
            ax.set_title(title, fontsize=fontsize[0])

        ax.grid(False)
        plt.legend(prop={'size': fontsize[1]})
        if save_fig is True:
            plt.savefig('{0}_{1}_queue_est.pdf'.format(config_id, alg_id), bbox_inches='tight')
            # clear the memory once generated
            plt.clf()
            plt.close()
        else:
            plt.draw()

    def plot_traveltime_for_scenario(self, rep, config_id, alg_id, true_grid=(5, 50),
                                     unit='metric', ylim=(0, 50),
                                     fig_size=(15, 8), fontsize=(40, 36, 34),
                                     save_fig=False, title=None):
        """
        This function plots the length of the queue
        :param rep: int, the replication ide
        :param config_id: string, the configuration id
        :param alg_id: string, the algorithm id
        :param unit: 'metric', 'imperial'; respectively 'm, s, m/s', and 'mile, hour, mph'
        :param ylim: The y axis limit
        :param fig_size: tuple, the figure size
        :param fontsize: tuple of three: title, label and legend, ticklabel size
        :param grid: (s, m), the grid for the true states.
        :param save_fig: True or False, save figure or not
        :param title: the title of the figure
        :return:
        """
        # read the result from file
        est_traveltime_file = self.__est_result_dir[rep] + '{0}_{1}_traveltime.txt'.format(config_id, alg_id)
        est_traveltime = np.genfromtxt(est_traveltime_file, delimiter=',')

        true_traveltime_file = self.__truestate_result_dir[rep] + \
                               'truestate_{0}s{1}m_true_traveltime.txt'.format(true_grid[0], true_grid[1])
        true_traveltime = np.genfromtxt(true_traveltime_file, delimiter=',')

        # smooth the true states
        window_width = 6
        smoothed_true_traveltime = deepcopy(true_traveltime)
        for i in range(0, len(true_traveltime)):
            if i < window_width:
                smoothed_true_traveltime[i] = np.mean(true_traveltime[0:i + 1])
            else:
                smoothed_true_traveltime[i] = np.mean(true_traveltime[i - window_width + 1:i + 1])

        # #==============================================================
        # visualization
        fig = plt.figure(figsize=fig_size, dpi=100)
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

        start_time = 60
        ax.plot(est_traveltime[start_time:] / 60.0, label='Estimated', linewidth=2)
        ax.plot(smoothed_true_traveltime[start_time:] / 60.0, label='True', linewidth=2)

        # set the x ticks
        # #==============================================================
        x_ticks = np.array([start_time, 360, 720, 1080, 1440, 1800])
        x_ticklabels = 5 * x_ticks / 60  # in min
        x_ticks -= start_time

        ax.set_xlabel('Time (min)', fontsize=fontsize[1])
        ax.set_ylabel('Travel time (min)', fontsize=fontsize[1])

        ax.set_xticks(x_ticks)
        ax.set_xticklabels(x_ticklabels, fontsize=fontsize[2])
        ax.set_xlim([x_ticks[0], x_ticks[-1]])
        ax.set_ylim(ylim)
        ax.tick_params(axis='both', which='major', labelsize=fontsize[2])

        if title is None:
            ax.set_title('Estimated travel time {0} {1}'.format(alg_id, config_id),
                         fontsize=fontsize[0])
        else:
            ax.set_title(title, fontsize=fontsize[0])

        ax.grid(False)

        plt.legend(prop={'size': fontsize[1]}, loc=2)

        if save_fig is True:
            plt.savefig('{0}_{1}_tt_est.pdf'.format(config_id, alg_id), bbox_inches='tight')
            # clear the memory once generated
            plt.clf()
            plt.close()
        else:
            plt.draw()

    # ================================================================================ #
    # visualizing true states
    # ================================================================================ #
    def plot_true_speed_field_partition(self, grid_res, rep, unit='metric',
                                        limit=(0, 40), queue_buffer=800,
                                        fig_size=(15, 8), fontsize=(40, 36, 34),
                                        save_fig=False, title=None):
        """
        This function plots the true speed profile in the specified unit with partition of the areas: freeflow,
            congested flow, and around the queue
        :param grid_res: (s, m), the grid for the true speed field.
        :param rep: the replication id
        :param unit: 'metric', 'imperial'; respectively 'm, s, m/s', and 'mile, hour, mph'
        :param limit: The limit of the colorbar in above units
        :param queue_buffer: the distance upstream and downstream of the queue, which is regarded as the area around
                the queue
        :param fig_size: tuple, the figure size
        :param fontsize: a tuple of three int, font size of : title, labels and legned, tick labels
        :param save_fig: True of False
        :param title: the figure tile
        :return: A figure profile with x-axis being the time, and y-axis being the space. (flow direction upwards)
        """
        # get the true speed from file
        true_speed_file = self.__truestate_result_dir[rep] + \
                          'truestate_{0}s{1}m_speed.txt'.format(grid_res[0], grid_res[1])

        speed_data = np.genfromtxt(true_speed_file, delimiter=',')
        speed_data = np.matrix(speed_data).T

        true_queue_file = self.__truestate_result_dir[rep] + \
                          'truestate_{0}s{1}m_queue.txt'.format(grid_res[0], grid_res[1])

        # ==========================================
        # get the true queue length from file
        queue = np.array(np.genfromtxt(true_queue_file, delimiter=','))

        if unit == 'metric':
            # all internal values are in metric, so plot directly
            speed = np.flipud(speed_data)
            unit_str = 'm/s'
        elif unit == 'imperial':
            speed = np.flipud(self.__metric2imperial(speed_data, 'speed'))
            unit_str = 'mph'
        else:
            raise Exception('Error: Unrecognized unit for plotting speed.')

        # convert the queue length into cells
        queue_cell = queue / grid_res[1]
        queue_cell = np.array([int(i) for i in queue_cell])
        buffer_cells = int(queue_buffer / grid_res[1])

        # upstream queue buffer
        us_queue_cell = queue_cell - buffer_cells
        us_queue_cell[us_queue_cell < 0] = 0
        us_queue_cell[us_queue_cell > speed_data.shape[0]] = speed_data.shape[0]

        # downstream queue buffer
        ds_queue_cell = queue_cell + buffer_cells
        ds_queue_cell[ds_queue_cell < 0] = 0
        ds_queue_cell[ds_queue_cell > speed_data.shape[0]] = speed_data.shape[0]

        fig = plt.figure(figsize=fig_size, dpi=100)
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        im = ax.imshow(speed, cmap=plt.get_cmap('jet_r'),
                       interpolation='nearest',
                       aspect='auto',
                       vmin=limit[0], vmax=limit[1])
        ax.autoscale(False)

        if title is None:
            ax.set_title('True speed ({0}) for Rep {1} on grid:{2}s{3}m'.format(unit_str, rep,
                                                                                grid_res[0], grid_res[1]),
                         fontsize=fontsize[0])
        else:
            ax.set_title(title, fontsize=fontsize[0])

        plt.xlabel('Time (s)')
        plt.ylabel('Space, traffic direction $\mathbf{\Rightarrow}$')
        cax = fig.add_axes([0.95, 0.25, 0.01, 0.5])
        fig.colorbar(im, cax=cax, orientation='vertical')

        # ==========================================
        # plot the queue buffer
        for t in range(0, speed.shape[1]):
            # The queue position
            ax.plot([t - 0.5, t + 0.5], [queue_cell[t] - 0.5, queue_cell[t] - 0.5], color='k', linewidth=2)

            # The upper bound of queue area
            ax.plot([t - 0.5, t + 0.5], [us_queue_cell[t] - 0.5, us_queue_cell[t] - 0.5], color='k', linewidth=2)

            # The lower bound of queue area
            ax.plot([t - 0.5, t + 0.5], [ds_queue_cell[t] - 0.5, ds_queue_cell[t] - 0.5], color='k', linewidth=2)

        if save_fig is True:
            plt.savefig('truestate_{0}s{1}m_speed.pdf'.format(grid_res[0], grid_res[1]))
            # clear the memory once generated
            plt.clf()
            plt.close()
        else:
            plt.draw()

    def plot_true_speed_for_rep(self, grid_res, rep, unit='metric', limit=(0, 40),
                                fig_size=(16, 8), fontsize=(40, 36, 34),
                                save_fig=False, title=None):
        """
        This function plots the true speed profile in the specified unit
        :param grid_res: (s, m), the grid for the true speed field.
        :param rep: the replication id
        :param unit: 'metric', 'imperial'; respectively 'm, s, m/s', and 'mile, hour, mph'
        :param limit: The limit of the colorbar in above units
        :param fig_size: tuple, the figure size
        :param fontsize: a tuple of three int, font size of : title, labels and legned, tick labels
        :param save_fig: True of False
        :param title: the figure tile
        :return: A figure profile with x-axis being the time, and y-axis being the space. (flow direction upwards)
        """
        # read the result from file
        true_speed_file = self.__truestate_result_dir[rep] + 'truestate_{0}s{1}m_speed.txt'.format(grid_res[0],
                                                                                                   grid_res[1])

        speed_data = np.genfromtxt(true_speed_file, delimiter=',')
        speed_data = np.matrix(speed_data).T

        if unit == 'metric':
            # all internal values are in metric, so plot directly
            speed = np.flipud(speed_data)
            unit_str = 'm/s'
        elif unit == 'imperial':
            speed = np.flipud(self.__metric2imperial(speed_data, 'speed'))
            unit_str = 'mph'
        else:
            raise Exception('Error: Unrecognized unit for plotting speed.')

        # trucate the first 5 mins
        start_time = int(5 * 60 / grid_res[0])

        fig = plt.figure(figsize=fig_size, dpi=100)
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        im = ax.imshow(speed[:, start_time:], cmap=plt.get_cmap('jet_r'),
                       interpolation='nearest',
                       aspect='auto',
                       vmin=limit[0], vmax=limit[1])
        ax.autoscale(False)

        if title is None:
            ax.set_title('True speed ({0}) for Rep {1} on grid:{2}s{3}m'.format(unit_str, rep,
                                                                                grid_res[0], grid_res[1]),
                         fontsize=fontsize[0])
        else:
            ax.set_title(title, fontsize=fontsize[0])

        plt.xlabel('Time (min)', fontsize=fontsize[1])
        plt.ylabel('Space (mile)', fontsize=fontsize[1])

        # =============================================
        # set x and y ticks
        x_ticks = np.array([start_time, 360, 720, 1080, 1440, 1800])
        x_ticklabels = 5 * x_ticks / 60  # in min
        x_ticks = x_ticks - start_time

        if self.workzone == 'I80':
            y_ticks = np.linspace(0, 160, 6).astype(int)
        else:
            y_ticks = np.linspace(0, 128, 5).astype(int)
        # y_ticks = y_ticks[:-1]      # remove the last one
        y_ticklabels = y_ticks[::-1] / 32

        ax.set_xticks(x_ticks)
        ax.set_xticklabels(x_ticklabels, fontsize=fontsize[2])
        ax.set_yticks(y_ticks)
        ax.set_yticklabels(y_ticklabels, fontsize=fontsize[2])
        ax.tick_params(axis='both', which='major', labelsize=fontsize[2])

        cax = fig.add_axes([0.95, 0.15, 0.01, 0.7])
        cbar = fig.colorbar(im, cax=cax, orientation='vertical')
        cbar.set_ticks([0, 20, 40, 60, 80])
        cbar.set_ticklabels([0, 20, 40, 60, 80])
        cax.tick_params(labelsize=30)

        if save_fig is True:
            plt.savefig('truestate_{0}s{1}m_speed.pdf'.format(grid_res[0], grid_res[1]),
                        bbox_inches='tight')
            # clear the memory once generated
            plt.clf()
            plt.close()
        else:
            plt.draw()

    # deprecated
    def plot_true_density_for_rep(self, grid_res, rep, unit='metric', limit=(0, 10),
                                  save_fig=False, title=None):
        """
        This function plots the true density profile in the specified unit
        :param grid_res: (s, m), the grid for the true speed field.
        :param rep: the replication id
        :param unit: 'metric', 'imperial'; respectively 'm, s, m/s', and 'mile, hour, mph'
        :param limit: The limit of the colorbar in above units
        :param save_fig: True of False
        :param title: the figure tile
        :return: A figure profile with x-axis being the time, and y-axis being the space. (flow direction upwards)
        """
        # read the result from file
        true_density_file = self.__truestate_result_dir[rep] + 'truestate_{0}s{1}m_density.txt'.format(grid_res[0],
                                                                                                       grid_res[1])

        density_data = np.genfromtxt(true_density_file, delimiter=',')
        density_data = np.matrix(density_data).T

        if unit == 'metric':
            # all internal values are in metric, so plot directly
            density = np.flipud(density_data)
            unit_str = 'veh/m'
        elif unit == 'imperial':
            density = np.flipud(self.__metric2imperial(density_data, 'density'))
            unit_str = 'veh/mile'
        else:
            raise Exception('Error: Unrecognized unit for plotting speed.')

        fig = plt.figure(figsize=(15, 8), dpi=100)
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        im = ax.imshow(density, cmap=plt.get_cmap('jet'),
                       interpolation='nearest',
                       aspect='auto',
                       vmin=limit[0], vmax=limit[1])
        ax.autoscale(False)

        if title is None:
            ax.set_title('True density ({0}), Rep:{1}, Grid:{2}s{3}m'.format(unit_str, rep, grid_res[0], grid_res[1]),
                         fontsize=24)
        else:
            ax.set_title(title, fontsize=24)

        plt.xlabel('Time (s)')
        plt.ylabel('Space, traffic direction $\mathbf{\Rightarrow}$')
        cax = fig.add_axes([0.95, 0.25, 0.01, 0.5])
        fig.colorbar(im, cax=cax, orientation='vertical')

        if save_fig is True:
            plt.savefig('truestate_{0}s{1}m_density.pdf'.format(grid_res[0], grid_res[1]))
            # clear the memory once generated
            plt.clf()
            plt.close()
        else:
            plt.draw()

    # deprecated
    def plot_true_speed_density_for_rep(self, grid_res, rep, unit='metric', speed_limit=(0, 40), density_limit=(0, 10)):
        """
        This function plots the true speed and density profile for replication
        :param grid_res: (s, m), the resolution of the grid
        :param rep: int, the replication id
        :param unit: 'metric'(m/s) and 'imperial' (mph)
        :param speed_limit: the color bar limit for speed field
        :param density_limit: th ecolor bar limit for the density field
        :return: Two figures
        """

        # ===============================================================
        # read the true speed data
        true_speed_file = self.__truestate_result_dir[rep] + 'truestate_{0}s{1}m_speed.txt'.format(grid_res[0],
                                                                                                   grid_res[1])
        speed_data = np.genfromtxt(true_speed_file, delimiter=',')
        speed_data = np.matrix(speed_data).T

        # read teh true density data
        true_density_file = self.__truestate_result_dir[rep] + 'truestate_{0}s{1}m_density.txt'.format(grid_res[0],
                                                                                                       grid_res[1])
        density_data = np.genfromtxt(true_density_file, delimiter=',')
        density_data = np.matrix(density_data).T

        # ===============================================================
        # change the unit
        if unit == 'metric':
            # all internal values are in metric, so plot directly
            speed = np.flipud(speed_data)
            density = np.flipud(density_data)
            unit_str = 'm/s'
        elif unit == 'imperial':
            speed = np.flipud(self.__metric2imperial(speed_data, 'speed'))
            density = np.flipud(self.__metric2imperial(density_data, 'density'))
            unit_str = 'mph'
        else:
            raise Exception('Error: Unrecognized unit for plotting speed.')

        # ===============================================================
        # remove a portion of the data
        speed[np.isnan(speed)] = -1
        density[np.isnan(density)] = -1
        invalid_data_idx = (speed <= 17.88) & (density <= 0.06215)
        print('invalid_data_idx shape: {0}'.format(invalid_data_idx.shape))
        speed[invalid_data_idx] = np.nan
        density[invalid_data_idx] = np.nan

        # ===============================================================
        # plot the speed
        fig = plt.figure(figsize=(15, 8), dpi=100)
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        im = ax.imshow(speed, cmap=plt.get_cmap('jet_r'),
                       interpolation='nearest',
                       aspect='auto',
                       vmin=speed_limit[0], vmax=speed_limit[1])
        ax.autoscale(False)

        ax.set_title('True speed, Rep:{0}, Grid:{1}s{2}m'.format(rep, grid_res[0], grid_res[1]))
        plt.xlabel('Time (min)')
        plt.ylabel('Space, traffic direction $\mathbf{\Rightarrow}$')
        cax = fig.add_axes([0.95, 0.25, 0.01, 0.5])
        fig.colorbar(im, cax=cax, orientation='vertical')
        plt.draw()

        # ===============================================================
        # plot the density
        fig = plt.figure(figsize=(15, 8), dpi=100)
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        im = ax.imshow(density, cmap=plt.get_cmap('jet'),
                       interpolation='nearest',
                       aspect='auto',
                       vmin=density_limit[0], vmax=density_limit[1])
        ax.autoscale(False)

        ax.set_title('True density, Rep:{0}, Grid:{1}s{2}m'.format(rep, grid_res[0], grid_res[1]))
        plt.xlabel('Time (min)')
        plt.ylabel('Space, traffic direction $\mathbf{\Rightarrow}$')
        cax = fig.add_axes([0.95, 0.25, 0.01, 0.5])
        fig.colorbar(im, cax=cax, orientation='vertical')

        plt.draw()

    def plot_true_queue_for_rep(self, grid_res, rep, unit='metric', limit=(0, 6600), fontsize=(40, 36, 34),
                                save_fig=False, title=None, fig_size=(10, 10)):
        """
        This function plots the true queue for the replication.
        :param grid_res: the grid for the true state
        :param rep: int, replication id
        :param unit: 'metric' (m) or 'imperial' (mile)
        :param limit: the ylim for the plot
        :param fontsize: a tuple of three int, font size for title, labels and legend, tick labels
        :param save_fig: True of False
        :param title: str, the tile of the figure
        :param fig_size: tuple, figure size
        :return: a plotted or saved figure
        """

        true_queue_file = self.__truestate_result_dir[rep] + \
                          'truestate_{0}s{1}m_queue.txt'.format(grid_res[0], grid_res[1])

        queue = np.genfromtxt(true_queue_file, delimiter=',')

        if unit == 'metric':
            unit_str = 'm'
        elif unit == 'imperial':
            queue = queue / 1609.34
            unit_str = 'mile'

        fig = plt.figure(figsize=fig_size, dpi=100)
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

        start_time = 60
        ax.plot(queue[start_time:], linewidth=2)

        ax.set_ylim(limit)
        ax.set_xlabel('Time (min)', fontsize=fontsize[1])
        ax.set_ylabel('Queue length (miles)', fontsize=fontsize[1])
        # ax.grid(True)

        # set x and y ticks
        x_ticks = np.array([start_time, 360, 720, 1080, 1440, 1800])
        x_ticklabels = 5 * x_ticks / 60  # in min
        x_ticks = x_ticks - start_time
        ax.set_xlim([x_ticks[0], x_ticks[-1]])
        ax.set_xticks(x_ticks)
        ax.set_xticklabels(x_ticklabels, fontsize=fontsize[2])
        ax.tick_params(axis='both', which='major', labelsize=fontsize[2])

        if title is None:
            ax.set_title(
                'True queueu length on {0}s{1}m, rep {2}, ({3})'.format(grid_res[0], grid_res[1], rep, unit_str),
                fontsize=fontsize[0])
        else:
            ax.set_title(title, fontsize=fontsize[0])

        if save_fig is True:
            plt.savefig('truestate_{0}s{1}m_queue.pdf'.format(grid_res[0], grid_res[1]),
                        bbox_inches='tight')
            # clear the memory once generated
            plt.clf()
            plt.close()
        else:
            plt.draw()

    def plot_true_traveltime_for_rep(self, grid_res, rep, unit='metric', limit=(0, 600),
                                     fig_size=(15, 8), fontsize=(40, 36, 34),
                                     save_fig=False, title=None):
        """
        This function plots the queue length
        :param grid_res: the grid for the true state
        :param rep: int, replication id
        :param unit: 'metric' (s) or 'imperial' (s)
        :param limit: the ylim for the plot
        :param fontsize: a tuple of three int, font size for title, labels and legend, tick labels
        :param save_fig: True of False
        :param title: str, the tile of the figure
        :param fig_size: tuple, figure size
        :return:
        """
        true_traveltime_file = self.__truestate_result_dir[rep] + 'truestate_{0}s{1}m_true_traveltime.txt' \
            .format(grid_res[0], grid_res[1])

        traveltime = np.genfromtxt(true_traveltime_file, delimiter=',')

        fig = plt.figure(figsize=fig_size, dpi=100)
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

        start_time = 60
        ax.plot(traveltime[start_time:] / 60.0, linewidth=2)

        ax.set_xlabel('Time (min)', fontsize=fontsize[1])
        ax.set_ylabel('Travel time (min)', fontsize=fontsize[1])

        # set x and y ticks
        x_ticks = np.array([start_time, 360, 720, 1080, 1440, 1800])
        x_ticklabels = 5 * x_ticks / 60  # in min
        x_ticks -= start_time

        ax.set_xticks(x_ticks)
        ax.set_xticklabels(x_ticklabels, fontsize=fontsize[2])
        ax.set_xlim([x_ticks[0], x_ticks[-1]])
        ax.tick_params(axis='both', which='major', labelsize=fontsize[2])

        if title is None:
            ax.set_title('True travel time on {0}s{1}m, rep {2}'.format(grid_res[0],
                                                                        grid_res[1], rep),
                         fontsize=fontsize[0])
        else:
            ax.set_title(title, fontsize=fontsize[0])
        if save_fig is True:
            plt.savefig('truestate_{0}s{1}m_tt.pdf'.format(grid_res[0], grid_res[1]), bbox_inches='tight')
            # clear the memory once generated
            plt.clf()
            plt.close()
        else:
            plt.draw()

    def plot_measured_traveltime_for_rep(self, grid_res, rep, unit='metric', limit=(0, 600),
                                         save_fig=False, title=None):
        """
        This function plots the measured travel time using bluetooth sensors
        :param grid_res: the grid for the true state
        :param rep: int, replication id
        :param unit: 'metric' (s) or 'imperial' (s)
        :param limit: the ylim for the plot
        :param save_fig: True of False
        :param title: str, the tile of the figure
        :return:
        """
        measured_traveltime_file = self.__truestate_result_dir[rep] + 'truestate_{0}s{1}m_measured_traveltime.txt' \
            .format(grid_res[0], grid_res[1])

        traveltime = np.genfromtxt(measured_traveltime_file, delimiter=',')

        fig = plt.figure(figsize=(15, 8), dpi=100)
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        ax.plot(traveltime)
        ax.set_ylabel('travel time (s)')

        if title is None:
            ax.set_title('Measured travel time on {0}s{1}m, rep {2}'.format(grid_res[0], grid_res[1], rep),
                         fontsize=24)
        else:
            ax.set_title(title, fontsize=24)
        if save_fig is True:
            plt.savefig('truestate_{0}s{1}m_tt.pdf'.format(grid_res[0], grid_res[1]))
            # clear the memory once generated
            plt.clf()
            plt.close()
        else:
            plt.draw()

    def plot_trueinst_and_measured_traveltime_for_rep(self, grid_res, replications, unit='metric', limit=(0, 600),
                                                      save_fig=False, title=None):
        """
        This function plots the instantaneous travel time computed from the true speed field and
            the measured travel time from Bluetooth sensors.
        :param grid_res: the grid for the true state
        :param replications: a list of replication ids (int)
        :param unit: 'metric' (m) or 'imperial' (mile)
        :param limit: the ylim for the plot
        :param save_fig: True of False
        :param title: str, the tile of the figure
        :return:
        """

        tt_meas_err = []
        tt_inst_err = []
        for rep in replications:

            measured_traveltime_file = self.__truestate_result_dir[rep] + \
                                       'truestate_{0}s{1}m_measured_traveltime.txt' \
                                           .format(grid_res[0], grid_res[1])

            truespeedfield_traveltime_file = self.__truestate_result_dir[rep] + \
                                             'truestate_{0}s{1}m_trueinst_traveltime.txt' \
                                                 .format(grid_res[0], grid_res[1])

            true_traveltime_file = self.__truestate_result_dir[rep] + 'truestate_{0}s{1}m_true_traveltime.txt' \
                .format(grid_res[0], grid_res[1])

            meas_traveltime = np.genfromtxt(measured_traveltime_file, delimiter=',')
            true_inst_traveltime = np.genfromtxt(truespeedfield_traveltime_file, delimiter=',')
            true_traveltime = np.genfromtxt(true_traveltime_file, delimiter=',')

            # ========================================================
            # compute the error
            t_start = 60  # throw away the first 5 min data due to the lack of true speed
            err_meas = np.abs(meas_traveltime[t_start:] - true_traveltime[t_start:])
            notNan = ~np.isnan(err_meas)
            L_tt_meas = err_meas[notNan].sum() / notNan.sum()
            tt_meas_err.append(L_tt_meas)

            err_inst = np.abs(true_inst_traveltime[t_start:] - true_traveltime[t_start:])
            notNan = ~np.isnan(err_inst)
            L_tt_inst = err_inst[notNan].sum() / notNan.sum()
            tt_inst_err.append(L_tt_inst)

            # ========================================================
            # visualization
            fig = plt.figure(figsize=(15, 8), dpi=100)
            ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
            ax.plot(meas_traveltime / 60.0, linewidth=2, color='r', label='measured')
            ax.plot(true_inst_traveltime / 60.0, linewidth=2, color='b', label='true instantaneous')
            ax.plot(true_traveltime / 60.0, linewidth=2, color='g', label='true')
            ax.set_ylabel('travel time for rep {0} (min)'.format(rep), fontsize=30)

            plt.legend()

            if title is None:
                ax.set_title('MAE: meas {0} s; inst {1} s'.format(L_tt_meas, L_tt_inst),
                             fontsize=30)
            else:
                ax.set_title(title, fontsize=30)

            if save_fig is True:
                plt.savefig('comparison_traveltime_rep{2}_{0}s{1}m_tt.pdf'.format(grid_res[0], grid_res[1],
                                                                                  rep), bbox_inches='tight')
                # clear the memory once generated
                plt.clf()
                plt.close()
            else:
                plt.draw()

        print('mean MAE:\n --- meas tt: {0} s\n --- inst tt: {1} s'.format(np.mean(tt_meas_err),
                                                                           np.mean(tt_inst_err)))

    # ================================================================================ #
    # plot the sensor measurements in stripes
    # ================================================================================ #
    def plot_speed_meas_for_config(self, rep, config_id, unit='metric', limit=(0, 40)):
        """
        This function plots the speed measurement as a thin stripe to illustrate how the algorithm works.
        It is built upon the plot_speed_for_scenario function. It basically removes the data not needed and convert the
        estimated speed in the cell of the speed sensor to a thin strip.
        - The sensor data files are timestamp(s), v (kph), count/agg_interval
        :param rep: int, the replciation id
        :param config_id: string, the configuration id
        :param unit: string, 'metric' or 'imperial'
        :param limit: the limit in corresponding unit
        :return: a figure
        """

        if config_id not in self.all_config.keys():
            raise Exception('Config {0} is not in loaded configurations'.format(config_id))

        # read the data from file [timestamp (s), v (kph), count/agg_interval]
        speed = OrderedDict()
        num_steps = None
        for sensor_id in self.all_config[config_id]['sensors'].keys():

            sensor_file = self.__sensor_data_dir[rep] + sensor_id + '.txt'

            if not exists(sensor_file):
                raise Exception('Error: sensor data for {0} does not exist.'.format(sensor_file))

            # timestamp(s), speed(kph), count(veh)
            data = np.genfromtxt(sensor_file, delimiter=',')

            # # If data is missing, speed value will be -1, replace with nan
            # nan_idx = (data < 0)
            # data[nan_idx] = np.nan
            data = self.__clean_data(data)

            # change the unit. 'metric': m/s; 'imperial': mph
            if unit == 'metric':
                speed[sensor_id] = data[:, 1] / 3.6
                unit_str = 'm/s'
            elif unit == 'imperial':
                speed[sensor_id] = data[:, 1] / 1.609
                unit_str = 'mph'
            else:
                raise Exception('Unrecognized unit {0}'.format(unit))

            # # if not car detected, speed value was generated as inf, change speed to maximum speed.
            # inf_idx = (data == np.inf)
            # data[inf_idx] = limit[1]
            print('num of inf in data:{0}'.format(sum(sum(np.isinf(data)))))

            # get the number of steps, if not the same, then error
            if num_steps is None:
                num_steps = len(speed[sensor_id])
                agg_interval = data[1, 0] - data[0, 0]  # in seconds
            elif num_steps != len(speed[sensor_id]):
                raise Exception('The sensor data length does not match.')

        # plot the measurement in a wide line
        fig = plt.figure(figsize=(20, 10))
        ax = fig.add_subplot(111)
        width_strip = 5
        loc_entrance = self.space_grid[0]
        loc_exit = self.space_grid[-1]

        sensor_ytick = {}
        sensor_ytick['tick'] = []
        sensor_ytick['label'] = []

        x_grid = np.linspace(0, num_steps * agg_interval, num_steps + 1)

        # set up color Map
        cNorm = colors.Normalize(vmin=limit[0], vmax=limit[1])
        scalarMap = cm.ScalarMappable(norm=cNorm, cmap=plt.get_cmap('jet_r'))
        scalarMap.set_array([])

        for sensor_id in self.all_config[config_id]['sensors'].keys():

            # absolute location regarding to the entrance
            abs_loc = self.__get_abs_loc(self.all_config[config_id]['sensors'][sensor_id]['section'],
                                         self.all_config[config_id]['sensors'][sensor_id]['distance'])

            sensor_ytick['tick'].append(abs_loc)
            sensor_ytick['label'].append(sensor_id)

            # plot small line segment with correponding color
            for cell in range(0, len(speed[sensor_id])):
                if ~np.isnan(speed[sensor_id][cell]):
                    color_val = scalarMap.to_rgba(speed[sensor_id][cell])
                else:
                    color_val = (1, 1, 1, 1)  # if nan, then set color as opaque white
                ax.plot([x_grid[cell], x_grid[cell + 1]],
                        [abs_loc, abs_loc],
                        color=color_val, linewidth=width_strip)

        ax.set_yticks(sensor_ytick['tick'])
        ax.set_yticklabels(sensor_ytick['label'])
        ax.set_xlim([0, num_steps * agg_interval])
        ax.set_ylim([-10, sensor_ytick['tick'][-1] + 10])
        plt.title('Speed measurement in configuration {0} ({1})'.format(config_id, unit_str), fontsize=24)
        plt.xlabel('Time (s)', fontsize=20)
        plt.ylabel('Traffic direction $\mathbf{\Rightarrow}$', fontsize=20)
        ax.grid()
        cax = fig.add_axes([0.95, 0.25, 0.01, 0.5])
        plt.colorbar(scalarMap, cax=cax, orientation='vertical')
        plt.draw()

    def plot_speed_meas_for_sensors(self, rep, sensor_ids, unit='metric', limit=(0, 40)):
        """
        This function plots the strip measurement for sensors. Similar to plot_speed_meas_for_config
            but without requiring a configuration.
        :param rep: the replication
        :param sensor_ids: a list of sensor ids to add to the plot
        :param unit: 'metric' or 'imperial'
        :param limit: limit of the plot to the corresponding unit
        :return:
        """

        # load the sensor configuration file to get the parameters for the sensors.
        sensor_paras = OrderedDict()
        with open(self.__sensor_data_log[rep]) as f:
            for line in f:
                items = line.strip().split(';')

                # check sensor id in this line
                line_sensor_id = items[0].split(':')[1]

                if line_sensor_id not in sensor_ids:
                    continue
                else:
                    sensor_id = line_sensor_id
                    sensor_paras[sensor_id] = OrderedDict()
                    # we only need section and distance
                    for item in items:
                        tup = item.split(':')
                        if tup[0] == 'section':
                            sensor_paras[sensor_id]['section'] = float(tup[1])
                        elif tup[0] == 'distance':
                            sensor_paras[sensor_id]['distance'] = float(tup[1])

        # read the data from file [timestamp (s), v (kph), count/agg_interval]
        speed = OrderedDict()
        num_steps = None
        for sensor_id in sensor_ids:

            sensor_file = self.__sensor_data_dir[rep] + sensor_id + '.txt'

            if not exists(sensor_file):
                raise Exception('Error: sensor data for {0} does not exist.'.format(sensor_file))

            # timestamp(s), speed(kph), count(veh)
            data = np.genfromtxt(sensor_file, delimiter=',')
            # replace inf, -1, with nan
            nan_idx = (data == np.inf) | (data < 0)
            data[nan_idx] = np.nan

            # change the unit. 'metric': m/s; 'imperial': mph
            if unit == 'metric':
                speed[sensor_id] = data[:, 1] / 3.6
                unit_str = 'm/s'
            elif unit == 'imperial':
                speed[sensor_id] = data[:, 1] / 1.609
                unit_str = 'mph'
            else:
                raise Exception('Unrecognized unit {0}'.format(unit))

            # get the number of steps, if not the same, then error
            if num_steps is None:
                num_steps = len(speed[sensor_id])
                agg_interval = data[1, 0] - data[0, 0]  # in seconds
            elif num_steps != len(speed[sensor_id]):
                raise Exception('The sensor data length does not match.')

        # plot the measurement as wide color coded lines
        fig = plt.figure(figsize=(20, 10))
        ax = fig.add_subplot(111)
        width_strip = 5
        loc_entrance = self.space_grid[0]
        loc_exit = self.space_grid[-1]

        sensor_ytick = {}
        sensor_ytick['tick'] = []
        sensor_ytick['label'] = []

        x_grid = np.linspace(0, num_steps * agg_interval, num_steps + 1)

        # set up color Map
        cNorm = colors.Normalize(vmin=limit[0], vmax=limit[1])
        scalarMap = cm.ScalarMappable(norm=cNorm, cmap=plt.get_cmap('jet_r'))
        scalarMap.set_array([])

        for sensor_id in sensor_ids:

            # absolute location regarding to the entrance
            abs_loc = self.__get_abs_loc(sensor_paras[sensor_id]['section'],
                                         sensor_paras[sensor_id]['distance'])

            sensor_ytick['tick'].append(abs_loc)
            sensor_ytick['label'].append(sensor_id)

            # plot small line segment with correponding color
            for cell in range(0, len(speed[sensor_id])):
                color_val = scalarMap.to_rgba(speed[sensor_id][cell])
                ax.plot([x_grid[cell], x_grid[cell + 1]],
                        [abs_loc, abs_loc],
                        color=color_val, linewidth=width_strip)

        ax.set_yticks(sensor_ytick['tick'])
        ax.set_yticklabels(sensor_ytick['label'])
        ax.set_xlim([0, num_steps * agg_interval])
        plt.title('Speed measurement from sensors ({0})'.format(unit_str), fontsize=24)
        plt.xlabel('Time (s)', fontsize=20)
        plt.ylabel('Traffic direction $\mathbf{\Rightarrow}$', fontsize=20)
        ax.grid()
        cax = fig.add_axes([0.95, 0.25, 0.01, 0.5])
        plt.colorbar(scalarMap, cax=cax, orientation='vertical')
        plt.draw()
        # =============================================
        # End approach for plotting using plt with line color
        # =============================================

    def plot_flow_meas_time_series_sensors(self, rep, sensor_id_list):
        """
        This function plots the flow measurement in time series in line plots
        :param rep: int, the replication id
        :param sensor_id_list: a list of sensor ids
        :return: a figure
        """

        flow = OrderedDict()

        fig = plt.figure(figsize=(18, 8.5), dpi=100)
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

        for sensor_id in sensor_id_list:

            sensor_file = self.__sensor_data_dir[rep] + sensor_id + '.txt'

            if not exists(sensor_file):
                raise Exception('Error: {0} does not exist.'.format(sensor_file))

            # timestamp(s), speed(kph), count(veh)
            data = np.genfromtxt(sensor_file, delimiter=',')
            # replace inf, -1, with nan
            nan_idx = (data == np.inf) | (data < 0)
            data[nan_idx] = np.nan

            agg_time = data[1, 0] - data[0, 0]
            flow[sensor_id] = data[:, 2] * 3600.0 / agg_time

            plt.plot(data[:, 0], flow[sensor_id], linewidth=2.0,
                     label='{0}'.format(sensor_id))

        plt.title('Flow measurement', fontsize=24)
        plt.xlabel('Time', fontsize=24)
        plt.ylabel('Flow (veh/hr)', fontsize=24)
        plt.legend(loc='best')
        plt.grid(True)

        plt.draw()

    def plot_speed_meas_time_series_sensors(self, rep, sensor_id_list, unit='metric', limit=(0, 40)):
        """
        This function plots the speed measurement in time series in line plots
        :param rep: int, the replication id
        :param sensor_id_list: a list of sensor ids
        :param unit: 'metric' (m/s) or 'imperial' (mph)
        :param limit: the y limit
        :return: a figure
        """

        speed = OrderedDict()

        fig = plt.figure(figsize=(18, 8.5), dpi=100)
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

        for sensor_id in sensor_id_list:

            sensor_file = self.__sensor_data_dir[rep] + sensor_id + '.txt'

            if not exists(sensor_file):
                raise Exception('Error: {0} does not exist.'.format(sensor_file))

            # timestamp(s), speed(kph), count(veh)
            data = np.genfromtxt(sensor_file, delimiter=',')
            # replace inf, -1, with nan
            nan_idx = (data == np.inf) | (data < 0)
            data[nan_idx] = np.nan

            if unit == 'imperial':
                agg_time = data[1, 0] - data[0, 0]
                speed[sensor_id] = data[:, 1] / 1.609
            elif unit == 'metric':
                # convert to m/s
                speed[sensor_id] = data[:, 1] * 1000.0 / 3600.0

            plt.plot(data[:, 0], speed[sensor_id], linewidth=2.0,
                     label='{0}'.format(sensor_id))

        plt.ylim(limit)
        plt.title('Speed measurement', fontsize=24)
        plt.xlabel('Time', fontsize=24)
        plt.ylabel('Speed (mph)', fontsize=24)
        plt.legend(loc='best')
        plt.grid(True)

        plt.draw()

    # ================================================================================ #
    # Utility functions
    # ================================================================================ #

    def __true_state_name(self, grid):
        """
        This function returns the file name that saves the trues states on (dt,dx) grid
        Should agree on this function.
        :param grid: (s, m) tuple, time space resolution
        :return: the string of the true density file
        """
        return 'truestate_{1}s{2}m_density.txt'.format(self.workzone, grid[0], grid[1])

    def __get_abs_loc(self, sec_id, dist):
        """
        This function returns the absolute post meter location of the location specified on the main freeway by
            section id and dist on the section
        :param sec_id: the section id on the freeway
        :param dist: the distant to the entrance of this section
        :return: float loc,
        """
        loc = 0.0

        for sec in self.workzone_topo['fwy_sec_order']:
            if sec_id == sec:
                loc = loc + dist
                return loc
            else:
                loc = loc + self.workzone_topo['sections'][sec]['length']

        raise Exception('Error: the network does not contain the section {0}'.format(sec_id))

    def __get_rel_loc(self, loc):
        """
        This function returns the relative location on the highway as (sect, dist)
        :param loc: the absolute location in meters
        :return: float tuple, (sect, dist)
        """
        dist = loc
        sect = 0

        for sec in self.workzone_topo['fwy_sec_order']:
            if self.workzone_topo['sections'][sec]['length'] >= dist:
                # found the section
                sect = sec
                return sect, dist
            else:
                dist -= self.workzone_topo['sections'][sec]['length']

        raise Exception('Error: location {0} is out of the network'.format(loc))

    def __isNumber(self, string):
        """
        This function returns true if the string can be converted to a float or int
        :param string: a string
        :return: True if the string can be converted to a float or int
        """
        try:
            float(string)
            return True
        except ValueError:
            return False

    def __assertType(self, string):
        """
        This function asserts the correct type of the string value:
            - int if the string is a int number
            - list of ints if the string is a comma separated list
            - float if the string can be converted to a float;
            - list of floats if the string is a comma separated list
            - string otherwise
        :param string: the input string
        :return: a float or the original string
        """
        converted_list = []
        items = string.strip().split(',')
        flag_float = False

        for v in items:

            # if int, convert to int
            try:
                int(v)
                converted_list.append(int(v))
                continue
            except ValueError:
                pass

            # if float, convert to int
            try:
                float(v)
                converted_list.append(float(v))
                flag_float = True
                continue
            except ValueError:
                pass

            # otherwise, just append the original string
            converted_list.append(v)

        if flag_float is True:
            # if one element is float, then convert all to float
            try:
                converted_list = [float(i) for i in converted_list]
            except ValueError:
                raise Exception('Error: List can not be converted to float: {0}.'.format(converted_list))

        if len(converted_list) == 1:
            return converted_list[0]
        else:
            return converted_list

    def __sensorDictToString(self, sensor_id, sensor_att):
        """
        This function convers the sensor key-value store to a string in the following format. Entries are separated by ;
        The entry name and value are separated by :
        :param sensor_id: the sensor id
        :param sensor_att: key-value store
        :return: a string in the above defined format
        """

        line = []
        # append the first id entry
        entry = ':'.join(['id', sensor_id])
        line.append(entry)

        # add the distance offset to the sensors
        if 'offset_dist' in sensor_att.keys():
            sensor_att['distance'] += sensor_att['offset_dist']

        # append other attribute entries
        for key in sensor_att.keys():

            if key != 'id':
                # Do Not repeat id entry
                # NOTE: the range list will be converted to a string with []
                entry = ':'.join([key, str(sensor_att[key])])
                line.append(entry)

        # join all entries and return
        return ';'.join(line)

    @staticmethod
    def __metric2imperial(value=np.zeros((1, 1)), option='speed'):
        """
        A utility function which converts the metric (m, s, m/s) to imperial (mile, hour, m/h)
        :param value: float, np.array, or np.matrix. the to be converted value
        :param option: 'speed', 'density'
        :return: converted value
        """

        if type(value) is float or type(value) is np.float64 \
                or type(value) is np.ndarray or type(value) is np.matrix:
            if option == 'speed':
                return value * 3600.0 / 1609.34
            elif option == 'density':
                return value * 1609.34
            elif option == 'distance':
                return value / 1609.34
            else:
                raise Exception('Error: Unrecognized unit conversion option.')
        else:
            print(type(value))
            raise Exception('Error: Unrecognized value type for unit conversion.')

    # ================================================================================ #
    # Calibrating the fundamental diagram
    # ================================================================================ #
    def calib_FD_on_grid(self, grid=None, replications=None):
        """
        This function calibrate the fundamental diagram using the true state data on grid
        :param grid: (s, m), the grid resolution of the true states data
        :param replications: the list of replications
        :return: a plot with all parameters
        """

        for rep in replications:
            # ===========================================================================================
            # read the density data
            true_density_file = self.__truestate_result_dir[rep] + 'truestate_{0}s{1}m_density.txt'.format(grid[0],
                                                                                                           grid[1])

            density_data = np.genfromtxt(true_density_file, delimiter=',')
            density_data = np.matrix(density_data).T

            true_speed_file = self.__truestate_result_dir[rep] + 'truestate_{0}s{1}m_speed.txt'.format(grid[0], grid[1])
            speed_data = np.genfromtxt(true_speed_file, delimiter=',')
            speed_data = np.matrix(speed_data).T

            # Use the data only at a few fixed locations
            num_sensors = 4
            locations = [int(round(i)) for i in np.linspace(10, density_data.shape[0] - 10, num_sensors)]
            print('Extracting data in rows:{0} out from {1} rows'.format(locations, density_data.shape[0]))

            # convert density from veh/m to veh/mile and speed from m/s to mph
            density_array = np.squeeze(np.array(density_data[locations, :].reshape(-1, ))) * 1609.34

            speed_array = np.squeeze(np.array(speed_data[locations, :].reshape(-1, ))) * 3600.0 / 1609.34
            flow_array = density_array * speed_array

            # =============================================================
            # remove the nan values
            valid_idx = ~np.isnan(speed_array)
            speed_array = speed_array[valid_idx]
            density_array = density_array[valid_idx]
            flow_array = flow_array[valid_idx]

            # =============================================================
            # remove the noisy point in a triangle where v < 20 mph & w >= w_noise from point (0,0) and (rho_noise,0))
            v_thres = 40
            rho_noise = 400
            w_noise = -12
            rho_thres = 135

            noise_idx = (flow_array <= v_thres * (density_array - 50)) & (
                flow_array <= w_noise * (density_array - rho_noise))
            noisy_density = density_array[noise_idx]
            noisy_flow = flow_array[noise_idx]

            # reassign weights to the density and speed
            # density_array, speed_array, flow_array = self.__weight_samples(density_array[~noise_idx],
            #                                                                speed_array[~noise_idx],
            #                                                                flow_array[~noise_idx])

            ff_index = (density_array <= rho_thres) & (~noise_idx)
            cg_index = (density_array > rho_thres) & (~noise_idx)
            rest_index = ~(ff_index | cg_index)

            # threshold for freeflow data mph
            # v_thres = 30
            # rho_thres = 100 # veh/mile
            # # ff_index = speed_array >= v_thres
            # # cg_index = (speed_array <= v_thres)
            # ff_index = (density_array <= rho_thres) | (speed_array >= v_thres)
            # cg_index = ~ff_index
            # rest_index = ~( ff_index | cg_index )

            ff_speed = speed_array[ff_index]
            ff_density = density_array[ff_index]

            ff_flow = flow_array[ff_index]
            cg_speed = speed_array[cg_index]
            cg_density = density_array[cg_index]
            cg_flow = flow_array[cg_index]
            rest_speed = speed_array[rest_index]
            rest_density = density_array[rest_index]
            rest_flow = flow_array[rest_index]

            # ===========================================================================================
            # fit a quadratic linear line to the data in the freeflow regime
            funcQuadFit = lambda vm_beta, rho: vm_beta[0] * rho - np.power(rho, 2) * vm_beta[0] / vm_beta[1]
            # funErr = lambda vm_beta, rho, q : funcQuadFit(vm_beta, rho) - q
            # vm_beta_init = [80, 600]     # initial guess of vm is 60

            # updated vresion, whic assumes beta is very large to get approximated TFD
            beta = 1000
            funcQuadFitVm = lambda vm, rho: vm * rho - np.power(rho, 2) * vm / beta
            funErr = lambda vm, rho, q: funcQuadFitVm(vm, rho) - q
            vm_init = [80]  # initial guess of vm is 60

            vm_est, success = optimize.leastsq(funErr, vm_init, args=(ff_density, ff_flow))

            all_vm = vm_est[0]
            # all_beta = vm_est[1]
            all_beta = beta
            print('vm:{0}; beta:{1}'.format(all_vm, all_beta))

            # fit a line to the congested regime
            rho_m = 500  # veh/m => 240.35 veh/mile
            # fit a linear line to the congested section through intercept (0.17, 0)
            cg_density = cg_density[:, np.newaxis]
            # print('shape of cg_flow: {0}; length:{1}'.format(cg_flow.shape, len(cg_flow)))
            wc, _, _, _ = np.linalg.lstsq(cg_density - rho_m, cg_flow)
            # print('end fitting cong')
            wc = wc[0]
            print('wc:{0}'.format(wc))
            print(
                'sqrt({0})'.format(
                    np.power(wc * all_beta - all_vm * all_beta, 2) - 4 * all_vm * (-wc * all_beta * rho_m)))

            # compute rho_c and rho_m
            rho_c = (-(wc * all_beta - all_vm * all_beta) -
                     np.sqrt(
                         np.power(wc * all_beta - all_vm * all_beta, 2) - 4 * all_vm * (-wc * all_beta * rho_m))) / (
                        2 * all_vm)
            # rho_c = 150
            q_max = funcQuadFit([all_vm, all_beta], rho_c)

            # ===========================================================================================
            # deprecated fitting
            # fit a linear line to the freeflow section
            # ff_density = ff_density[:, np.newaxis]
            # vf, _, _, _ = np.linalg.lstsq( ff_density, ff_flow )
            # print('vf:{0}, with std:{1}'.format(vf, np.std(ff_speed)))

            # ===========================================================================================
            # fit a model to congested section
            # ------------------------------------------------------------
            # second order fit
            # coe = np.polyfit( cg_density, cg_flow, 2 )
            # wc = coe[0]
            # ------------------------------------------------------------

            # ------------------------------------------------------------
            # second order polyfit with forced intercept at (rho_max,0)
            # center data
            # x_c = cg_density - rho_max
            # zero_index= x_c==0
            # y_d = cg_flow[~zero_index]/x_c[~zero_index]
            # coe = np.polyfit( x_c[~zero_index], y_d, 1)
            # # flow = a*(desnity-0.17)^2 + b*(desnity-0.17) + 0
            # a = coe[0]
            # b = coe[1]
            # ------------------------------------------------------------
            # second order polyfit with forced intercept at (0,0), (rho_max, 0)
            # x_s = -cg_density/rho_max + 1
            # zero_index = cg_density == 0
            # y_x = cg_flow[~zero_index]/cg_density[~zero_index]
            # coe = np.polyfit( x_s[~zero_index], y_x, 1 )
            # b = coe[0]
            # a = -b/rho_max

            # print 'wc:{0}'.format(wc)

            fig_window = plt.figure(figsize=(15, 8), dpi=100)
            fig = fig_window.add_subplot(111)

            # scatter freeflow
            plt.scatter(ff_density, ff_flow, color='g')
            dens = np.linspace(0, rho_c, 100)
            plt.plot(dens, funcQuadFit([all_vm, all_beta], dens), 'r-', linewidth=2.0)

            # scatter congestion
            plt.scatter(cg_density, cg_flow, color='k')
            dens = np.linspace(rho_c, rho_m, 100)
            plt.plot(dens, wc * (dens - rho_m), 'r-', linewidth=2.0)

            # plot rest points
            plt.scatter(noisy_density, noisy_flow, color='b')

            plt.title('Fundamental diagram grid {0}sx{1}m for rep {2}'.format(grid[0], grid[1], rep), fontsize=24)
            plt.xlabel('Traffic density (veh/mile)', fontsize=24)
            plt.ylabel('Traffic flow (veh/hr)', fontsize=24)

            text_str = r'freeflow: $q = v_m\rho - v_m\rho^2/\beta$' + '\n' \
                                                                      r'congflow: $q = w(\rho - \rho_m)$' + '\n' + \
                       r' $v_m$=   {0} mph ({1} m/s)'.format(np.round(all_vm, 2),
                                                             np.round(all_vm * 1609.34 / 3600.0, 2)) + '\n' + \
                       r' $\beta$=    {0} veh/mile ({1} veh/m)'.format(np.round(all_beta, 2),
                                                                       np.round(all_beta / 1609.34, 4)) + '\n' + \
                       r' $w$=    {0} mph ({1} m/s)'.format(np.round(wc, 2),
                                                            np.round(wc * 1609.34 / 3600.0, 2)) + '\n' + \
                       r' $\rho_c$=   {0} veh/mile ({1} veh/m)'.format(np.round(rho_c, 2),
                                                                       np.round(rho_c / 1609.34, 4)) + '\n' + \
                       r' $\rho_m$=   {0} veh/mile ({1} veh/m)'.format(np.round(rho_m, 2),
                                                                       np.round(rho_m / 1609.34, 4)) + '\n' + \
                       r' $q_m$=   {0} veh/hr ({1} veh/s)'.format(np.round(q_max, 2), np.round(q_max / 3600.0, 4))

            anchored_text = AnchoredText(text_str, loc=1)
            fig.add_artist(anchored_text)

            plt.grid(True)

            plt.draw()

            # Flow speed
            fig = plt.figure(figsize=(10, 7), dpi=100)
            ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
            plt.scatter(ff_speed, ff_flow, color='g')
            plt.grid(True)
            plt.scatter(cg_speed, cg_flow, color='k')
            plt.xlabel('speed (mph)')
            plt.ylabel('Flow (veh/hr)')
            plt.title('Flow-speed on Grid {0} rep {1}'.format(grid, rep))

            # now plot the q-k
            # fig = plt.figure( figsize=(18,10), dpi=100 )
            # ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
            # plt.scatter( density_list, flow_list )
            # plt.grid(True)
            # plt.xlabel('Density (veh/m)')
            # plt.ylabel('Flow (veh/s)')
            # plt.title('Flow-Density on Grid {0}'.format(grid))
            #
            # fig = plt.figure( figsize=(18,10), dpi=100 )
            # ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
            # plt.scatter( density_list, speed_list )
            # plt.grid(True)
            # plt.xlabel('Density (veh/m)')
            # plt.ylabel('Speed (m/s)')
            # plt.title('Speed-Density on Grid {0}'.format(grid))


            plt.draw()

    def calib_FD_using_sensor_data(self, data_files, agg_interval, plot_meas=False):
        """
        This function calibrate the fundamental diagram using virtual RTMS data generated by different replciations
        :param data_files: a list of data files
            data_file: timestamp(s), speed (kph), count
        :param agg_interval: int, seconds, aggregate the data to agg_interval to reduce quantization error
        :param plot_meas: True or False, plot the measurements,
        :return:
        """
        all_speed = []
        all_flow = []

        for data_f in data_files:

            if plot_meas is True:
                fig = plt.figure(figsize=(10, 7), dpi=100)
                ax_speed = fig.add_axes([0.1, 0.1, 0.8, 0.8])
                plt.hold(True)

                fig = plt.figure(figsize=(10, 7), dpi=100)
                ax_flow = fig.add_axes([0.1, 0.1, 0.8, 0.8])
                plt.hold(True)

            # ===========================================================================================
            # read the sensor data file
            if not exists(data_f):
                raise Exception('Error: sensor data for {0} does not exist.'.format(data_f))

            # timestamp(s), speed(kph), count(veh)
            data = np.genfromtxt(data_f, delimiter=',')
            detection_cycle = data[1, 0] - data[0, 0]

            # time stamp is shifted, hence the last row is not effective
            data = data[:-1, :]
            data[:, 0] = data[:, 0] + detection_cycle

            # ============================================
            # aggregate the data into aggregation intervals
            no_veh_idx = np.isinf(data[:, 1])
            data[no_veh_idx, 1] = 0.0
            data[no_veh_idx, 2] = 0.0
            missing_data_idx = np.isnan(data[:, 1]) | (data[:, 1] < 0)
            data[missing_data_idx, 1] = np.nan
            data[missing_data_idx, 2] = np.nan

            agg_steps = agg_interval / detection_cycle
            agg_data = []
            for i in range(0, int(data.shape[0] / agg_steps)):
                start_idx = int(i * agg_steps)
                end_idx = int((i + 1) * agg_steps)
                avg_speed = np.sum(data[start_idx:end_idx, 1] * data[start_idx:end_idx, 2]) / \
                            np.sum(data[start_idx:end_idx, 2])
                cum_count = np.sum(data[start_idx:end_idx, 2])
                agg_data.append([data[end_idx - 1, 0], avg_speed, cum_count])

            agg_data = np.array(agg_data)

            if plot_meas is True:
                # visualize the smoothing
                ax_speed.plot(data[:, 0], data[:, 1] / 1.609, color='b',
                              linewidth=2, label='30 s')
                plt.hold(True)
                ax_speed.plot(agg_data[:, 0], agg_data[:, 1] / 1.609, color='r',
                              linewidth=2, label='{0} s'.format(agg_interval))

                ax_flow.plot(data[:, 0], data[:, 2] * 3600.0 / detection_cycle, color='b',
                             linewidth=2, label='30 s')
                plt.hold(True)
                ax_flow.plot(agg_data[:, 0], agg_data[:, 2] * 3600.0 / agg_interval, color='r',
                             linewidth=2, label='{0} s'.format(agg_interval))

            print('aggregated {0} measurements to {1}'.format(data.shape[0], agg_data.shape[0]))
            data = agg_data

            # ============================================
            # clean data, remove inf, -1, and nan values, using the speed entry
            valid_idx = ~((np.isnan(data[:, 1]) | np.isinf(data[:, 1])) | (data[:, 1] < 0))
            cleaned_data = data[valid_idx, :]
            print('Got clean data points {0}/{1}'.format(cleaned_data.shape[0],
                                                         data.shape[0]))

            # ============================================
            # convert speed from kph to mph, and count to flow
            cleaned_data[:, 1] = cleaned_data[:, 1] / 1.609
            cleaned_data[:, 2] = cleaned_data[:, 2] * 3600.0 / agg_interval

            # ============================================
            # append data points
            all_speed = np.concatenate([all_speed, cleaned_data[:, 1]])
            all_flow = np.concatenate([all_flow, cleaned_data[:, 2]])

            if plot_meas is True:
                ax_flow.set_title('flow {0}'.format(data_f.split('/')[-1]))
                ax_speed.set_title('speed {0}'.format(data_f.split('/')[-1]))

        # =============================================================
        # data set for calibration
        speed_array = all_speed.reshape(1, all_speed.size).squeeze()
        flow_array = all_flow.reshape(1, all_flow.size).squeeze()
        density_array = flow_array / speed_array

        # =============================================================
        # remove the noisy outlier point in a triangle
        # where v < 20 mph & w >= w_noise from point (0,0) and (rho_noise,0))
        v_thres = 40
        flow_thres = 500
        rho_noise = 300
        w_noise = -12
        rho_thres = 150

        noise_idx = (speed_array <= 1) | ((density_array < rho_noise) & (speed_array < v_thres))
        noisy_density = density_array[noise_idx]
        noisy_speed = speed_array[noise_idx]
        noisy_flow = flow_array[noise_idx]

        ff_index = (speed_array > v_thres) & (~noise_idx)
        cg_index = (speed_array <= v_thres) & (~noise_idx)
        rest_index = ~(ff_index | cg_index)

        ff_speed = speed_array[ff_index]
        ff_density = density_array[ff_index]
        ff_flow = flow_array[ff_index]

        cg_speed = speed_array[cg_index]
        cg_density = density_array[cg_index]
        cg_flow = flow_array[cg_index]

        # ===========================================================================================
        # fit a quadratic linear line to the data in the free flow regime
        funcQuadFit = lambda vm_beta, rho: vm_beta[0] * rho - np.power(rho, 2) * vm_beta[0] / vm_beta[1]

        # updated vresion, whic assumes beta is very large to get approximated TFD
        beta = 1000
        funcQuadFitVm = lambda vm, rho: vm * rho - np.power(rho, 2) * vm / beta
        funErr = lambda vm, rho, q: funcQuadFitVm(vm, rho) - q
        vm_init = [80]  # initial guess of vm

        vm_est, success = optimize.leastsq(funErr, vm_init, args=(ff_density, ff_flow))

        all_vm = vm_est[0]
        # all_beta = vm_est[1]
        all_beta = beta
        print('vm:{0}; beta:{1}'.format(all_vm, all_beta))

        # ====================================================================
        # Fit a linear function to the congested regime
        preset_rhom = 500
        # funcLinearCong = lambda w_rhom, rho: w_rhom[0]*(rho-w_rhom[1])
        funcLinearCong = lambda w_rhom, rho: w_rhom[0] * (rho - preset_rhom)
        funErrCong = lambda w_rhom, rho, q: funcLinearCong(w_rhom, rho) - q

        w_rhom_init = [-10, 500.0]
        est_w_rhom, success = optimize.leastsq(funErrCong, w_rhom_init, args=(cg_density, cg_flow))
        wc = est_w_rhom[0]
        rho_m = preset_rhom
        print('wc:{0}'.format(wc))
        print(
            'sqrt({0})'.format(np.power(wc * all_beta - all_vm * all_beta, 2) - 4 * all_vm * (-wc * all_beta * rho_m)))

        # compute rho_c and rho_m
        rho_c = (-(wc * all_beta - all_vm * all_beta) -
                 np.sqrt(np.power(wc * all_beta - all_vm * all_beta, 2) - 4 * all_vm * (-wc * all_beta * rho_m))) / (
                    2 * all_vm)
        # rho_c = 150
        q_max = funcQuadFit([all_vm, all_beta], rho_c)

        # ====================================================================
        fig_window = plt.figure(figsize=(15, 8), dpi=100)
        fig = fig_window.add_subplot(111)

        # scatter freeflow
        plt.scatter(ff_density, ff_flow, color='g')
        dens = np.linspace(0, rho_c, 100)
        plt.plot(dens, funcQuadFit([all_vm, all_beta], dens), 'r-', linewidth=2.0)

        # scatter congestion
        plt.scatter(cg_density, cg_flow, color='k')
        dens = np.linspace(rho_c, rho_m, 100)
        plt.plot(dens, wc * (dens - rho_m), 'r-', linewidth=2.0)

        # plot rest points
        plt.scatter(noisy_density, noisy_flow, color='b')

        plt.title('Fundamental diagram, agg {0} s'.format(agg_interval), fontsize=24)
        plt.xlabel('Traffic density (veh/mile)', fontsize=24)
        plt.xlim([0, 800])
        plt.ylabel('Traffic flow (veh/hr)', fontsize=24)

        text_str = r'freeflow: $q = v_m\rho - v_m\rho^2/\beta$' + '\n' \
                                                                  r'congflow: $q = w(\rho - \rho_m)$' + '\n' + \
                   r' $v_m$=   {0} mph ({1} m/s)'.format(np.round(all_vm, 2),
                                                         np.round(all_vm * 1609.34 / 3600.0, 2)) + '\n' + \
                   r' $\beta$=    {0} veh/mile ({1} veh/m)'.format(np.round(all_beta, 2),
                                                                   np.round(all_beta / 1609.34, 4)) + '\n' + \
                   r' $w$=    {0} mph ({1} m/s)'.format(np.round(wc, 2), np.round(wc * 1609.34 / 3600.0, 2)) + '\n' + \
                   r' $\rho_c$=   {0} veh/mile ({1} veh/m)'.format(np.round(rho_c, 2),
                                                                   np.round(rho_c / 1609.34, 4)) + '\n' + \
                   r' $\rho_m$=   {0} veh/mile ({1} veh/m)'.format(np.round(rho_m, 2),
                                                                   np.round(rho_m / 1609.34, 4)) + '\n' + \
                   r' $q_m$=   {0} veh/hr ({1} veh/s)'.format(np.round(q_max, 2), np.round(q_max / 3600.0, 4))

        anchored_text = AnchoredText(text_str, loc=1)
        fig.add_artist(anchored_text)

        plt.grid(True)
        plt.draw()

        # now plot the flow vs speed
        fig = plt.figure(figsize=(10, 7), dpi=100)
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        plt.scatter(ff_speed, ff_flow, color='g')
        plt.grid(True)
        plt.hold(True)

        plt.scatter(cg_speed, cg_flow, color='k')
        plt.grid(True)
        plt.xlabel('Speed (mph)')
        plt.ylabel('Flow (veh/hr)')
        plt.title('Flow speed sensor, agg {0} s'.format(agg_interval))

        plt.draw()

    # deprecated
    def run_fwdsim_EnKF(self):
        """
        This function runs the EnKF forward simulation to see if the model is propagating properly.
        :return:
        """
        # run the estimators on each replication dataset
        for rep in self.active_replications:

            print('Status: Running forward simulation... ')

            for config_id in self.all_config.keys():

                # a configuration may potentially have multiple algorithms
                for alg_id in self.all_config[config_id]['algorithms'].keys():

                    # get the data for this configuration
                    flow_data, speed_data, traveltime_data, sensors = self.__prepare_data_for_estimation(config_id, rep)

                    # overwrite the data
                    flow_data = None
                    speed_data = None
                    traveltime_data = None

                    if self.all_config[config_id]['algorithms'][alg_id]['type'] == 'enkf':
                        # print 'Status: running estimator {0}...'.format(alg_id)
                        time0 = time.time()

                        estimator = self.__run_EnKF_estimator(flow_data, speed_data, traveltime_data, sensors,
                                                              self.all_config[config_id]['algorithms'][alg_id])
                        time1 = time.time()
                        print('      --- Took {0:.2f} s.'.format(time1 - time0))

                    # save the generated result into corresponding result_dir, and update result_record_file
                    result_name_prefix = self.__truestate_result_dir[rep] + 'fwdsim_enkf'
                    # save the transposed density and speed: num_steps x num_cells on the finest grid
                    np.savetxt(result_name_prefix + "_density.txt", estimator.est_density.T, delimiter=",")
                    np.savetxt(result_name_prefix + "_speed.txt", estimator.est_speed.T, delimiter=",")

                    np.savetxt(result_name_prefix + "_queue.txt",
                               np.matrix(estimator.est_queue).reshape((len(self.time_grid) - 1, 1)),
                               delimiter=",")
                    np.savetxt(result_name_prefix + "_traveltime.txt",
                               np.matrix(estimator.est_traveltime).reshape((len(self.time_grid) - 1, 1)),
                               delimiter=",")

                    # update the configuration_log
                    with open(self.__result_log[rep], 'a+') as f_log:
                        f_log.write('{0}&{1}\n'.format(config_id, alg_id))
                    f_log.close()

    def __weight_samples(self, density_array, speed_array, flow_array, density_range=(0, 560), num_bins=56):
        """
        In FD calibration, the collected density and speed data may cluster to a certain density, resulting in inaccurate
        calibration of the FD. For instance, if 49% of freeflow data points are clustering around (50veh/mile, 60 mph)
        and 49% of congested points are around (400 veh/mile, 2.5 mph), then the line fitting will be terrible, e.g.,
        would not capture the other points
        This function investigate the distribution of the points and re-weight the points by adding more points to less
        represented density bins.
        :param density_array: the density data point, array
        :param flow_array: the corresponding speed data points, array
        :param density_range: the range of the density
        :param num_bins: number of bins for distribution
        :return: the weighted samples
        """
        bins = np.linspace(density_range[0], density_range[1], num_bins + 1)
        num_pts = np.zeros(num_bins)

        for b in range(0, num_bins):
            num_pts[b] = np.sum((density_array >= bins[b]) & (density_array < bins[b + 1]))

        # get the maximum number of point in bins
        max_num_pts = int(np.max(num_pts))

        # Save data points in array
        density_array_weighted = np.array([])
        speed_array_weighted = np.array([])
        flow_array_weighted = np.array([])
        for b in range(0, num_bins):

            pt_idx_in_bin = (density_array >= bins[b]) & (density_array < bins[b + 1])

            # e.g. expand 5000 pt to 1300 pt;
            # first duplicate 5000 to 10000 pt
            for i in range(0, max_num_pts / int(num_pts[b])):
                density_array_weighted = np.concatenate([density_array_weighted,
                                                         density_array[pt_idx_in_bin]])
                speed_array_weighted = np.concatenate([speed_array_weighted,
                                                       speed_array[pt_idx_in_bin]])
                flow_array_weighted = np.concatenate([flow_array_weighted,
                                                      flow_array[pt_idx_in_bin]])

            # for the rest 3000 pt, random select 3000 pionts
            rand_idx = random.sample(range(0, len(density_array)), max_num_pts % int(num_pts[b]))
            density_array_weighted = np.concatenate([density_array_weighted,
                                                     density_array[rand_idx]])
            speed_array_weighted = np.concatenate([speed_array_weighted,
                                                   speed_array[rand_idx]])
            flow_array_weighted = np.concatenate([flow_array_weighted,
                                                  flow_array[rand_idx]])

        return density_array_weighted, speed_array_weighted, flow_array_weighted

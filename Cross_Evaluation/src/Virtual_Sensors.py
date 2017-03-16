import numpy as np
import sys
import os
from os.path import exists
import time
import math
from collections import OrderedDict
import random
import warnings
from scipy.interpolate import InterpolatedUnivariateSpline
warnings.filterwarnings("error")

"""
This class is used to generate the true states, and the virtual sensor data.
methods:
    generate_virtual_sensor_data
    generate_true_states_data
"""
__author__ = 'Juan Carlos Martinez'


class Virtual_Sensors:

    def __init__(self, work_zone, log_dir, data_dir, sections, main_grid_order, replications,
                 aimsun_start_dur_step=(55800,9000,0.8), space_start_end=None):
        """

        :param work_zone: str, 'I57', 'I80'
        :param log_dir: the top level directory of the logs.
        :param data_dir: the top level directory of the actual data, default E:\\Workzone\\
        :param sections:
        :param main_grid_order:
        :param aimsun_start_dur_step: [start, dur, step ] the start timestamp, duration and step length in seconds
        :param space_start_end: the absolute start and end location of the space domain under investigation
        :return:
        """

        print('\n')
        print('Constructing virtual_sensors class...')
        
        self.__work_zone = work_zone
        self.__sections = sections
        self.__main_grid_order = main_grid_order
        self.__space_start = space_start_end[0]
        self.__space_end = space_start_end[1]
        self.replications = replications
        self.__num_lanes = 2  # Hard coded

        self.__start_timestamp = aimsun_start_dur_step[0]
        self.__end_timestamp = self.__start_timestamp + aimsun_start_dur_step[1]
        
        # default parameters for the replication and the sensor types
        from default_parameters import default_parameters
        self.__default_parameters = default_parameters
        self.__default_parameters['time_step'] = aimsun_start_dur_step[2]

        # Make the directories a dict with keys as the replication id
        self.__trajectories_file = OrderedDict()
        self.__virtual_sensors_data_folder = OrderedDict()
        self.__true_state_data_folder = OrderedDict()
        self.__virtual_sensors_to_generate_file = OrderedDict()
        self.__virtual_sensors_generated_file = OrderedDict()

        if sys.platform == 'win32':
            for rep in self.replications:
                self.__trajectories_file[rep] = data_dir + 'Trajectory_data\\' + '{0}_rep{1}.csv'.format(self.__work_zone, rep)
                self.__virtual_sensors_data_folder[rep] = data_dir + 'Virtual_sensor_data\\' + self.__work_zone + '\\' + 'rep{0}'.format(rep) + '\\'
                self.__true_state_data_folder[rep] = data_dir + 'Estimation_results\\{0}\\rep{1}\\truestate\\'.format(self.__work_zone, rep)
                self.__virtual_sensors_to_generate_file[rep] = log_dir + 'Virtual_sensor_data\\' +  self.__work_zone + '_' + 'rep{0}'.format(rep)  + '_to_generate.txt'
                self.__virtual_sensors_generated_file[rep] =  log_dir + 'Virtual_sensor_data\\' +  self.__work_zone + '_' + 'rep{0}'.format(rep)  + '_generated.txt'
        elif sys.platform == 'darwin' or sys.platform == 'linux2':
            for rep in self.replications:
                self.__trajectories_file[rep] = data_dir + 'Trajectory_data/' + '{0}_rep{1}.csv'.format(self.__work_zone, rep)
                self.__virtual_sensors_data_folder[rep] = data_dir + 'Virtual_sensor_data/' + self.__work_zone + '/' + 'rep{0}'.format(rep) + '/'
                self.__true_state_data_folder[rep] = data_dir + 'Estimation_results/{0}/rep{1}/truestate/'.format(self.__work_zone, rep)
                self.__virtual_sensors_to_generate_file[rep] = log_dir + 'Virtual_sensor_data/' +  self.__work_zone +'_'+ 'rep{0}'.format(rep)  + '_to_generate.txt'
                self.__virtual_sensors_generated_file[rep] =  log_dir + 'Virtual_sensor_data/' +  self.__work_zone + '_' + 'rep{0}'.format(rep)  + '_generated.txt'
        else:
            raise Exception('Error: run "import sys; sys.platform" to check the platform. ')

        # do Not load data until generating virtual sensor data and the true states.
        self.__raw_data = None
        self.__organized_data = None

    def generate_virtual_sensors_data(self):
        """
        this is a public function to generate the sensor data specified when class
        was constructed
        :return:
        """

        time0 = time.time()
        print('Status: Generating virtual sensors data...')

        for rep in self.replications:

            if not exists(self.__virtual_sensors_to_generate_file[rep]):
                print('Status: no additional sensor data need to be generated for Replication {0}.'.format(rep))
                continue

            print('Status: generating sensor data for replication {0}...'.format(rep))

            if self.__organized_data is None:
                time0 = time.time()
                print('Status: Loading and organizing trajectory data...')
                self.__load_and_organize()
                time1 = time.time()
                print('Status: Loaded and organized trajectory data...')
                elapsed_time = '{0:.2f}'.format(time0 - time1)
                print('Elapsed time: ' + elapsed_time + ' s')

            with open(self.__virtual_sensors_to_generate_file[rep]) as sensors_configuration:

                print('Status: opened sensor to generate file:\n {0}'.format(
                    self.__virtual_sensors_to_generate_file[rep]))

                for sensor in sensors_configuration:

                    sensor_parameters = self.__parse_sensor_line(sensor)
                    time_s = time.time()

                    # call appropriate helper functions
                    if (sensor_parameters['type'] == 'icone' or
                            sensor_parameters['type'] == 'radar' or
                            sensor_parameters['type'] == 'rtms' or
                            sensor_parameters['type'] == 'vol_gt'):
                        self.__build_volume_sensor(sensor_parameters)
                    elif (sensor_parameters['type'] == 'bluetooth' or
                            sensor_parameters['type'] == 'tt_gt'):
                        self.__build_travel_time_sensor(sensor_parameters)

                    time_e = time.time()
                    print('Status: -- generated {0} data. Took {1:.2f}.'.format(
                        sensor_parameters['id'], time_e - time_s))
                    
        time1 = time.time()
        elapsed_time = '{0:.2f}'.format(time1 - time0)
        print('Elapsed time: ' + elapsed_time + ' s')

    def generate_true_states_data(self, grid_resolution, rep, queue_threshold, plot=False):
        """
        this function generates the true states
        :param grid_resolution: (seconds, meters)
        :param rep: replication number
        :param queue_threshold: [m]
        :param plot: True or False
        :return:
        """

        if self.__organized_data is None:
            time0 = time.time()
            print('Status: Loading and organizing trajectory data...')
            self.__load_and_organize()
            time1 = time.time()
            print('Status: Loaded and organized trajectory data...')
            elapsed_time = '{0:.2f}'.format(time0 - time1)
            print('Elapsed time: ' + elapsed_time + ' s')

        time0 = time.time()
        print('\n')
        print('Generating true states data...')
        self.__build_true_flow_density_speed_states(grid_resolution, rep, queue_threshold)  # call for private function
        time1 = time.time()
        elapsed_time = '{0:.2f}'.format(time1 - time0)
        print('Elapsed time: ' + elapsed_time + ' s')

    def __parse_sensor_line(self, sensor):
        """
        this function parses the sensor file whose data need to be generated
        :param sensor: the sensor line
        :return:
        """
        # parse, load and update sensor parameters
        sensor_parameters = OrderedDict()

        # custom parameters (different from default)
        sensor_parameters['custom'] = OrderedDict()
        for couple in sensor.split(';'):
            category = couple.split(':', 1)[0].rstrip()
            value = couple.split(':', 1)[1].rstrip()
            if value[0] == '[':
                value = [float(i) for i in value.strip('[]').split(',')]
            elif value.find('.') != -1:
                value = float(value)
            elif value.isdigit():
                value = int(value)
            elif value == 'True':   # the occlusion
                value = True
            elif value == 'False':
                value = False

            # parameters that are specific to the sensor type
            if (category != 'type' and category != 'id' and category != 'section' and category != 'distance' and
                    category != 'section_up' and category != 'section_down' and category != 'distance_up' and
                    category != 'distance_down'):
                sensor_parameters['custom'][category] = value
            # parameters that all sensors have
            else:
                sensor_parameters[category] = value

        # update 'default' sensor parameters for radar, rtms, and ideal sensors.
        sensor_parameters['default'] = self.__default_parameters[sensor_parameters['type']].copy()
        if sensor_parameters['custom']:
            for category in sensor_parameters['custom']:
                sensor_parameters['default'][category] = sensor_parameters['custom'][category]

        return self.__compute_sensor_offsets(sensor_parameters)

    def __compute_sensor_offsets(self, sensor_parameters):
        """
        this function computes the offsets of a sensor based on absolute distances
        :param sensor_parameters: dictionary holding sensor parameters
        :return:
        """

        # compute absolute device location (offset),
        sensor_parameters['device_abs_location'] = (self.__sections[sensor_parameters['section']]['offset'] +
                                                    sensor_parameters['distance'] +
                                                    sensor_parameters['default']['offset_dist'])
        # compute device section
        section_idx = self.__main_grid_order.index(sensor_parameters['section'])
        while (sensor_parameters['device_abs_location'] >
                (self.__sections[self.__main_grid_order[section_idx]]['offset'] +
                 self.__sections[self.__main_grid_order[section_idx]]['length'])):
            section_idx += 1
        sensor_parameters['device_section'] = self.__main_grid_order[section_idx]

        # compute sensing area boundaries for each lane
        l = sensor_parameters['default']['L']
        s = sensor_parameters['default']['s']
        alpha = np.radians(sensor_parameters['default']['alpha'])
        theta = np.radians(sensor_parameters['default']['theta'])
        sensor_parameters['sensing_area'] = OrderedDict()
        for lane_idx in range(1, self.__num_lanes + 1):
            sensor_parameters['sensing_area'][lane_idx] = OrderedDict()
            # offsets
            b = (s + lane_idx * l - l / 2) * (np.tan(alpha + theta) - np.tan(alpha))
            a = (s + lane_idx * l - l / 2) * (np.tan(alpha) - np.tan(alpha - theta))

            # section and absolute distance of upstream  and downstream boundary
            sensor_parameters['sensing_area'][lane_idx]['upstream'] = OrderedDict()
            sensor_parameters['sensing_area'][lane_idx]['upstream']['abs_distance'] \
                = sensor_parameters['device_abs_location'] - sensor_parameters['default']['offset_dist'] - b
            sensor_parameters['sensing_area'][lane_idx]['downstream'] = OrderedDict()
            sensor_parameters['sensing_area'][lane_idx]['downstream']['abs_distance'] \
                = sensor_parameters['device_abs_location'] - sensor_parameters['default']['offset_dist'] + a

            section_idx = self.__main_grid_order.index(sensor_parameters['section'])
            while (sensor_parameters['sensing_area'][lane_idx]['upstream']['abs_distance'] <
                    (self.__sections[self.__main_grid_order[section_idx]]['offset'])):
                section_idx -= 1
            sensor_parameters['sensing_area'][lane_idx]['upstream']['section'] = self.__main_grid_order[section_idx]

            section_idx = self.__main_grid_order.index(sensor_parameters['section'])
            while (sensor_parameters['sensing_area'][lane_idx]['downstream']['abs_distance'] >
                    (self.__sections[self.__main_grid_order[section_idx]]['offset'] +
                     self.__sections[self.__main_grid_order[section_idx]]['length'])):
                section_idx += 1
            sensor_parameters['sensing_area'][lane_idx]['downstream']['section'] = self.__main_grid_order[section_idx]

        return sensor_parameters

    def __set_section_offsets(self):
        """
        this function computes the offsets for each of the sections in the network
        :return:
        """

        # set the offset along the main grid (including tapers)
        total_offset = 0
        for section in self.__main_grid_order:
            self.__sections[section]['offset'] = total_offset
            total_offset += self.__sections[section]['length']

            # print('{0} Section offset: {1}'.format(section, self.__sections[section]['offset']))

        # set the offsets for the ramps
        for section in self.__sections:
            if 'connections' in self.__sections[section]:
                for up_section in self.__sections[section]['connections']['upstream']:
                    if 'offset' not in self.__sections[up_section]:
                        self.__sections[up_section]['offset'] = self.__sections[section]['offset'] - \
                                                                self.__sections[up_section]['length']
                for down_section in self.__sections[section]['connections']['downstream']:
                    if 'offset' not in self.__sections[down_section]:
                        self.__sections[down_section]['offset'] = self.__sections[section]['offset'] + \
                                                                  self.__sections[section]['length']

        # for section in self.__sections:
        #     print('offset of section {0}: {1}'.format(section, self.__sections[section]['offset']))

    # @profile
    def __load_and_organize(self):
        """
        this function loads and organized the trajectories data specified in the constructor
        remark: AIMSUN must be configured as 'metric' when generating the trajectory data file
        this function loads the trajectory data, with the columns as follows:
            did: replication id
            oid: vehicle id
            ent: counter of the vehicle trajectory
            sectionID: the section id of the current traj point
            lane index: the lane index
            xCoord: the x coordinate
            yCoord: the y coordinate
            timeSta: (s) the time stamp
            speed: the speed of the vehicle at the current point
            travelledDistance: the travelled distance of the current traj point
            acceleration: the acceleration rate
        note:
            when AIMSUN unit specified as 'metric', the units are:
             timeSta (s), speed (kph), travelledDistance (m); acceleration (?)
            when AIMSUN unit specified as 'imperial', the units are:
             timeSta (s), speed (mph), travelledDistance (m); acceleration (?)
        :return:
        """

        # all data will be organized here
        self.__organized_data = OrderedDict()

        # set section offsets based on main road
        self.__set_section_offsets()

        # keep a list of the sections that are nodes
        nodes_list = []
        for section in self.__sections:
            if 'connections' in self.__sections[section]:
                nodes_list.append(section)

        # initialize some variables for larger scope
        route_offset = None
        prev_section = None
        prev_line = None

        # iterate over every replication
        for rep in self.replications:

            # each replication is a key
            self.__organized_data[rep] = OrderedDict()

            # read the trajectories file line by line
            with open(self.__trajectories_file[rep], 'r', 1) as traj_file:
                for line in traj_file:
                    line = line.split(',')

                    section = int(line[3])
                    # indexes are as follows (variables not generated
                    # for better computation time)
                    # get data from line
                    # vehicle = int(line[1])
                    # ent = int(line[2])
                    # lane_idx = int(line[4])
                    # timestamp = float(line[7])
                    # speed = float(line[8])
                    # travelled_distance = float(line[9])
                    # acceleration = float(line[10])

                    # only consider segment under estimation
                    # if section not in self.__sections.keys():
                    #    continue

                    # each vehicle is a key
                    if int(line[1]) not in self.__organized_data[rep]:
                        self.__organized_data[rep][int(line[1])] = OrderedDict()
                        # offset for section distance of each vehicle
                        route_offset = 0.0
                        # keep track of previous section traversed by vehicle
                        prev_section = section

                    # if section is upstream of a node, reconstruct trajectory
                    # on the node
                    if section != prev_section:

                        # update route_offset
                        route_offset += self.__sections[prev_section]['length']

                        # find the node section being crossed
                        for node in nodes_list:
                            if (prev_section in self.__sections[node]['connections']['upstream'] and
                                    section in self.__sections[node]['connections']['downstream']):
                                node_section = node

                        # lane at node is assumed the same as upstream
                        # node_lane_idx = lane_idx

                        # get data from previous line
                        # prev_speed = float(prev_line[8])
                        # these are not directly referenced by prev_line because
                        # profiler results suggest this is faster
                        prev_timestamp = float(prev_line[7])
                        prev_travelled_distance = float(prev_line[9])

                        # compute length left to traverse at previous section
                        prev_section_remaining_length = route_offset - prev_travelled_distance

                        # compute general data for node trajectory reconstruction
                        # speed used for trajectory reconstruction
                        node_speed_m_s = (float(line[9]) - prev_travelled_distance)/(float(line[7]) - prev_timestamp)
                        # recorded speed
                        node_speed_km_h = (float(line[8]) + float(prev_line[8]))/2
                        # node_acceleration = 0.0
                        node_timestamps = np.arange(prev_timestamp + self.__default_parameters['time_step'],
                                                    float(line[7]), self.__default_parameters['time_step'])

                        # compute data for each timestamp reconstructed
                        for node_timestamp in node_timestamps:

                            node_timestamp = float('{0:.2f}'.format(node_timestamp))
                            # if difference between node_timestamp and timestamp is less than half time step,
                            # consider it to be a rounding error and pass
                            if abs(float(line[7]) - node_timestamp) < self.__default_parameters['time_step']/2:
                                pass
                            else:
                                # distance on node since last timestamp (if negative, the vehicle
                                # is still at the previous section)
                                node_section_distance = float('{0:.2f}'.format(
                                    node_speed_m_s*(node_timestamp - prev_timestamp) - prev_section_remaining_length))
                                # travelled distance at current timestamp
                                node_travelled_distance = float('{0:.2f}'.format(
                                    prev_travelled_distance + node_speed_m_s*(node_timestamp - prev_timestamp)))

                                if node_section_distance < 0:
                                    # vehicle is still on previous section
                                    selected_section = prev_section
                                    selected_section_distance = float('{0:.2f}'.format(
                                        self.__sections[prev_section]['length'] + node_section_distance))
                                else:
                                    # vehicle is in node
                                    selected_section = node_section
                                    selected_section_distance = float('{0:.2f}'.format(node_section_distance))

                                # main offset based on selected section
                                main_offset = selected_section_distance + self.__sections[selected_section]['offset']

                                # record reconstructed data
                                # -1 is set as the ent identifier for reconstructed data
                                # [ent, section, section_distance, travelled_distance, lane_idx, speed,
                                # acceleration, main_offset]
                                self.__organized_data[rep][int(line[1])][node_timestamp] \
                                    = [-1, selected_section, selected_section_distance, node_travelled_distance,
                                       int(line[4]), node_speed_km_h, 0.0, main_offset]

                        # update route_offset
                        route_offset += self.__sections[node_section]['length']
                        # update previous section
                        prev_section = section

                    # compute section distance based on the vehicle route offset
                    section_distance = float(line[9]) - route_offset

                    # section distances over the section length are truncated to avoid
                    # possible inconsistencies
                    if section_distance > self.__sections[section]['length']:
                        section_distance = self.__sections[section]['length']

                    # offset with respect to main road and its origin
                    # main_offset = section_distance + self.__sections[section]['offset']

                    # record data
                    # [ent, section, section_distance, travelled_distance, lane_idx, speed, acceleration, main_offset]
                    self.__organized_data[rep][int(line[1])][float(line[7])] = \
                        [int(line[2]), section, section_distance, float(line[9]), int(line[4]),
                         float(line[8]), float(line[10]), section_distance + self.__sections[section]['offset']]
                    # keep track of previous line
                    prev_line = line

    # @profile
    def __build_volume_sensor(self, sensor_parameters):
        """
        this is a private function that builds a volume sensor,
        saves the data and logs the sensor
        :param sensor_parameters: dictionary holding sensor parameters
        :return:
        """

        # iterate over every replication
        for replication in self.__organized_data:
            
            # build crosses and readings
            sensor_crosses = self.__build_sensor_crosses(sensor_parameters, replication)
            sensor_readings = self.__build_volume_sensor_readings(sensor_crosses, sensor_parameters)

            # set sensor file
            sensor_file = self.__virtual_sensors_data_folder[replication] + str(sensor_parameters['id']) + '.txt'
            if not os.path.exists(os.path.dirname(sensor_file)):
                os.makedirs(os.path.dirname(sensor_file))

            # save virtual sensor data
            with open(sensor_file, 'w+') as output_file:
                for interval in sensor_readings:
                    line = str([interval, float('%.2f' % sensor_readings[interval]['sms']),
                                sensor_readings[interval]['count']]).strip('[]').replace(' ', '')
                    output_file.write(line + '\n')
            # log virtual sensor
            with open(self.__virtual_sensors_generated_file[replication], 'a+') as output_file:
                line = self.__sensordict2string(str(sensor_parameters['id']), sensor_parameters)
                output_file.write(line + '\n')

    # @profile
    def __build_sensor_crosses(self, sensor_parameters, replication):
        """
        this function computes the sensor crosses matrix
        :param sensor_parameters: dictionary containing the parameters of the sensor being crossed
        :param replication: replication id
        :return sensor_crosses: matrix of sensor crosses
        """

        # stack that will hold crosses information
        sensor_crosses = np.atleast_2d(np.empty((0, 5)))

        # extract road and sensor information
        half_veh_length = self.__default_parameters['veh_length']/2

        # iterate over every vehicle
        for vehicle in self.__organized_data[replication]:
            # binary search
            timestamps = list(self.__organized_data[replication][vehicle].keys())
            vehicle_data = list(self.__organized_data[replication][vehicle].values())
            min_idx = 0
            max_idx = len(vehicle_data) - 1

            # discard if vehicle never reaches sensor location (at max l
            if (sensor_parameters['sensing_area'][self.__num_lanes]['upstream']['abs_distance'] <
                    vehicle_data[min_idx][7] or
                    sensor_parameters['sensing_area'][self.__num_lanes]['downstream']['abs_distance'] >
                    vehicle_data[max_idx][7]):
                continue

            # used to break from while loop
            # vehicle can only cross once
            crossed = False

            # loop breaks when max_idx and min_idx are adjacent
            while (max_idx - min_idx) > 1 and not crossed:

                # get info from middle index
                mid_idx = (max_idx - min_idx)//2 + min_idx
                lane_idx = vehicle_data[mid_idx][4]
                offset_distance = vehicle_data[mid_idx][7]

                # not sure what this does (ask Yanning)
                if (min_idx + 1 == mid_idx and mid_idx + 1 == max_idx and
                    (vehicle_data[min_idx][4] != vehicle_data[mid_idx][4] or
                     vehicle_data[mid_idx][4] != vehicle_data[max_idx][4])):
                    break

                # check if the upstream sensing boundary has been crossed
                upstream_cross_before = (offset_distance + half_veh_length <
                                         sensor_parameters['sensing_area'][lane_idx]['upstream']['abs_distance'])
                upstream_cross_after = (sensor_parameters['sensing_area'][lane_idx]['upstream']['abs_distance'] <=
                                        vehicle_data[mid_idx + 1][7] + half_veh_length)

                # if upstream sensing boundary was crossed, get crossing data and
                # iterate until downstream sensing boundary is crossed
                if upstream_cross_before and upstream_cross_after:

                    # discard cars that switched lanes right on sensor boundaries
                    if lane_idx != vehicle_data[mid_idx + 1][4]:
                        break

                    # check that vehicle section is same as sensor section
                    # section = vehicle_data[mid_idx][1]
                    if (vehicle_data[mid_idx][1] !=
                            sensor_parameters['sensing_area'][lane_idx]['upstream']['section'] and
                        vehicle_data[mid_idx+1][1] !=
                            sensor_parameters['sensing_area'][lane_idx]['downstream']['section']):
                        break

                    timestamp = timestamps[mid_idx]
                    distance_to_upstream = (sensor_parameters['sensing_area'][lane_idx]['upstream']['abs_distance'] -
                                            offset_distance - half_veh_length)
                    cross_speed_m_s = (vehicle_data[mid_idx + 1][7] - offset_distance)/(timestamps[mid_idx + 1] -
                                                                                        timestamp)

                    enter_timestamp = float('{0:.2f}'.format(timestamp + distance_to_upstream/cross_speed_m_s))
                    enter_speed_km_h = vehicle_data[mid_idx][5]
                    enter_lane_idx = lane_idx

                    for idx in range(mid_idx, len(timestamps)-1):

                        try:
                            lane_idx = vehicle_data[idx][4]
                            # offset_distance = vehicle_data[idx][7]

                            downstream_cross_before = (vehicle_data[idx][7] - half_veh_length <
                                                       sensor_parameters['sensing_area'][lane_idx]['downstream'][
                                                           'abs_distance'])
                            downstream_cross_after = (sensor_parameters['sensing_area'][lane_idx]['downstream'][
                                                           'abs_distance'] <=
                                                      vehicle_data[idx + 1][7] - half_veh_length)

                            # if downstream sensing boundary was crossed, get crossing data
                            if downstream_cross_before and downstream_cross_after:

                                # timestamp = timestamps[idx]
                                distance_to_downstream = (sensor_parameters['sensing_area'][lane_idx]['downstream'][
                                                           'abs_distance'] -
                                                          vehicle_data[idx][7] + half_veh_length)
                                cross_speed_m_s = ((vehicle_data[idx + 1][7] - vehicle_data[idx][7]) /
                                                   (timestamps[idx + 1] - timestamps[idx]))
                                exit_timestamp = float('{0:.2f}'.format(timestamps[idx] +
                                                                        distance_to_downstream/cross_speed_m_s))

                                # exit_speed_km_h = vehicle_data[idx][5]
                                # exit_lane_idx = lane_idx

                                recorded_speed = (enter_speed_km_h + vehicle_data[idx][5])/2
                                recorded_lane_idx = min(enter_lane_idx, lane_idx)

                                cross_data = [enter_timestamp, exit_timestamp, vehicle, recorded_lane_idx,
                                              recorded_speed]

                                sensor_crosses = np.vstack((sensor_crosses, cross_data))
                                # to break from while loop
                                crossed = True
                                break

                            # if there is a lane switch at downstream edge, consider that car crossed but
                            # do not record data (car was effectively missed)
                            elif (not downstream_cross_before) and downstream_cross_after:
                                # to break from while loop
                                crossed = True
                                break

                        # to avoid indexing past end of list
                        # car is discarded (didn't actually finish the crossing)
                        except IndexError:
                            crossed = True
                            break


                    crossed = True
                    break

                elif upstream_cross_before:
                    min_idx = mid_idx
                else:
                    max_idx = mid_idx

        return sensor_crosses

    def __build_volume_sensor_readings(self, crosses, sensor_parameters):
        """
        this function takes the sensor crossings and computes sensor readings
        :param crosses: matrix of sensor crossings
        :param sensor_parameters: dictionary of sensor paramaters
        :return:
        """

        # obtain noise parameters
        random.seed()
        noise_type, p_occlusion_accept, occlusion_config = \
            sensor_parameters['default']['noise_type'], \
            sensor_parameters['default']['p_occlusion_accept'], \
            sensor_parameters['default']['occlusion']

        aggregation_sec, awake_sec = sensor_parameters['default']['aggregation_sec'], sensor_parameters['default']['awake_sec']
        v_range, v_threshold = sensor_parameters['default']['v_range'], sensor_parameters['default']['v_threshold']
        p_missing_ff, p_missing_cf = sensor_parameters['default']['p_missing_ff'], sensor_parameters['default']['p_missing_cf']
        if noise_type == 'relative':
            v_bias_ff, v_bias_cf = sensor_parameters['default']['v_bias_ff'], sensor_parameters['default']['v_bias_cf']
            v_accuracy_p_ff, v_accuracy_p_cf = sensor_parameters['default']['v_accuracy_p_ff'], sensor_parameters['default']['v_accuracy_p_cf']
        elif noise_type == 'absolute':
            v_noise_sigma_ff, v_noise_sigma_cf = sensor_parameters['default']['v_noise_sigma_ff'], sensor_parameters['default']['v_noise_sigma_cf']
            v_noise_mu_ff, v_noise_mu_cf = sensor_parameters['default']['v_noise_mu_ff'], sensor_parameters['default']['v_noise_mu_cf']


        curr_r_occlusion_accept = random.uniform(p_occlusion_accept[0], p_occlusion_accept[1])/100
        curr_r_missing_ff, curr_r_missing_cf = random.uniform(p_missing_ff[0], p_missing_ff[1])/100, random.uniform(p_missing_cf[0], p_missing_cf[1])/100

        if noise_type == 'relative':
            curr_v_bias_ff, curr_v_bias_cf = random.uniform(v_bias_ff[0], v_bias_ff[1]), random.uniform(v_bias_cf[0], v_bias_cf[1])
            curr_v_accuracy_r_ff, curr_v_accuracy_r_cf = random.uniform(v_accuracy_p_ff[0], v_accuracy_p_ff[1])/100, random.uniform(v_accuracy_p_cf[0], v_accuracy_p_cf[1])/100
        if noise_type == 'absolute':
            curr_v_noise_sigma_ff, curr_v_noise_sigma_cf = random.uniform(v_noise_sigma_ff[0], v_noise_sigma_ff[1]), random.uniform(v_noise_sigma_cf[0], v_noise_sigma_cf[1])
            curr_v_noise_mu_ff, curr_v_noise_mu_cf = random.uniform(v_noise_mu_ff[0], v_noise_mu_ff[1]), random.uniform(v_noise_mu_cf[0], v_noise_mu_cf[1])

        sensor = OrderedDict()
        intervals = aggregation_sec*(np.arange(math.floor(self.__start_timestamp/aggregation_sec),math.floor(self.__end_timestamp/aggregation_sec + 1)) -
                                     math.floor(self.__start_timestamp/aggregation_sec))

        for interval in intervals:
            sensor[interval] = OrderedDict()
            sensor[interval]['crossed?'] = False
            sensor[interval]['inverse_speed_sum'] = 0
            sensor[interval]['count'] = 0
            sensor[interval]['flow'] = 0
            sensor[interval]['sms'] = float('inf')
            sensor[interval]['density'] = 0

        if crosses.any():
            max_lane_idx = int(np.max(crosses[:,3]))
            # iterate over every lane_idx available
            for lane_idx in range (1, max_lane_idx + 1):
                for veh_cross in crosses[crosses[:,3] == float(lane_idx), :]:
                    # the interval is set based on the enter_timestamp
                    enter_timestamp = veh_cross[0]
                    if (enter_timestamp % aggregation_sec) <= awake_sec:
                        # check for occlusion
                        clear_vision = True
                        if lane_idx != 1 and curr_r_occlusion_accept < 1:
                            clear_vision = self.__occlusion_check(crosses, veh_cross, curr_r_occlusion_accept, occlusion_config)
                        if clear_vision:
                            interval = aggregation_sec*(math.floor(enter_timestamp/aggregation_sec) -
                                                        math.floor(self.__start_timestamp/aggregation_sec))
                            if not sensor[interval]['crossed?']:
                                sensor[interval]['crossed?'] = True
                            speed = veh_cross[4]
                            # add noise
                            if noise_type == 'relative':
                                if speed >= v_threshold:
                                    r_missing, v_bias, v_accuracy_r = curr_r_missing_ff, \
                                                                      curr_v_bias_ff, curr_v_accuracy_r_ff
                                else:
                                    r_missing, v_bias, v_accuracy_r = curr_r_missing_cf, \
                                                                      curr_v_bias_cf, curr_v_accuracy_r_cf

                                # noisy_speed = abs(random.gauss(speed + v_bias,(1 - v_accuracy_r)*speed))
                                noisy_speed = abs(random.gauss(speed + v_bias, v_accuracy_r*speed))
                            elif noise_type == 'absolute':
                                if speed >= v_threshold:
                                    r_missing, v_noise_mu, v_noise_sigma = curr_r_missing_ff, \
                                                                           curr_v_noise_mu_ff, curr_v_noise_sigma_ff
                                else:
                                    r_missing, v_noise_mu, v_noise_sigma = curr_r_missing_cf, \
                                                                           curr_v_noise_mu_cf, curr_v_noise_sigma_cf
                                noisy_speed = abs(speed + random.gauss(v_noise_mu,v_noise_sigma))
                            retain_data_check = random.uniform(0, 1)

                            # update data
                            if v_range[0] <= noisy_speed and noisy_speed <= v_range[1] and retain_data_check >= r_missing:
                                if noisy_speed:
                                    sensor[interval]['inverse_speed_sum'] += 1/noisy_speed
                                sensor[interval]['count'] += 1
                                sensor[interval]['flow'] = int(sensor[interval]['count']*3600/aggregation_sec)
                                if sensor[interval]['inverse_speed_sum']:
                                    # print( sensor[interval]['inverse_speed_sum'] )
                                    sensor[interval]['sms'] = float('{0:.10f}'.format(sensor[interval]['count']/sensor[interval]['inverse_speed_sum']))
                                    # print( sensor[interval]['sms'] )
                                    sensor[interval]['density'] = float('{0:.2f}'.format(sensor[interval]['flow']/sensor[interval]['sms']))

        # tell difference between zero count and no veh read because of drop rate
        for interval in sensor:
            if sensor[interval]['crossed?'] and not sensor[interval]['count']:
                sensor[interval]['count'] = -1
                sensor[interval]['flow'] = -1
                sensor[interval]['density'] = -1
                sensor[interval]['sms'] = -1

        # drop the data based on the data missing parameter
        for interval in sensor:
            # if freeflow
            if sensor[interval]['sms'] != -1:

                if sensor[interval]['sms'] >= v_threshold:
                    r_missing = curr_r_missing_ff
                    count_bias = random.uniform(sensor_parameters['default']['c_bias_p_ff'][0],
                                                sensor_parameters['default']['c_bias_p_ff'][1])/100.0
                    count_std = random.uniform(sensor_parameters['default']['c_sigma_p_ff'][0],
                                               sensor_parameters['default']['c_sigma_p_ff'][1])/100.0
                # if congested flow
                else:
                    r_missing = curr_r_missing_cf
                    count_bias = random.uniform(sensor_parameters['default']['c_bias_p_cf'][0],
                                                sensor_parameters['default']['c_bias_p_cf'][1])/100.0
                    count_std = random.uniform(sensor_parameters['default']['c_sigma_p_cf'][0],
                                               sensor_parameters['default']['c_sigma_p_cf'][1])/100.0

                retain_data_check = random.uniform(0,1)
                if retain_data_check <= r_missing:
                    # throw away data
                    sensor[interval]['sms'] = -1
                    sensor[interval]['count'] = -1
                    sensor[interval]['flow'] = -1
                    sensor[interval]['density'] = -1

                # add count bias and noise into the data set
                if sensor[interval]['sms'] != -1:

                    tmp_count = random.gauss(sensor[interval]['count']*(1+count_bias),
                                             sensor[interval]['count']*count_std) 
                    # put a zero threshold and update the value
                    tmp_count = np.max([tmp_count, 0])

                    sensor[interval]['count'] = tmp_count
                    sensor[interval]['flow'] = sensor[interval]['count']*3600/aggregation_sec

        return sensor

    @staticmethod
    def __occlusion_check(crosses, veh_cross, curr_r_occlusion_accept, occlusion_config):
        """
        this function determines if a vehicle crossing is occluded
        :param crosses: crosses matrix
        :param veh_cross: vehicle cross
        :param curr_r_occlusion_accept: occlusion acceptance
        :param occlusion_config: occlusion configuration
        :return boolean: True if clear vision, False otherwise
        """

        clear_vision = True

        if occlusion_config is True:

            lane_idx = int(veh_cross[3])
            veh_cross_time = veh_cross[1] - veh_cross[0]
            occlusion_acceptance = curr_r_occlusion_accept*veh_cross_time
            # Check if occlusion occurs on any of the inner lanes
            for inner_lane_idx in range(1, lane_idx):
                # Update the clear_vision boolean by checking that there is no overlap on any of the lanes. This happens when:
                # current vehicle exit_time < enter_time of any vehicle with a lower lane_idx (e2 < s1) OR
                # current vehicle enter_time > exit_time of any vehicle with a lower lane_index (e1 < s2)
                clear_vision = clear_vision and np.logical_or(veh_cross[1] < crosses[crosses[:,3] == float(inner_lane_idx),0] - occlusion_acceptance,
                                                              veh_cross[0] > crosses[crosses[:,3] == float(inner_lane_idx),1] + occlusion_acceptance
                                                              ).all()
        return clear_vision

    def __build_travel_time_sensor(self, sensor_parameters, save2file=True, agg_sec=float('nan'), true_state=False):
        """
        this function builds a travel time sensor pair
        :param sensor_parameters: sensor pair parameters
        :param save2file: True if used for virtual sensors, False if used for true states
        :return:
        """

        # get paramters from both sensors (upstream and downstream)
        up_sensor_parameters = sensor_parameters.copy()
        up_sensor_parameters['section'] = sensor_parameters['section_up']
        up_sensor_parameters['distance'] = sensor_parameters['distance_up']
        up_sensor_parameters = self.__compute_sensor_offsets(up_sensor_parameters)

        down_sensor_parameters = sensor_parameters.copy()
        down_sensor_parameters['section'] = sensor_parameters['section_down']
        down_sensor_parameters['distance'] = sensor_parameters['distance_down']
        down_sensor_parameters = self.__compute_sensor_offsets(down_sensor_parameters)

        all_pair_readings = OrderedDict()

        # iterate
        for replication in self.__organized_data:

            # build crosses at upstream and downstream sensors
            sensor_up_crosses = self.__build_sensor_crosses(up_sensor_parameters, replication)
            sensor_down_crosses = self.__build_sensor_crosses(down_sensor_parameters, replication)

            # build the sensor readings based on the sensor crossings
            pair_readings = self.__build_travel_time_sensor_readings(sensor_up_crosses, sensor_down_crosses,
                                                                     sensor_parameters, replication, true_state)

            # used for virtual sensors
            if save2file is True:

                sensor_file = self.__virtual_sensors_data_folder[replication] + str(sensor_parameters['id']) + '.txt'

                if not os.path.exists(os.path.dirname(sensor_file)):
                    os.makedirs(os.path.dirname(sensor_file))

                with open(sensor_file, 'w+') as output_file:
                    for interval in pair_readings:
                        # if travel time measurement exist
                        if ('travel_time' in pair_readings[interval].keys() and not
                                np.isnan(pair_readings[interval]['travel_time'])):
                            value = pair_readings[interval]['travel_time']
                        else:
                            value = np.nan
                        line = ','.join([str(i) for i in [interval, value]])
                        output_file.write(line + '\n')

                with open(self.__virtual_sensors_generated_file[replication], 'a+') as output_file:
                    line = self.__sensordict2string(str(sensor_parameters['id']), sensor_parameters)
                    output_file.write(line + '\n')

            # used for true states
            else:
                travel_time = []
                latest_reading = np.nan
                intervals = list(pair_readings.keys())
                for interval in intervals:
                    # add interval if there are readings
                    if not np.isnan(pair_readings[interval]['travel_time']):
                        # update the latest reading
                        latest_reading = pair_readings[interval]['travel_time']
                        travel_time.append(latest_reading)
                    else:
                        # no update
                        travel_time.append(latest_reading)
                # keep in dictionary
                all_pair_readings[replication] = {}
                all_pair_readings[replication]['tt'] = travel_time
                all_pair_readings[replication]['tt_beg_idx'] = intervals[0]//agg_sec

        return all_pair_readings

    def __build_travel_time_sensor_readings(self, sensor_up_crosses, sensor_down_crosses, sensor_parameters,
                                            replication, true_state):
        """
        compute the travel time sensor readings based on the sensor crosses upstream and downstream
        :param sensor_up_crosses:
        :param sensor_down_crosses:
        :param sensor_parameters:
        :param replication:
        :param true_state: boolean that determines if the sensors are true state or not
        :return:
        """

        # obtain noise parameters
        random.seed()
        aggregation_sec = sensor_parameters['default']['aggregation_sec']
        p_penetration = sensor_parameters['default']['p_penetration']
        t_noise_sigma = sensor_parameters['default']['t_noise_sigma']
        t_noise_mu = sensor_parameters['default']['t_noise_mu']
        t_sigma, t_mu = random.uniform(t_noise_sigma[0],t_noise_sigma[1]), random.uniform(t_noise_mu[0], t_noise_mu[1])
        curr_r_penetration = random.uniform(p_penetration[0],p_penetration[1])/100.0
        
        pair = OrderedDict()
        # add the travel time of each penetrated vehicle
        for vehicle in self.__organized_data[replication]:
            if random.uniform(0, 1) <= curr_r_penetration:
                up_timestamp = sensor_up_crosses[sensor_up_crosses[:, 2] == vehicle, 0]
                down_timestamp = sensor_down_crosses[sensor_down_crosses[:, 2] == vehicle, 0]
                if up_timestamp and down_timestamp:
                    noisy_up_timestamp = np.max([up_timestamp+random.gauss(t_mu, t_sigma), self.__start_timestamp])
                    noisy_down_timestamp = np.max([down_timestamp+random.gauss(t_mu, t_sigma), self.__start_timestamp])

                    # definition of travel time is different for true state and for measurements
                    if true_state:
                        interval = aggregation_sec*(math.floor(noisy_up_timestamp/aggregation_sec) -
                                                    math.floor(self.__start_timestamp/aggregation_sec))
                    else:
                        interval = aggregation_sec * (math.floor(noisy_down_timestamp / aggregation_sec) -
                                                      math.floor(self.__start_timestamp / aggregation_sec))

                    if interval not in pair:
                        pair[interval] = OrderedDict()
                        pair[interval]['up_times_list'] = []
                        pair[interval]['down_times_list'] = []
                    pair[interval]['up_times_list'].append(noisy_up_timestamp)
                    pair[interval]['down_times_list'].append(noisy_down_timestamp)

        # compute the means of vehicle travel times for each interval
        for interval in pair:
            # return the mean travel time
            if len(pair[interval]['up_times_list']):
                pair[interval]['travel_time'] = float('%.2f' % (np.mean( np.array(pair[interval]['down_times_list']) -
                                                                         np.array(pair[interval]['up_times_list']))))
            else:
                # otherwise not travel time recorded
                pair[interval]['travel_time'] = np.nan

        return pair

    def __build_true_flow_density_speed_states(self, resolution, replication, queue_threshold):

        # this function generates matrices that hold the true traffic state on
        # a road network based on Edie's definitions
        # :param resolution: tuple [aggregation in seconds, aggregation distance in meters]
        # :param replication: replication number
        # :param plot: Boolean (True if plots are requested, False otherwise)
        # #:return:

        # extract resolution and compute cell area
        agg_sec = resolution[0]
        agg_dist = resolution[1]
        cell_area = agg_sec*agg_dist

        # get the space cell grids only for the section investigation
        # start at self.__space_start
        space_cell_bounds = np.arange(self.__space_start, self.__space_end + agg_dist, agg_dist)

        # get the time cell grids only for the time domain under investigation
        # start at time stamp 0
        time_cell_bounds = (np.arange(self.__start_timestamp, self.__end_timestamp + agg_sec, agg_sec) -
                            self.__start_timestamp)

        # initialize an empty matrix for the space sum and time sum on each cell
        # for n cell boundaries on an axis, there are n-1 cells on that axis
        spaces_sum_mat = np.zeros((len(space_cell_bounds)-1, (len(time_cell_bounds)-1)))
        times_sum_mat = np.zeros((len(space_cell_bounds)-1, (len(time_cell_bounds)-1)))

        # iterate over every vehicle
        for vehicle in self.__organized_data[replication]:

            # get timestamps, sections and distances in arrays
            # normalize timestamps
            timestamps = np.array(list(self.__organized_data[replication][vehicle].keys())) - self.__start_timestamp
            vehicle_data = list(self.__organized_data[replication][vehicle].values())
            sections = [i[1] for i in vehicle_data]
            abs_distances = np.array([i[7] for i in vehicle_data])

            # get idxs of the sections that match the sections on the main grid
            main_idxs = []
            for idx in range(0, len(sections)):
                if sections[idx] in self.__main_grid_order:
                    main_idxs.append(idx)
            main_idxs = np.array(main_idxs)

            # get space and time bounds based on sections visited
            if len(main_idxs):

                # extract timestamps and absolute distances
                main_timestamps = timestamps[[main_idxs]]
                main_abs_distances = abs_distances[[main_idxs]]

                # find idxs that are not duplicate
                keep_idxs = [0]
                for i in range(1,len(main_abs_distances)):
                    if main_abs_distances[i] - main_abs_distances[i-1] > 0.1:
                        keep_idxs.append(i)
                keep_idxs = np.array(keep_idxs)

                # extract timestamps and absolute that are not duplicate
                main_timestamps = main_timestamps[[keep_idxs]]
                main_abs_distances = main_abs_distances[[keep_idxs]]

                # find idxs that are not duplicate
                keep_idxs = [0]
                for i in range(1,len(main_timestamps)):
                    if main_timestamps[i] - main_timestamps[i-1] > 0.1:
                        keep_idxs.append(i)
                keep_idxs = np.array(keep_idxs)

                # extract timestamps and absolute that are not duplicate
                main_timestamps = main_timestamps[[keep_idxs]]
                main_abs_distances = main_abs_distances[[keep_idxs]]

                # obtain cell bounds that are relevant to current vehicle
                veh_space_cell_bounds = [i for i in space_cell_bounds if (main_abs_distances[-1] + agg_dist >= i >=
                                                                          main_abs_distances[0] - agg_dist)]
                veh_time_cell_bounds = [i for i in time_cell_bounds if (main_timestamps[-1] + agg_sec >= i >=
                                                                        main_timestamps[0] - agg_sec)]

                # continue only if vehicle visited area on time and space under study
                if len(veh_space_cell_bounds) and len(veh_time_cell_bounds):

                    # interpolate/extrapolate on space given time (extrapolation is done by default)
                    space_spl = InterpolatedUnivariateSpline(main_timestamps, main_abs_distances, k=1)
                    # interpolate/extrapolate on time given space (extrapolation is done by default)
                    time_spl = InterpolatedUnivariateSpline(main_abs_distances, main_timestamps, k=1)

                    # get tuples of crossing at time and space boundaries
                    time_at_space_bounds = self.__interpolate_space_crosses(time_spl, veh_space_cell_bounds)
                    space_at_time_bounds = self.__interpolate_time_crosses(space_spl, veh_time_cell_bounds)

                    # sort crosses and remove nans
                    crosses = (np.sort(time_at_space_bounds + space_at_time_bounds, axis=0)).tolist()
                    crosses = [cross for cross in crosses if (not math.isnan(cross[0]) and not math.isnan(cross[1]) and
                                                              cross[0] >= 0 and cross[1] >= 0)]
                    prev_cross = crosses[0]
                    # update values on spaces_sum_mat and times_sum_mat
                    for i in range(1, len(crosses)):
                        # break once last boundaries are crossed
                        if crosses[i][0] >= time_cell_bounds[-1] or crosses[i][1] >= space_cell_bounds[-1]:
                            break

                        # get values and indices to put in matrix
                        delta_t = crosses[i][0] - prev_cross[0]
                        t_idx = int((prev_cross[0]-time_cell_bounds[0])//agg_sec)
                        delta_x = crosses[i][1] - prev_cross[1]
                        x_idx = int((prev_cross[1] - space_cell_bounds[0])//agg_dist)

                        # update matrix
                        times_sum_mat[x_idx, t_idx] = times_sum_mat[x_idx, t_idx] + delta_t
                        spaces_sum_mat[x_idx, t_idx] = spaces_sum_mat[x_idx, t_idx] + delta_x

                        # keep track of last cross
                        prev_cross = crosses[i]

        # compute flow, density based on edie's definitions
        spaces_sum_mat = np.matrix(spaces_sum_mat)
        times_sum_mat = np.matrix(times_sum_mat)
        q_mat = np.true_divide(spaces_sum_mat, cell_area)   # veh/s
        k_mat = np.true_divide(times_sum_mat, cell_area)    # veh/m

        with np.errstate(divide='ignore', invalid='ignore'):
            v_mat = np.true_divide(spaces_sum_mat, times_sum_mat)  # m/s

        # this is really slow but a quick hack for the time being
        for t_idx in range(1, len(time_cell_bounds)-1):
            for x_idx in range(1, len(space_cell_bounds)-1):
                if np.isnan(v_mat[x_idx, t_idx]):
                    v_mat[x_idx, t_idx] = v_mat[x_idx, t_idx - 1]

                # save the true states on the appropriate files
        true_density_file = (self.__true_state_data_folder[replication] +
                             'truestate_{0}s{1}m_density.txt'.format(agg_sec, agg_dist))
        true_speed_file = (self.__true_state_data_folder[replication] +
                           'truestate_{0}s{1}m_speed.txt'.format(agg_sec, agg_dist))
        true_queue_file = (self.__true_state_data_folder[replication] +
                           'truestate_{0}s{1}m_queue.txt'.format(agg_sec, agg_dist))
        true_traveltime_file = (self.__true_state_data_folder[replication] +
                                'truestate_{0}s{1}m_true_traveltime.txt'.format(agg_sec, agg_dist))
        measured_traveltime_file = (self.__true_state_data_folder[replication] +
                                    'truestate_{0}s{1}m_measured_traveltime.txt'.format(agg_sec, agg_dist))

        # matrix for queue
        dim_space, dim_time = v_mat.shape
        queue_mat = np.zeros((1, dim_time))

        with np.errstate(invalid='ignore'):
            for step in range(0, dim_time):
                index = v_mat[:, step] <= queue_threshold

                if sum(index) <= 2:  # use 2 to suppress false alarms
                    # if less or equal then 2 cells are in congestion, it may be caused by noise.
                    queue_mat[0, step] = 0
                else:
                    queue_mat[0, step] = \
                        agg_dist*(dim_space - np.argmax(index))


        # generate parameters for perfect bluetooth sensors
        sensor_paras = OrderedDict()
        sensor_paras['id'] = 'BTperfect'
        sensor_paras['type'] = 'bluetooth'
        i = 0
        section_up = self.__main_grid_order[i]
        while self.__space_start > self.__sections[section_up]['offset'] + self.__sections[section_up]['length']:
            i += 1
            section_up = self.__main_grid_order[i]
        sensor_paras['section_up'] = section_up
        sensor_paras['distance_up'] = self.__space_start - self.__sections[section_up]['offset']
        j = -1
        section_down = self.__main_grid_order[j]
        while self.__space_end < self.__sections[section_down]['offset']:
            j -= 1
            section_down = self.__main_grid_order[j]
        sensor_paras['section_down'] = section_down
        sensor_paras['distance_down'] = self.__space_end - self.__sections[section_down]['offset']

        sensor_paras['default'] = {}
        sensor_paras['default']['alpha'] = 0
        sensor_paras['default']['theta'] = 0
        sensor_paras['default']['s'] = 1
        sensor_paras['default']['L'] = 3.5
        sensor_paras['default']['offset_dist'] = 0

        sensor_paras['default']['aggregation_sec'] = agg_sec
        sensor_paras['default']['p_penetration'] = [100,100]
        sensor_paras['default']['t_noise_sigma'] = [0.0,0.0]
        sensor_paras['default']['t_noise_mu'] = [0.0,0.0]

        # generate true and measured travel times
        true_tt_mat = np.full(queue_mat.shape, np.nan)
        true_pair_readings = self.__build_travel_time_sensor(sensor_paras,
                                                             save2file=False, agg_sec=agg_sec, true_state=True)
        true_tt_mat_with_time = np.array(true_pair_readings[replication]['tt'])
        true_tt_beg_idx = true_pair_readings[replication]['tt_beg_idx']
        true_tt_mat[0, true_tt_beg_idx:len(true_tt_mat_with_time)+true_tt_beg_idx] = np.transpose(true_tt_mat_with_time)

        measured_tt_mat = np.full(queue_mat.shape, np.nan)
        measured_pair_readings = self.__build_travel_time_sensor(sensor_paras,
                                                                 save2file=False, agg_sec=agg_sec, true_state=False)
        measured_tt_mat_with_time = np.array(measured_pair_readings[replication]['tt'])
        measured_tt_beg_idx = measured_pair_readings[replication]['tt_beg_idx']
        measured_tt_mat[0, measured_tt_beg_idx:len(measured_tt_mat_with_time) + measured_tt_beg_idx] = \
            np.transpose(measured_tt_mat_with_time)

        np.savetxt(true_speed_file, v_mat.T, delimiter=',')
        np.savetxt(true_density_file, k_mat.T, delimiter=',')
        np.savetxt(true_queue_file, queue_mat.T, delimiter=',')
        np.savetxt(true_traveltime_file, true_tt_mat.T, delimiter=',')
        np.savetxt(measured_traveltime_file, measured_tt_mat.T, delimiter=',')

    @staticmethod
    def __interpolate_space_crosses(time_spl, space_cell_bounds):

        """
        This function interpolates on space given a time spline
        :param time_spl: time spline that takes space
        :param space_cell_bounds: space cell bounds
        :return vector of tuples time_space_cross
        """

        time_space_cross = []
        for space_bound in space_cell_bounds:
            time_space_cross.append([time_spl(space_bound), space_bound])
        return time_space_cross

    @staticmethod
    def __interpolate_time_crosses(space_spl, time_cell_bounds):

        """
        This function interpolates on time given a space spline
        :param space_spl: space spline that takes time
        :param time_cell_bounds: time cell bounds
        :return vector of tuples time_space_cross
        """

        time_space_cross = []
        for time_bound in time_cell_bounds:
            time_space_cross.append([time_bound, space_spl(time_bound)])
        return time_space_cross

    @staticmethod
    def __sensordict2string( sensor_id, sensor_att):
        """
        This function converts the sensor key-value store to a string in the following format.
        Entries are separated by ; The entry name and value are separated by :, first entry is id:
        :param sensor_id: the sensor id
        :param sensor_att: key-value store of the sensors.
                general_keys: id, type, section, distance...
                default: Dict, first copied from default parameters and then overwritten by custom values. Use this dict for parameters
                custom: Dict, the set of parameters that have been written to default.
        :return: a string in the above defined format
        """

        line = []
        # append the first id entry
        entry = ':'.join( ['id', sensor_id ] )
        line.append(entry)

        # append other attribute entries
        for key in sensor_att.keys():

            if key != 'id':
                # Do Not repeat id entry

                if key == 'custom':
                    # skip the custom entries since they have been written to default sets.
                    continue
                elif key == 'default':
                    # add all parameters in default dict
                    for para in sensor_att['default']:
                        entry = ':'.join( [ para, str(sensor_att['default'][para]) ] )
                        line.append(entry)
                else:
                    # other keys, such as type, distance...
                    entry = ':'.join( [ key, str(sensor_att[key]) ] )
                    line.append(entry)

        # join all entries and return
        return ';'.join( line )

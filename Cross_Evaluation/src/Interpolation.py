from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np

from Estimator import Estimator

__author__ = 'Yanning Li'
"""
This class implements the spatial smoothing and filtering algorithm used in practice.
"""


class Interpolation(Estimator):
    def __init__(self, time_grid=None, space_grid=None,
                 loc_onramp=None, loc_offramp=None, sensors=None,
                 queue_threshold=17.88, missing_data='fill'):
        """
        Constructor of interpolation methods
        :param time_grid: list, time grid points
        :param space_grid: list, time grid points including first (maybe nonzero) grid.
        :param loc_onramp: ordered list of floats, each float is the absolute location on the freeway
        :param loc_offramp: ordered list of floats, each float is the absolute location on the freeway
        :param sensors: dict,
                sensors[s_id]: with keys: type, loc, cell, flow_std, speed_std, aggregation_sec
        :param queue_threshold: meters.
        :param missing_data: 'fill': fill in the data from adjacent sensors
                             'blank': fill with nan
                             'no_update': fill with previous available data
        :return:
        """

        Estimator.__init__(self, time_grid, space_grid, loc_onramp, loc_offramp,
                           sensors, queue_threshold)

        self.num_cells = len(space_grid) - 1
        self.sim_steps = len(time_grid) - 1
        self.len_cell = space_grid[1] - space_grid[0]
        self.step_dur = time_grid[1] - time_grid[0]
        self.queue_threshold = queue_threshold
        self.missing_data = missing_data

        self.speed_data = None

    def constant_interpolation(self):
        """
        This function simply assumes the speed measurement of one detector is the constant
        speed ranging halfway from its upstream adjacent detector to halfway to its downstream
        adjacent detector.
        :return:
        """

        t_grid = np.array(self.time_grid[1:])

        t_start = 0
        for meas in range(0, len(self.speed_data['time'])):

            # get effective measurement based on missing_data option
            sensor_to_section = self.__get_sensor_coverage_midpoint(self.speed_data['sensors'],
                                                                    self.speed_data['data'][:, meas])

            # get the interval time
            t_now = self.speed_data['time'][meas]
            time_index = (t_start < t_grid) & (t_grid <= t_now)

            # map the effective measurement to sections
            for s_id in sensor_to_section.keys():
                idx = self.speed_data['sensors'].index(s_id)

                # start and end cell of the constant measurement
                sec_start = sensor_to_section[s_id][0]
                sec_end = sensor_to_section[s_id][1]

                self.est_speed[sec_start:sec_end + 1, time_index] = self.__get_meas(idx, meas)

            # update t_start of the interval
            t_start = t_now

        # TODO, convert the speed to density

        # compute the snap-shot-speed-based travel time.
        for step in range(0, self.sim_steps):
            self.est_traveltime[step] = sum(self.len_cell / self.est_speed[:, step])

        # compute the end of the queue
        # simply consider the first cell that has a lower speed than the threshold as the end of queue
        for step in range(0, self.sim_steps):
            # index of the cells that has lower speed
            index = (self.est_speed[:, step] <= self.queue_threshold)

            if sum(index) <= 2:
                # all uncongested
                self.est_queue[step] = 0
            else:
                self.est_queue[step] = \
                    self.len_cell * (self.num_cells - np.argmax(index))

    def linear_interpolation(self):
        """
        This function linearly interpolates the speed in space between two detectors.
        The first cell to the first cell with a sensor will take the value of the measurement of the sensor
        :return:
        """
        # for each time step
        t_grid = np.array(self.time_grid[1:])
        t_start = 0
        for meas in range(0, len(self.speed_data['time'])):

            # get the effective measurement
            speed_sensor_cells, speed_sensor_ids = self.__get_effective_sensor_measures(self.speed_data['sensors'],
                                                                                        self.speed_data['data'][:,
                                                                                        meas])

            # get the time index span that should be mapped to: [t_start, t_now]
            t_now = self.speed_data['time'][meas]
            time_index = (t_start < t_grid) & (t_grid <= t_now)
            dt_agg = t_now - t_start

            # map each sensor into corresponding section linearly
            sec_start_cell = 0
            for i in range(0, len(speed_sensor_cells)):

                s_id = speed_sensor_ids[i]

                if i == 0:
                    # the first section takes the constant measurement from the sensor
                    sec_end_cell = speed_sensor_cells[i] - 1
                    v_sensor_index = self.speed_data['sensors'].index(speed_sensor_ids[i])

                    self.est_speed[sec_start_cell:sec_end_cell + 1, time_index] = self.__get_meas(v_sensor_index, meas)
                    sec_start_cell = sec_end_cell + 1

                else:
                    sec_end_cell = speed_sensor_cells[i] - 1
                    # linear interpolation between each pair of sensors
                    v_i_index = self.speed_data['sensors'].index(speed_sensor_ids[i])
                    v_i_1_index = self.speed_data['sensors'].index(speed_sensor_ids[i - 1])
                    v_i = self.__get_meas(v_i_index, meas)
                    v_i_1 = self.__get_meas(v_i_1_index, meas)

                    # dv may be nan in the case if missing_data == nan, or no_update if one sensor never worked.

                    # print('v_i: {0}; v_i-1:{1}'.format(v_i, v_i_1))
                    if np.isinf(v_i) and np.isinf(v_i_1):
                        dv = np.inf
                    else:
                        dv = float(v_i - v_i_1) \
                             / (speed_sensor_cells[i] - speed_sensor_cells[i - 1])

                    # if the flow is inf, which means no car detected
                    if np.isinf(dv):
                        self.est_speed[sec_start_cell: sec_end_cell + 1, time_index] = np.inf
                    elif not np.isnan(dv):

                        if dv != 0:
                            self.est_speed[sec_start_cell: sec_end_cell + 1, time_index] = \
                                np.matrix(np.arange(v_i_1 + dv / 2, v_i, dv)).T * \
                                np.matrix(np.ones((1, int(round(dt_agg / self.step_dur)))))
                        else:

                            self.est_speed[sec_start_cell: sec_end_cell + 1, time_index] = \
                                np.matrix(np.ones((sec_end_cell + 1 - sec_start_cell, 1), float)) * \
                                np.matrix(np.ones((1, int(round(dt_agg / self.step_dur))))) * v_i
                    else:
                        self.est_speed[sec_start_cell: sec_end_cell + 1, time_index] = np.nan

                    sec_start_cell = sec_end_cell + 1

            # last section use a constant speed from the last sensor
            v_sensor_index = self.speed_data['sensors'].index(speed_sensor_ids[-1])
            sec_end_cell = self.num_cells - 1
            self.est_speed[sec_start_cell: sec_end_cell + 1, time_index] = self.__get_meas(v_sensor_index, meas)

            # update the start time of the interval
            t_start = t_now


            # compute the snap-shot-speed-based travel time.
            # for step in range(0, self.sim_steps):
            # self.est_traveltime[step] = sum( self.len_cell/self.est_speed[:,step] )

            # compute the end of the queue
            # simply consider the first cell that has a lower speed than the threshold as the end of queue
            # for step in range(0, self.sim_steps):
            #     # index of the cells that has lower speed
            #     index = (self.est_speed[:,step] <= self.queue_threshold)
            #
            #     if sum(index) <= 2:
            #         # all un-congested
            #         self.est_queue[step] = 0
            #     else:
            #         self.est_queue[step] = \
            #             self.len_cell*( self.num_cells - np.argmax( index ) )

    def __get_sensor_coverage_midpoint(self, data_sensor_ids, data):
        """
        This funciton maps a point measurement from a sensor to a neighboring section.
        :param data_sensor_ids: the list of sensor ids in the order of the data
        :param data: np.matrix, len(data_sensor_ids) x 1, the raw measurement, may have np.inf, -1, or np.nan
        :return: sensor_to_section, dict
                keys: sensor_id: [start_cell, end_cell]
        """

        speed_sensors = []

        if self.missing_data == 'fill':
            # only record the sensors with effective data.
            for i in range(0, len(data_sensor_ids)):
                if np.isnan(data[i, 0]) or np.isinf(data[i, 0]) or data[i, 0] < 0:
                    # data is missing
                    pass
                else:
                    speed_sensors.append((self.sensors[data_sensor_ids[i]]['cell'], data_sensor_ids[i]))

        if self.missing_data == 'blank' or self.missing_data == 'no_update' or len(speed_sensors) == 0:
            # simply map the measurement to neighbouring section as if no data is missing
            # or if missing_data == fill but no data from any sensor.

            for s_id in data_sensor_ids:
                speed_sensors.append((self.sensors[s_id]['cell'], s_id))

        if self.missing_data != 'fill' and self.missing_data != 'blank' and self.missing_data != 'no_update':
            raise Exception('Error: unrecognized missing_data option')

        speed_sensors = sorted(speed_sensors)
        speed_sensor_cells = [t[0] for t in speed_sensors]
        speed_sensor_ids = [t[1] for t in speed_sensors]

        # Use a key-value store to save how the speed measurement maps to a section
        sensor_to_section = OrderedDict()

        sec_start_cell = 0
        for i in range(0, len(speed_sensor_cells) - 1):
            sec_end_cell = (speed_sensor_cells[i] + speed_sensor_cells[i + 1] - 1) / 2
            sensor_to_section[speed_sensor_ids[i]] = [sec_start_cell, sec_end_cell]

            # the start of the next section
            sec_start_cell = sec_end_cell + 1

        sensor_to_section[speed_sensor_ids[-1]] = [sec_start_cell, self.num_cells - 1]

        return sensor_to_section

    def __get_effective_sensor_measures(self, data_sensor_ids, data):
        """
        This function checks the data for nan values and returns the effective non-nan values.
        :param data_sensor_ids: the list of sensor ids in the order of the data
        :param data: np.matrix, len(data_sensor_ids) x 1, the raw measurement, may have np.inf, -1, or np.nan
        :return: sensor_to_section, dict
                keys: sensor_id: [start_cell, end_cell]
        """

        speed_sensors = []

        if self.missing_data == 'fill':
            # only record the sensors with effective data.
            for i in range(0, len(data_sensor_ids)):
                if np.isnan(data[i, 0]) or np.isinf(data[i, 0]) or data[i, 0] < 0:
                    # data is missing
                    pass
                else:
                    speed_sensors.append((self.sensors[data_sensor_ids[i]]['cell'], data_sensor_ids[i]))

        if self.missing_data == 'blank' or self.missing_data == 'no_update' or len(speed_sensors) == 0:
            # simply map the measurement to neighbouring section as if no data is missing
            # or if missing_data == fill but no data from any sensor, then fill all with nan

            for s_id in data_sensor_ids:
                speed_sensors.append((self.sensors[s_id]['cell'], s_id))

        if self.missing_data != 'fill' and self.missing_data != 'blank' and self.missing_data != 'no_update':
            raise Exception('Error: unrecognized missing_data option')

        speed_sensors = sorted(speed_sensors)
        speed_sensor_cells = [t[0] for t in speed_sensors]
        speed_sensor_ids = [t[1] for t in speed_sensors]

        return speed_sensor_cells, speed_sensor_ids

    def __get_meas(self, sensor_idx, meas_idx):
        """
        This function returns the sensor measurement based on the missing_data option.
        It returns the latest non-nan value for no_update missing data option, and the current value for other options
        :param sensor_idx: the index in the speed_data for the sensor
        :param meas_idx: the current measurement index
        :return:
        """
        if self.missing_data == 'no_update':
            # find the latest non-nan data
            temp_meas = meas_idx
            while temp_meas != 0 and np.isnan(self.speed_data['data'][sensor_idx, temp_meas]):
                temp_meas -= 1
            return self.speed_data['data'][sensor_idx, temp_meas]
        else:
            return self.speed_data['data'][sensor_idx, meas_idx]

    def plot_estimated_speed(self, V_MAX=40):
        """
        This function plots the estimated speed profile in the entire time horizon
        :param V_MAX: m/s, this is the maximal speed that to normalize the speed map
        :return:
        """
        # The speed estimate is in the speed_estimates matrix.
        speed_est = np.flipud(self.est_speed[0:self.num_cells, :])

        fig = plt.figure(figsize=(10, 10), dpi=100)
        im = plt.imshow(speed_est, cmap=plt.get_cmap('jet'),
                        interpolation='nearest',
                        vmin=0, vmax=V_MAX)

        plt.title('Estimated speed profile')
        cax = fig.add_axes([0.02, 0.25, 0.01, 0.4])
        fig.colorbar(im, cax=cax, orientation='vertical')
        plt.draw()

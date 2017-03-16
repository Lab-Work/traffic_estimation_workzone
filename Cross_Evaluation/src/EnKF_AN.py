import numbers
import sys
from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np

from EnKF import EnKF
from Estimator import Estimator

__author__ = 'Yanning Li'
"""
This is the EnKF using additive noise to the estimate density state.
The states are the densities and the qin and qout.
The onramp and offramp are included in the larger additive noise in the corresponding cell.
The measurements are the velocity and flows at cell boundaries.
"""


class EnKF_AN(Estimator, EnKF):
    """
    This is the Ensemble filter for traffic estimation using the CTM scheme.
    This EnKF is based on the quadratic linear fundamental diagram
    This class inherit the EnKF and the Estimator class, it inherits the architecture of Estimator and overwrite the state update
    and observation class.

    Here is an example of using this class

    filter = EnKF_workzone(.....)
    filter.set_meas_data(speed_data, flow_data, BT_data)
    filter.run_filter()

    # or you can get the estimate by
    filter.get_density_est()
    """

    def __init__(self, time_grid=None, space_grid=None,
                 sensors=None,
                 loc_onramp=None, loc_offramp=None,
                 vm_cells=None, beta_cells=None, rhoc_cells=None, wc_cells=None,
                 num_ensembles=0,
                 std_model_noise=None, queue_threshold=17.88,
                 init_rho=0, init_qin=0.5, init_qout=0.0):
        """
        The constructor of the EnKF estimator for the speed, travel time, and queue in a work zone.
        :param time_grid: [], second, the finest time grid including starting and ending time.
        :param space_grid: [], meter, the finest space grid including entrance and exit positions.
        :param sensors:
        :param loc_onramp: [], meters, locations of onramps
        :param loc_offramp: [], meters, locations of offramps
        :param vm_cells: [], m/s, freeflow speed in each space cell
        :param wc_cells: [], m/s, >0, congestion wave propagation speed in each cell
        :param num_ensembles: int, number of ensembles
        :param std_model_noise: Dict, 'vf', 'w', 'qmax', 'qin' 'qout' 'onramp','offramp'
                    Q['cell'] = std, the additive noise for normal cell
                    Q['oncell'] = std, the larger additive noise for onramp and offramp cell
                    Q['offcell'] = std, the larger additive noise for onramp and offramp cell
                    Q['qin'] = float, std
                    Q['qout'] = float, std
        :param queue_threshold: m/s, the threshold to determine congestion. If speed in a cell is below, then congested.
        :return:
        """

        self.__debug = False
        self.__debug_entrance_sensor = 'IDEALLoc100m'
        self.__debug_exit_sensor = 'IDEALLoc8300m'

        # initialize the superclass Estimator
        Estimator.__init__(self, time_grid, space_grid,
                           loc_onramp, loc_offramp,
                           sensors,
                           queue_threshold)

        # build the index for the system state
        self.x_index, dim_state = self.__build_state_index()

        # initialize the super class
        EnKF.__init__(self, dim_state, num_ensembles)

        # y_index, and dim_obs, which will be dynamically updated upon arrival of each new data
        self.y_index = None
        self.dim_obs = None

        # keep track of the flow between cells for each ensemble which will be used to construct the observation
        self.__f_flow = {}
        self.__f_flow['time'] = np.array(self.time_grid[1:])
        self.__f_flow['data'] = OrderedDict()
        for i in range(0, self.num_ensembles):
            self.__f_flow['data'][i] = []

        # keep track of the speed between cells for each ensemble which will be used to construct the observation
        self.__f_speed = {}
        self.__f_speed['time'] = np.array(self.time_grid[1:])
        self.__f_speed['data'] = OrderedDict()
        for i in range(0, self.num_ensembles):
            self.__f_speed['data'][i] = []

        # save all the estimated states here
        self.est_state_all = np.matrix(np.zeros((self.dim_state, self.num_steps), float))

        # =================================================
        # Add additive noise to state.
        self.Q = OrderedDict()
        # initialize with all cell var
        self.Q = np.diag(np.ones(dim_state) * (std_model_noise['cell'] ** 2))

        # print('onramps:{0}; offramps:{1}'.format(self.cell_onramp, self.cell_offramp))
        # add onramp larger noise
        if self.cell_onramp is not None:
            for on_cell in self.cell_onramp:
                if 0 <= on_cell <= self.num_cells:
                    idx = self.x_index['density'][on_cell]
                    self.Q[idx, idx] = std_model_noise['oncell'] ** 2
        # add offramp larger noise
        if self.cell_offramp is not None:
            for off_cell in self.cell_offramp:
                if 0 <= off_cell <= self.num_cells:
                    idx = self.x_index['density'][off_cell]
                    self.Q[idx, idx] = std_model_noise['offcell'] ** 2
        # add qin variance
        idx = self.x_index['qin']
        self.Q[idx, idx] = std_model_noise['qin'] ** 2
        # add qout variance
        idx = self.x_index['qout']
        self.Q[idx, idx] = std_model_noise['qout'] ** 2

        # self.Q = std_model_noise
        # if np.size( self.Q['vm'] ) == 1:
        #     # if it was a single value, then it was specified as std, not var (which = std^2)
        #     self.Q['vm'] = np.diag( np.ones( self.num_cells )*(self.Q['vm']**2) )
        # if np.size( self.Q['beta'] ) == 1:
        #     self.Q['beta'] = np.diag( np.ones( self.num_cells )*(self.Q['beta']**2) )
        # if np.size( self.Q['rhoc'] ) == 1:
        #     self.Q['rhoc'] = np.diag( np.ones( self.num_cells )*(self.Q['rhoc']**2) )
        # if np.size( self.Q['wc'] ) == 1:
        #     self.Q['wc'] = np.diag( np.ones( self.num_cells )*(self.Q['wc']**2) )
        #
        # if self.loc_onramp is not None and np.size(self.Q['onramp']) == 1:
        #     self.Q['onramp'] = np.diag( np.ones(len(loc_onramp))*(self.Q['onramp']**2) )
        # if self.loc_offramp is not None and np.size(self.Q['offramp']) == 1:
        #     self.Q['offramp'] = np.diag( np.ones(len(loc_offramp))*(self.Q['offramp']**2) )


        # =================================================
        # save the fundamental diagram for each cell
        # vm parameter
        if isinstance(vm_cells, numbers.Number):
            self.vm_cells = np.ones((self.num_cells, 1)) * float(vm_cells)
        else:
            self.vm_cells = np.array(vm_cells).astype(float)
            self.vm_cells = self.vm_cells.reshape((self.num_cells, 1))

        # beta parameter
        if isinstance(beta_cells, numbers.Number):
            self.beta_cells = np.ones((self.num_cells, 1)) * float(beta_cells)
        else:
            self.beta_cells = np.array(beta_cells).astype(float)
            self.beta_cells = self.beta_cells.reshape((self.num_cells, 1))

        # rhoc parameter
        if isinstance(rhoc_cells, numbers.Number):
            self.rhoc_cells = np.ones((self.num_cells, 1)) * float(rhoc_cells)
        else:
            self.rhoc_cells = np.array(rhoc_cells).astype(float)
            self.rhoc_cells = self.rhoc_cells.reshape((self.num_cells, 1))

        # wc parameter
        if isinstance(wc_cells, numbers.Number):
            self.wc_cells = np.ones((self.num_cells, 1)) * float(wc_cells)
        else:
            self.wc_cells = np.array(wc_cells).astype(float)
            self.wc_cells = self.wc_cells.reshape((self.num_cells, 1))

        # other use ful parameters
        self.qmax_cells = self.vm_cells * self.rhoc_cells - \
                          self.vm_cells * (self.rhoc_cells ** 2) / self.beta_cells

        self.rhomax_cells = - self.qmax_cells / self.wc_cells + self.rhoc_cells

        # =======================================================================
        self.init_rho = init_rho
        self.init_qin = init_qin
        self.init_qout = init_qout

        # =======================================================================
        # FOR DEBUGGING
        # recored the forecast and analysis value for qin and qout
        if self.__debug:
            self.qin_f = []
            self.qin_a = []
            self.qin_obs = []
            self.qout_f = []
            self.qout_a = []
            self.qout_obs = []

    def __build_state_index(self):
        """
        In this nonlinear CTM system, the state variable is defined as
            x = [rho_1, ..., rho_I, q_in, q_out]
            where r_i, and f_i are the on off ramp flows. This function build a index for easy access to the state
            variables. Note, there is less than I r_i or f_i since not all cells have an on off ramp in it.
        :return:
            x_index, dict. keys: 'density', 'qin', 'qout', 'onramp', 'offramp'
                x_index['density', or 'onramp', or 'offramp'] = dict: [float_keys] float_key is the cell id, gives
                    the index of the corresponding variable. E.g. onramp r_i of cell i index: x_index['onramp'][i]
                x_index['qin','qout'] = int, the absolute index for qin and qout
            dim_state: int, the dimension of the state
        """

        # the index for the system state
        # [rho_i, q_in, q_out, r_i, f_i]
        x_index = {}

        # add the density index
        x_index['density'] = OrderedDict()
        for i in range(0, self.num_cells):
            x_index['density'][i] = i
        dim_state = self.num_cells

        # add the upstream boundary flow
        x_index['qin'] = dim_state
        x_index['qout'] = dim_state + 1
        dim_state += 2

        # add on ramp variables
        # x_index['onramp'] = OrderedDict()
        # if self.cell_onramp is not None:
        #     # if onramp exist in the network, otherwise skip
        #     for cell_id in self.cell_onramp:
        #         # add the absolute index into the state index dictionary
        #         # r_i index  = self.x_index{'onramp'][cell_i]
        #         x_index['onramp'][cell_id] = dim_state
        #         dim_state += 1
        #
        #
        # # add off ramp state variables
        # x_index['offramp'] = OrderedDict()
        # if self.cell_offramp  is not None:
        #     for cell_id in self.cell_offramp:
        #         # add the absolute index
        #         x_index['offramp'][cell_id] = dim_state
        #         dim_state += 1

        return x_index, dim_state

    def __build_obs_index(self, flow_data, speed_data, traveltime_data):
        """
        This function builds an observation index based on the sensors at this time instance.
        It also builds a Dict which maps the location to the observation index, which will be used to extract data.
            y = [ q_i, r_i, f_i, v_i, tau_j],
            where q_i, and v_i are sensors installed at the upstream boundary of cell i, r_i and f_i are the onramp and
            off ramps. tau_j is the travel time for bluetooth sensor pair j

        :return:
            y_index: dict, keys: 'flow', 'speed', 'qin', 'qout', 'onramp', 'offramp', 'travel_time'
                y_index['flow',or 'speed', or 'onramp', or 'offramp'] gives dict with keys being cell ID, giving index
                    in the observation array
                y_index['travel_time']['3_89'] give the travel time measurement index in the observation array. key is
                    a string composed by '{0}_{1}'.format(tuple(0), tuple(1))
            dim_obs: int, the dimension of the observation array
        """
        y_index = OrderedDict()
        dim_obs = 0
        y_obs = []
        cov_obs = []

        # assume the measurement noise of each sensor is independent with std
        std = []

        # add q_i, r_i, or f_i
        if flow_data is not None:
            for s_id in flow_data['sensors']:
                if 'flow' not in y_index.keys():
                    y_index['flow'] = OrderedDict()
                # map the sensor id to the index
                y_index['flow'][s_id] = dim_obs
                dim_obs += 1

                # add the flow std of this sensor to construct covariance
                std.append(float(self.sensors[s_id]['flow_std']))

            y_obs += flow_data['data']

        # add v_i
        if speed_data is not None:
            for s_id in speed_data['sensors']:
                if 'speed' not in y_index.keys():
                    y_index['speed'] = OrderedDict()
                y_index['speed'][s_id] = dim_obs
                dim_obs += 1

                # add the flow std of this sensor to construct covariance
                std.append(float(self.sensors[s_id]['speed_std']))

            y_obs += speed_data['data']

        # add \tau_i
        if traveltime_data is not None:
            for s_id in traveltime_data['sensors']:
                if 'travel_time' not in y_index.keys():
                    y_index['travel_time'] = OrderedDict()
                y_index['travel_time'][s_id] = dim_obs
                dim_obs += 1

                # add the flow std of this sensor to construct covariance
                std.append(float(self.sensors[s_id]['traveltime_std']))

            y_obs += traveltime_data['data']

        # Assuming independent variables. Then each entry is the variance = std^2
        if len(std) != 0:
            cov_obs = np.matrix(np.diag(std)) * np.matrix(np.diag(std))

        return y_index, dim_obs, np.matrix(y_obs).reshape((dim_obs, 1)), cov_obs

    def __get_eff_meas(self, cur_time):
        """
        This function returns the effective measurement (those that are not -1 or NaN) at current time from self.flow_data,
        self.speed_data, and self.traveltime_data
        :param time: the timestamp of the measurement, must be on time grid
        :return: eff_flow_data, eff_speed_data, eff_traveltime_data as Dict, None if no data available
        """
        eff_flow_data = None
        eff_speed_data = None
        eff_traveltime_data = None

        # get the effective speed sensor data
        if self.speed_data is not None:
            if cur_time in self.speed_data['time']:
                t_index = self.speed_data['time'].index(cur_time)

                for s_index in range(0, self.speed_data['data'].shape[0]):

                    if ~np.isnan(self.speed_data['data'][s_index, t_index]) and \
                            ~np.isinf(self.speed_data['data'][s_index, t_index]) and \
                                    self.speed_data['data'][s_index, t_index] >= 0:

                        if eff_speed_data is None:
                            eff_speed_data = {'sensors': [], 'data': []}

                        eff_speed_data['sensors'].append(self.speed_data['sensors'][s_index])
                        eff_speed_data['data'].append(self.speed_data['data'][s_index, t_index])

        # get the effective flow sensor data
        if self.flow_data is not None:
            if cur_time in self.flow_data['time']:
                t_index = self.flow_data['time'].index(cur_time)
                # check all sensors in this time slice, and throw away those sensors with no reading
                for s_index in range(0, len(self.flow_data['sensors'])):

                    # check if effective speed data exists, if not skip
                    if eff_speed_data is None or self.flow_data['sensors'][s_index] not in eff_speed_data['sensors']:
                        continue

                    if ~np.isnan(self.flow_data['data'][s_index, t_index]) and \
                            ~np.isinf(self.flow_data['data'][s_index, t_index]) and \
                                    self.flow_data['data'][s_index, t_index] >= 0:

                        if eff_flow_data is None:
                            eff_flow_data = {'sensors': [], 'data': []}

                        eff_flow_data['sensors'].append(self.flow_data['sensors'][s_index])

                        # Radar sensors have severe occlusion issue.
                        # In data, it was found to be x2 lower. Here compensate.
                        if 'RADAR' in self.flow_data['sensors'][s_index]:
                            # find out the state
                            sensor_id = self.flow_data['sensors'][s_index]

                            for s_idx, s_id in enumerate(eff_speed_data['sensors']):

                                if s_id == sensor_id:
                                    s_speed = eff_speed_data['data'][s_idx]
                                    if s_speed <= 17.88:
                                        # in free flow, occlusion is not that bad
                                        coe = 2.0
                                    else:
                                        coe = 1.8

                            eff_flow_data['data'].append(self.flow_data['data'][s_index, t_index] * coe)
                            # print('multiplied radar {0} by {1} at speed {2}'.format(self.flow_data['sensors'][s_index],
                            #                                                         coe, s_speed))

                        # ICONE has worse occlusion and sleep cycle issue. It was found in the simulated data
                        # the icone flow is around 3~4 times lower.
                        elif 'ICONE' in self.flow_data['sensors'][s_index]:
                            # find out the state
                            sensor_id = self.flow_data['sensors'][s_index]
                            for s_idx, s_id in enumerate(eff_speed_data['sensors']):
                                if s_id == sensor_id:
                                    s_speed = eff_speed_data['data'][s_idx]
                                    if s_speed <= 17.88:
                                        # in free flow, occlusion is not that bad
                                        coe = 4.0
                                    else:
                                        coe = 3.0

                            eff_flow_data['data'].append(self.flow_data['data'][s_index, t_index] * coe)
                            # print('multiplied radar {0} by {1} at speed {2}'.format(self.flow_data['sensors'][s_index],
                            #                                                         coe, s_speed))
                        else:
                            eff_flow_data['data'].append(self.flow_data['data'][s_index, t_index])

        # get the effective travel time sensor data
        if self.traveltime_data is not None:
            if cur_time in self.traveltime_data['time']:
                t_index = self.traveltime_data['time'].index(cur_time)

                for s_index in range(0, self.traveltime_data['data'].shape[0]):
                    if ~np.isnan(self.traveltime_data['data'][s_index, t_index]) and \
                            ~np.isinf(self.traveltime_data['data'][s_index, t_index]) and \
                                    self.traveltime_data['data'][s_index, t_index] >= 0:

                        if eff_traveltime_data is None:
                            eff_traveltime_data = {'sensors': [], 'data': []}

                        eff_traveltime_data['sensors'].append(self.traveltime_data['sensors'][s_index])
                        eff_traveltime_data['data'].append(self.traveltime_data['data'][s_index, t_index])

        return eff_flow_data, eff_speed_data, eff_traveltime_data

    # @profile
    def run_batch_filter(self):
        """
        Call this function to run the estimation for a batch data, assume all data are synchronized
        enKF.set_initial_ensembles( X_a )

        for step in range(0, self.num_steps)

            obs = self.all_obs[:, step]
            est_state = enKF.run_EnKF(obs)

        :return:
        """
        if self.speed_data is None and self.flow_data is None and self.traveltime_data is None:
            print(
            'Warning: The measurement data must be set before running the batch filter: use function self.set_meas_data()')

        # =======================================================================
        # the initial ensembles, which should have been set externally
        X_init = np.matrix(np.zeros((self.dim_state, self.num_ensembles)))
        print(
        'Setting initial ensembles: rho {0}; qin {1}; qout {2}'.format(self.init_rho, self.init_qin, self.init_qout))
        for ens in range(0, self.num_ensembles):
            X_init[self.x_index['density'][0]:
            self.x_index['density'][self.num_cells - 1], ens] = self.init_rho
            X_init[self.x_index['qin'], ens] = self.init_qin
            X_init[self.x_index['qout'], ens] = self.init_qout

            # print('setted qin {0}; qout {1}'.format(X_init[self.x_index['qin'], ens], X_init[self.x_index['qout'], ens]  ))
            # add noise to each ensemble
            X_init[:, ens] += np.matrix(np.random.multivariate_normal(
                np.zeros(self.dim_state), self.Q)).reshape((self.dim_state, 1))

        self.set_initial_ensembles(X_init)

        # =======================================================================
        # DEBUG
        # save the qin and qout in the corresponding probe data
        # save the initial state
        if self.__debug:
            self.qin_f.append(np.squeeze(np.array(self.X_f[self.x_index['qin'], :])).tolist())
            self.qin_a.append(np.squeeze(np.array(self.X_a[self.x_index['qin'], :])).tolist())
            self.qin_obs.append(np.nan)

            self.qout_f.append(np.squeeze(np.array(self.X_f[self.x_index['qout'], :])).tolist())
            self.qout_a.append(np.squeeze(np.array(self.X_a[self.x_index['qout'], :])).tolist())
            self.qout_obs.append(np.nan)

        # The enKF runs at the finest time grid
        # for each step, update the system
        for step in range(0, self.num_steps):

            # update status
            sys.stdout.write('\r')
            sys.stdout.write('Status: filtering step {0}/{1}'.format(step, self.num_steps))
            sys.stdout.flush()
            # print('Status: filtering step {0}'.format(step))

            cur_time = (step + 1) * self.dur_steps

            # get the effective measurement
            eff_flow, eff_speed, eff_traveltime = self.__get_eff_meas(cur_time)

            # build the observation index
            self.y_index, self.dim_obs, y_obs, cov_noise = self.__build_obs_index(eff_flow, eff_speed, eff_traveltime)

            # update the estimate for this step
            est_state = self.update_estimate(y_obs, cov_noise, cur_time)

            # =======================================================================
            # DEBUG
            # save the qin and qout in the corresponding probe data
            # save the updated state
            if self.__debug:
                self.qin_f.append(np.squeeze(np.array(self.X_f[self.x_index['qin'], :])).tolist())
                self.qin_a.append(np.squeeze(np.array(self.X_a[self.x_index['qin'], :])).tolist())
                if 'flow' in self.y_index.keys() and self.__debug_entrance_sensor in self.y_index['flow'].keys():
                    self.qin_obs.append(y_obs[self.y_index['flow'][self.__debug_entrance_sensor]])
                    # print('y_index[flow]:{0}'.format(self.y_index['flow'].keys()))
                    # print('y_obs[ y_index[flow][entrance] ]:{0}'.format(
                    #     y_obs[ self.y_index['flow'][self.__debug_entrance_sensor]],
                    #     self.__debug_entrance_sensor))
                else:
                    self.qin_obs.append(np.nan)

                self.qout_f.append(np.squeeze(np.array(self.X_f[self.x_index['qout'], :])).tolist())
                self.qout_a.append(np.squeeze(np.array(self.X_a[self.x_index['qout'], :])).tolist())
                if 'flow' in self.y_index.keys() and self.__debug_exit_sensor in self.y_index['flow'].keys():
                    self.qout_obs.append(y_obs[self.y_index['flow'][self.__debug_exit_sensor]])
                else:
                    self.qout_obs.append(np.nan)
            # =======================================================================
            # save the estimated state
            self.est_state_all[:, step] = est_state

            # decouple and save into self.est_density, self.est_speed, self.est_queue, self.est_traveltime
            self.est_density[:, step] = est_state[0:self.num_cells, 0]

            # the speed is computed using the fundamental diagram
            for cell_id in range(0, self.num_cells):
                # use the static FD at this step
                self.est_speed[cell_id, step] = self.__rho2v(self.vm_cells[cell_id, 0], self.beta_cells[cell_id, 0],
                                                             self.rhoc_cells[cell_id, 0], self.wc_cells[cell_id, 0],
                                                             self.est_density[cell_id, step])

            # REMARK: the queue and travel time a post-processed from the speed field.
            # They are computed in cross_evaluation class for all algorithms
            # the queue length starts from the first cell with speed below queue_threshold to the end of road
            # index = (self.est_speed[:, step] <= self.queue_threshold)
            #
            # # filter out the outliers
            # index_smoothed = deepcopy(index)
            # outlier_max = 3
            # counter = 0
            # for i in range(0, len(index)):
            #
            #     if index[i] == True:
            #         # trigger the coutner
            #         counter += 1
            #     elif index[i] == False and counter != 0:
            #         if counter <= outlier_max:
            #             # found outliers
            #             index_smoothed[ i-counter : i ] = False
            #             # reset counter
            #             counter = 0
            #
            #     # if i != 0 and i != len(index)-1:
            #     #     if sum( index[i-1:i+3] ) >=2:
            #     #         index_smoothed[i] = True
            #     #     else:
            #     #         index_smoothed[i] = False
            #     # elif i == 0:
            #     #     if sum(index[0: 5] ) >= 3:
            #     #         index_smoothed[i] = True
            #     #     else:
            #     #         index_smoothed[i] = False
            #     # elif i == len(index)-1:
            #     #     if sum(index[ i-4 :len(index)]) >= 3:
            #     #         index_smoothed[i] = True
            #     #     else:
            #     #         index_smoothed[i] = False
            #
            # if sum(index_smoothed) <= 3: # use 4 to suppress false alarms
            #     # if less or equal then 2 cells are in congestion, it may be caused by noise.
            #     self.est_queue[step] = 0
            # else:
            #     # if step > 105 and step < 115:
            #         # print(sum(index_smoothed))
            #         # print(index_smoothed)
            #         # print(index)
            #
            #     self.est_queue[step] = \
            #         self.len_cells*( self.num_cells - np.argmax(index_smoothed) )
            # # try:
            # #     first_cong_cell_id = [x[0] for x in enumerate( self.est_speed[:,step] ) if x[1] < self.queue_threshold][0]
            # # except IndexError:
            # #     # no congested cell
            # #     first_cong_cell_id = self.num_cells
            # # # the estimated queue length
            # # self.est_queue[step] = self.len_cells*( self.num_cells - first_cong_cell_id )
            #
            # # the travel time estimate is computed by summing up the travel time in each cell
            # self.est_traveltime[step] = np.sum(self.len_cells/self.est_speed[:,step])


            # =======================================================================
            # DEBUG
            # plot the update
            if self.__debug:
                plot_len = 19
                # qin
                if False:
                    if not np.isnan(self.qin_obs[-1]):
                        fig1 = plt.figure(figsize=(10, 5), dpi=100)
                        ax1 = fig1.add_subplot(111)
                        positions_f = np.arange(0, len(self.qin_f)) - 0.1
                        positions_a = np.arange(0, len(self.qin_a)) + 0.1
                        positions_obs = np.arange(0, len(self.qin_obs))
                        # predicted as red
                        bp = ax1.boxplot(self.qin_f[-plot_len:],
                                         positions=positions_f[-plot_len:], widths=0.15,
                                         patch_artist=False)
                        for box in bp['boxes']:
                            # change outline color
                            box.set(color='#FF4633', linewidth=1)
                            # change fill color
                            # box.set( facecolor = '#FF4633' )
                        # corrected as green
                        bp = ax1.boxplot(self.qin_a[-plot_len:],
                                         positions=positions_a[-plot_len:], widths=0.15, patch_artist=False)
                        for box in bp['boxes']:
                            # change outline color
                            box.set(color='#07891B', linewidth=1)
                            # change fill color
                            # box.set( facecolor = '#07891B' )
                        # measurement as blue
                        ax1.scatter(positions_obs[-plot_len:], self.qin_obs[-plot_len:], color='b', marker='o', s=40,
                                    label='Observation')
                        ax1.set_title('qin')
                        # x_ticks = np.arange(0, len(self.qin_f))
                        # ax1.set_xticks(x_ticks[-plot_len:])
                        plt.show()

                # qout
                if False:
                    if not np.isnan(self.qout_obs[-1]):
                        fig2 = plt.figure(figsize=(10, 5), dpi=100)
                        ax2 = fig2.add_subplot(111)
                        positions_f = np.arange(0, len(self.qout_f)) - 0.1
                        positions_a = np.arange(0, len(self.qout_a)) + 0.1
                        positions_obs = np.arange(0, len(self.qout_obs))
                        # predicted as red
                        bp = ax2.boxplot(self.qout_f[-plot_len:], positions=positions_f[-plot_len:], widths=0.18,
                                         patch_artist=True)
                        for box in bp['boxes']:
                            # change outline color
                            box.set(color='#7570b3', linewidth=1)
                            # change fill color
                            box.set(facecolor='#FF4633')
                        # corrected as green
                        bp = ax2.boxplot(self.qout_a[-plot_len:], positions=positions_a[-plot_len:], widths=0.18,
                                         patch_artist=True)
                        for box in bp['boxes']:
                            # change outline color
                            box.set(color='#7570b3', linewidth=1)
                            # change fill color
                            box.set(facecolor='#07891B')
                        # measurement as blue
                        ax2.scatter(positions_obs[-plot_len:], self.qout_obs[-plot_len:], color='b', marker='o', s=30,
                                    label='Observation')
                        ax2.set_title('qout')
                        # x_ticks = np.arange(0, len(self.qout_f))
                        # ax2.set_xticks(x_ticks[-plot_len:])

                        plt.show()

        # plot the estimated qin and qout
        if self.__debug:
            if True:
                qin = np.squeeze(np.array(self.est_state_all[self.x_index['qin'], :]))
                qin_meas = np.array(self.qin_obs)[1:]
                print(len(qin), len(qin_meas))
                fig1 = plt.figure(figsize=(10, 5), dpi=100)
                ax1 = fig1.add_subplot(111)
                t = np.arange(len(qin))
                ax1.plot(t, qin, 'r-', label='Estimated')
                not_nan = ~np.isnan(qin_meas)
                ax1.plot(t[not_nan], qin_meas[not_nan], 'b', label='Measured')
                ax1.legend()
                ax1.grid(True)
                ax1.set_title('qin')

                plt.draw()

            if True:
                qout = np.squeeze(np.array(self.est_state_all[self.x_index['qout'], :]))
                qout_meas = np.array(self.qout_obs)[1:]
                fig2 = plt.figure(figsize=(10, 5), dpi=100)
                ax2 = fig2.add_subplot(111)
                t = np.arange(len(qout))
                ax2.plot(t, qout, 'r-', label='Estimated')
                not_nan = ~np.isnan(qout_meas)
                ax2.plot(t[not_nan], qout_meas[not_nan], 'b', label='Measured')
                ax2.set_title('qout')
                ax2.legend()
                ax2.grid(True)
                plt.draw()

    # @profile
    def f_model(self, x_a, e_id):
        """
        The state evolution equation. This function propagate the state for one ensemble by one step.
        :param x_a: the previous estimate
        :param e_id: the ensemble id
        :return: x_f for the next time step.
        """

        # ======================================================== #
        # initialize x_f with same dimension as x_a
        x_f = np.matrix(np.zeros(x_a.shape))

        # ======================================================== #
        # First compute the flow and velocity between cells using estimates x_a, given the previous state estimate
        # dim: 1 x self.num_cells+1, flow[0] is the qin, and flow[self.num_cells] is qout
        vm = self.vm_cells
        beta = self.beta_cells
        rhoc = self.rhoc_cells
        wc = self.wc_cells
        qmax = np.multiply(vm, rhoc) - np.multiply(vm, np.multiply(rhoc, rhoc)) / beta
        rhom = rhoc - qmax / wc

        flow = []
        speed = []
        # ---------------------------------
        # the flow in the upstream boundary
        # first_flow = self.__receiving_flow(wc[0,0], rhoc[0,0], rhom[0,0], qmax[0,0],
        #                                               x_a[ self.x_index['density'][0],0])
        # flow.append( np.min(  [ x_a[ self.x_index['qin'], 0], first_flow ] ))
        # # print('inflow:{0}, {1}'.format( x_a[ self.x_index['qin'], 0], first_flow ))
        # if x_a[ self.x_index['qin'], 0] <= first_flow:
        #     # In freeflow condition, then determine speed by the inflow
        #     speed.append( self.__q2v_ff(vm[0,0],beta[0,0], x_a[ self.x_index['qin'], 0]) )
        # else:
        #     # congested flow, then determine by the density in the first cell
        #     speed.append( self.__rho2v(vm[0,0], beta[0,0], rhoc[0,0], wc[0,0],  x_a[ self.x_index['density'][0], 0]) )

        # Updated qin to be the actual inflow to the first cell
        flow.append(x_a[self.x_index['qin'], 0])
        speed.append(self.__rho2v(vm[0, 0], beta[0, 0], rhoc[0, 0], wc[0, 0], x_a[self.x_index['density'][0], 0]))

        for i in range(1, self.num_cells):
            f_sending = self.__sending_flow(vm[i - 1, 0], beta[i - 1, 0], rhoc[i - 1, 0], qmax[i - 1, 0],
                                            x_a[self.x_index['density'][i - 1], 0])
            f_receiving = self.__receiving_flow(wc[i, 0], rhoc[i, 0], rhom[i, 0], qmax[i, 0],
                                                x_a[self.x_index['density'][i], 0])
            if f_sending <= f_receiving:
                flow.append(f_sending)
            else:
                flow.append(f_receiving)

            # flow.append( np.min( [f_sending, f_receiving] ) )

            # get the velocity between cells
            if f_sending < f_receiving:
                # freeflow; if S < R, use v_upstream
                select_i = i - 1
                speed.append(self.__rho2v(vm[select_i, 0], beta[select_i, 0], rhoc[select_i, 0], wc[select_i, 0],
                                          x_a[self.x_index['density'][select_i], 0]))
            elif f_sending > f_receiving:
                # congested; # if S < R, use v_downstream
                select_i = i
                speed.append(self.__rho2v(vm[select_i, 0], beta[select_i, 0], rhoc[select_i, 0], wc[select_i, 0],
                                          x_a[self.x_index['density'][select_i], 0]))
            else:
                # S== R, then use v(rho_c)
                select_i = i - 1
                speed.append(self.__rho2v(vm[select_i, 0], beta[select_i, 0], rhoc[select_i, 0], wc[select_i, 0],
                                          rhoc[select_i, 0]))

        # ---------------------------------
        # the flow in the downstream boundary
        last_i = self.num_cells - 1
        # last_sending = self.__sending_flow(vm[last_i,0],beta[last_i,0],rhoc[last_i,0],qmax[last_i,0],
        #                                                                         x_a[ self.x_index['density'][last_i], 0])
        # flow.append( np.min([ x_a[self.x_index['qout'], 0], last_sending]) )
        # # the velocity in the downstream boundary
        # if last_sending <= x_a[self.x_index['qout'], 0]:
        #     # freeflow, determine by the density in the last cell
        #     speed.append( self.__rho2v(vm[last_i,0], beta[last_i,0], rhoc[last_i,0], wc[last_i,0],
        #                                 x_a[ self.x_index['density'][last_i], 0]) )
        # else:
        #     # congested, determine by the outflow
        #     speed.append( self.__q2v_cf(wc[last_i,0],rhom[last_i,0], x_a[self.x_index['qout'], 0] ) )

        # updated qout to be the actual outflow
        flow.append(x_a[self.x_index['qout'], 0])
        speed.append(self.__rho2v(vm[last_i, 0], beta[last_i, 0], rhoc[last_i, 0], wc[last_i, 0],
                                  x_a[self.x_index['density'][last_i], 0]))

        # append the flow and vel in this step to all forecast_flow and forecast_speed
        self.__f_flow['data'][e_id].append(flow)
        self.__f_speed['data'][e_id].append(speed)

        # ======================================================== #
        # State Propagation
        # Propagate density for step k, considering the on/off ramp
        for i in range(0, self.num_cells):
            x_f[self.x_index['density'][i], 0] = \
                x_a[self.x_index['density'][i], 0] + \
                (flow[i] - flow[i + 1]) * self.dur_steps / self.len_cells

            # # check the on and offramps
            # if i not in self.x_index['onramp'].keys() and i not in self.x_index['offramp'].keys():
            #
            # elif i in self.x_index['onramp'].keys() and i not in self.x_index['offramp'].keys():
            #     # The update here is a bit complicated. The onramp flow here specifies the demand.
            #     # The actual onramp inflow depends on the supply in this cell.
            #     # q_on = min( onramp_demand, Receiving(i)-flow[i] )
            #     q_on = np.min([ x_a[ self.x_index['onramp'][i], 0],
            #                     self.__receiving_flow(wc[i,0],rhoc[i,0],rhom[i,0],qmax[i,0],
            #                                           x_a[ self.x_index['density'][i], 0]) - flow[i] ])
            #     x_f[ self.x_index['density'][i], 0] = \
            #         x_a[ self.x_index['density'][i], 0 ] +\
            #         ( flow[i] - flow[i+1] + q_on)*self.dur_steps/self.len_cells
            # elif i not in self.x_index['onramp'].keys() and i in self.x_index['offramp'].keys():
            #     # The offramp state specifies the supply.
            #     # Hence the actual offramp flow depends on the demand in the cell
            #     # q_off = min( offramp_supply, Sending(i) - flow[i+1] )
            #     q_off = np.min([ x_a[ self.x_index['offramp'][i], 0],
            #                      self.__sending_flow(vm[i,0],beta[i,0],rhoc[i,0],qmax[i,0],
            #                                          x_a[ self.x_index['density'][i], 0]) - flow[i+1] ])
            #     x_f[ self.x_index['density'][i], 0] = \
            #         x_a[ self.x_index['density'][i], 0 ]+\
            #         ( flow[i] - flow[i+1] - q_off)*self.dur_steps/self.len_cells
            # elif i in self.x_index['onramp'].keys() and i in self.x_index['offramp'].keys():
            #     q_on = np.min([ x_a[ self.x_index['onramp'][i], 0],
            #                     self.__receiving_flow(wc[i,0],rhoc[i,0],rhom[i,0],qmax[i,0],
            #                                           x_a[ self.x_index['density'][i], 0]) - flow[i] ])
            #     q_off = np.min([ x_a[ self.x_index['offramp'][i], 0],
            #                      self.__sending_flow(vm[i,0],beta[i,0],rhoc[i,0],qmax[i,0],
            #                                          x_a[ self.x_index['density'][i], 0]) - flow[i+1] ])
            #     x_f[ self.x_index['density'][i], 0] = \
            #         x_a[ self.x_index['density'][i], 0 ]+\
            #         ( flow[i] - flow[i+1] + q_on - q_off)*self.dur_steps/self.len_cells

        # propagate the boundary flow using random walk.
        # add saturation
        q_in = x_a[self.x_index['qin'], 0]
        x_f[self.x_index['qin'], 0] = float(np.min([np.max([q_in, 0.0]),
                                                    self.qmax_cells[0, 0]]))
        q_out = x_a[self.x_index['qout'], 0]
        x_f[self.x_index['qout'], 0] = float(np.min([np.max([q_out, 0.0]),
                                                     self.qmax_cells[self.num_cells - 1, 0]]))

        # ======================================================== #
        # add noise to forecaste state
        x_f += np.matrix(np.random.multivariate_normal(
            np.zeros(self.dim_state), self.Q)).reshape((self.dim_state, 1))

        # # update the onramp flows by a random walk
        # if self.loc_onramp is not None:
        #     num_onramps = len(self.loc_onramp)
        #     e_onramp = np.matrix(np.random.multivariate_normal(np.zeros(num_onramps),
        #                                                                 self.Q['onramp'])).reshape((num_onramps,1))
        #
        #     for i in range(0,len(self.cell_onramp)):
        #         cell_id = self.cell_onramp[i]
        #         x_f[ self.x_index['onramp'][cell_id], 0] = \
        #             np.min( [ np.max( [ x_a[ self.x_index['onramp'][cell_id], 0 ] + e_onramp[i] , 0.0]),
        #                       self.qmax_cells[cell_id]])
        #
        #
        #
        # # update the offramp flow by a random walk
        # if self.loc_offramp is not None:
        #     num_offramps = len(self.loc_offramp)
        #     e_offramp = np.matrix(np.random.multivariate_normal(np.zeros(num_offramps),
        #                                                                 self.Q['offramp'])).reshape((num_offramps,1))
        #
        #     for i in range(0, len(self.cell_offramp)):
        #         cell_id = self.cell_offramp[i]
        #         x_f[ self.x_index['offramp'][cell_id], 0] = \
        #             np.min( [ np.max( [ x_a[ self.x_index['offramp'][cell_id], 0 ] + e_offramp[i], 0.0 ]),
        #                      self.qmax_cells[cell_id]])

        return x_f

    # @profile
    def h_obs(self, cur_time, e_id):
        """
        This is the observation function, only for one ensemble
        :param cur_time: Different sensors have different aggregation interval
        :param e_id: the ensemble id
        :return: y_f
        """
        end_time = cur_time

        y_f = np.matrix(np.zeros((self.dim_obs, 1), float))

        # the predicted value should be the aggregated value during the evolution process.
        # the forecast flow
        if 'flow' in self.y_index.keys():
            for s_id in self.y_index['flow'].keys():

                # if it is on freeway, use self.f_flow
                if s_id in self.sensors['freeway']:
                    cell_id = self.sensors[s_id]['cell']

                    # forecast flow data during time interval
                    # avoid doing matrix transformation to speed code up.
                    # current length of data
                    len_window = len(self.__f_flow['data'][e_id])

                    start_time = end_time - self.sensors[s_id]['aggregation_sec']

                    # get the index of the data instances to aggregate
                    index_start = np.searchsorted(self.__f_flow['time'][0:len_window], start_time, side='right')
                    index_end = np.searchsorted(self.__f_flow['time'][0:len_window], end_time, side='left')

                    # time_index = (start_time < self.__f_flow['time'][0:len_window]) & \
                    #              (self.__f_flow['time'][0:len_window]<= end_time)
                    # # get the index
                    # ele_index = [i for i,b in enumerate(time_index) if b]

                    flow_sum = 0.0
                    flow_num = 0
                    for i in range(index_start, index_end + 1):
                        flow_sum += self.__f_flow['data'][e_id][i][cell_id]
                        flow_num += 1

                    # tmp_flow = np.matrix( self.__f_flow['data'][e_id] ).T
                    # data_interval = tmp_flow[cell_id, time_index]

                    # save the average data as the forecast data
                    if flow_num != 0:
                        y_f[self.y_index['flow'][s_id]] = flow_sum / flow_num
                    else:
                        raise Exception('Error: take mean over empty array.')

                # if it is on onramp, use est_state_all[ x_index['onramp'][cell] ]
                elif s_id in self.sensors['onramp']:
                    cell_id = self.sensors[s_id]['cell']
                    # we can just reuse the same time grid to find the index in est_state_all
                    # __f_flow['time'] include all th time points. Extract the times point till now

                    start_time = end_time - self.sensors[s_id]['aggregation_sec']
                    time_index = (start_time < self.__f_flow['time']) & (self.__f_flow['time'] <= end_time)

                    # forecast onramp data during the time interval as the mean estimated value
                    # TODO, it should be depending on each ensemble, but here we are just using the mean estimated state.
                    data_interval = self.est_state_all[self.x_index['onramp'][cell_id], time_index]

                    # save the average data as the forecast data
                    y_f[self.y_index['flow'][s_id]] = data_interval.mean(1)

                # if it is on onramp, use est_state_all[ x_index['onramp'][cell] ]
                elif s_id in self.sensors['offramp']:
                    cell_id = self.sensors[s_id]['cell']
                    # we can just reuse the same time grid to find the index in est_state_all
                    start_time = end_time - self.sensors[s_id]['aggregation_sec']
                    time_index = (start_time < self.__f_flow['time']) & (self.__f_flow['time'] <= end_time)

                    # forecast onramp data during the time interval
                    # TODO, it should be depending on each ensemble
                    data_interval = self.est_state_all[self.x_index['offramp'][cell_id], time_index]

                    # save the average data as the forecast data
                    y_f[self.y_index['flow'][s_id]] = data_interval.mean(1)

                else:
                    raise Exception('Error: Flow sensors must locate on freeway grid or on/off ramps.')

        # the forecast velocity
        if 'speed' in self.y_index.keys():
            for s_id in self.y_index['speed'].keys():
                # the speed sensor must be on the freeway
                if s_id in self.sensors['freeway']:
                    cell_id = self.sensors[s_id]['cell']

                    len_window = len(self.__f_speed['data'][e_id])

                    start_time = end_time - self.sensors[s_id]['aggregation_sec']

                    # get the index of the data instances to aggregate
                    index_start = np.searchsorted(self.__f_speed['time'][0:len_window], start_time, side='right')
                    index_end = np.searchsorted(self.__f_speed['time'][0:len_window], end_time, side='left')

                    # time_index = (start_time < self.__f_speed['time'][0:len_window]) & \
                    #              (self.__f_speed['time'][0:len_window]<= end_time)
                    # ele_index = [i for i,b in enumerate(time_index) if b]

                    speed_sum = 0.0
                    speed_num = 0
                    for i in range(index_start, index_end + 1):
                        speed_sum += self.__f_speed['data'][e_id][i][cell_id]
                        speed_num += 1

                    # get the data during the interval
                    # tmp_speed = np.matrix( self.__f_speed['data'][e_id] ).T
                    # len_window = tmp_speed.shape[1]
                    # data_interval = tmp_speed[cell_id, time_index]

                    # the forecast velocity is the average value
                    if speed_num != 0:
                        y_f[self.y_index['speed'][s_id]] = speed_sum / speed_num
                    else:
                        raise Exception('Error: take mean over empty array')

                else:
                    raise Exception('Error: Speed sensors must locate on the freeway grid.')

        # the forecast travel time
        # TODO: need to update the bluetooth observation equation, now it is an average of the snap shot travel time during the interval
        if 'travel_time' in self.y_index.keys():
            for s_id in self.y_index['travel_time'].keys():

                start_cell, end_cell = self.sensors[s_id]['cell']
                # still use the __f_flow to get the time grid index
                start_time = end_time - self.sensors[s_id]['aggregation_sec']
                time_index = (start_time < self.__f_flow['time']) & (self.__f_flow['time'] <= end_time)

                # get the traffic density data from the interval
                data_interval = self.est_state_all[
                                self.x_index['density'][start_cell]: self.x_index['density'][start_cell],
                                time_index]

                travel_time = []
                speed_cells = []
                # compute a snap-shot travel time using each density profile
                for j in range(0, data_interval.shape[1]):
                    for i in range(0, data_interval.shape[0]):
                        cell_id = start_cell + i
                        vel = self.__rho2v(self.vm_cells[cell_id, 0], self.beta_cells[cell_id, 0],
                                           self.rhoc_cells[cell_id, 0], self.wc_cells[cell_id, 0],
                                           data_interval[i, j])
                        if vel == 0:
                            raise Exception('Error: Got negative speed.')
                        else:
                            speed_cells.append(float(vel))

                    speed_cells = np.array(speed_cells)

                    # get the length of cells
                    L_cells = self.len_cells[start_cell: end_cell]

                    # append this snap shot travel time
                    travel_time.append(np.sum(L_cells / speed_cells))

                # compute the mean travel time as the forecast travel time
                travel_time = np.array(travel_time)
                y_f[self.y_index['travel_time'][s_id]] = travel_time.mean()

        return y_f

    def __sending_flow(self, vm, beta, rhoc, qmax, rho):
        """
        This function computes the sending flow of a cell given the density, with saturation at qmax using QLFD
        :param vm: para for freeflow regime
        :param beta: para for freeflow regime
        :param rhoc: para for freeflow regime
        :param qmax: para for freeflow regime
        :param rho: the density to compute
        :return: float, the flow
        """
        if rho < 0.0:
            return 0.0
        elif rho < rhoc:
            return float(vm * rho - vm * (rho ** 2) / beta)
        else:
            return float(qmax)

    def __receiving_flow(self, w, rhoc, rhom, qmax, rho):
        """
        This function computes the receiving flow of a cell given the density, with saturation at qmax using QLFD
        :param w: parameter for congested regime
        :param rhom: parameter for congested regime
        :param rhoc: parameter for congested regime
        :param qmax: parameter for congested regime
        :param rho: density
        :return: float, flow
        """
        if rho > rhoc:
            return float(w * (rho - rhom))
        else:
            return float(qmax)

    def __rho2v(self, vm, beta, rhoc, w, rho):
        """
        This function computes the speed of a cell given the density, with saturation at qmax using QLFD
        :param vm: para for freeflow regime
        :param beta: para for freeflow regime
        :param rhoc: para for freeflow regime
        :param qmax: para for freeflow regime
        :param rho: the density to compute
        :return: float, the flow
        """
        if rho < 0:
            return float(vm)
        elif rho <= rhoc:
            return float(vm - vm * rho / beta)
        else:
            rhom = rhoc - (vm * rhoc - vm * (rhoc ** 2) / beta) / w
            # print('rho {0}; rhoc {1}'.format(rho, rhoc))
            return float(w * (rho - rhom) / rho)

    def __q2v_ff(self, vm, beta, q):
        """
        This function converts the flow to speed in the freeflow regime
        :param vm: para
        :param beta: para
        :param q: para
        :return: speed
        """
        return float((vm * beta - np.sqrt(np.power(vm * beta, 2) - 4 * vm * beta * q)) / (2 * vm))

    def __q2v_cf(self, w, rhom, q):
        """
        This function converts the flow to speed in the congested flow regime
        :param w: para, should be negative
        :param rhom: para
        :param q: para
        :return: speed
        """
        return float(q / (rhom + q / w))

    def plot_normalized_density_est(self):
        """
        This function plots the normalized estimated density profile in the entire time horizon
        :return:
        """
        # The true density is in self.x_all, num_cells x sim_steps.
        # the first num_cells row is the density
        dens = np.flipud(self.est_density[0:self.num_cells, :])

        # normalize the density [0, rho_max] to [0,1] for all cells
        rho_max = np.flipud(np.matrix(self.rhomax_cells).reshape((self.num_cells, 1)))

        norm_dens = dens / rho_max

        fig = plt.figure(figsize=(10, 10), dpi=100)
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        im = ax.imshow(norm_dens, cmap=plt.get_cmap('jet'),
                       interpolation='nearest',
                       vmin=0, vmax=1)
        ax.autoscale(False)

        # ==============================================================
        # plot the sensors, and on/off ramps
        # plot onramps
        onramp_ytick = {}
        onramp_ytick['tick'] = []
        onramp_ytick['label'] = []
        if self.cell_onramp is not None:
            for cell_id in self.cell_onramp:
                # the onramp will be drawn as a line in the middle of the cell
                val = self.num_cells - 1 - cell_id
                onramp_ytick['tick'].append(val)
                onramp_ytick['label'].append('Onramp {0}'.format(cell_id))
                ax.plot([-0.5, self.num_steps - 0.5], [val, val], color='k', linewidth=1)

        # plot offramps
        offramp_ytick = {}
        offramp_ytick['tick'] = []
        offramp_ytick['label'] = []
        if self.cell_offramp is not None:
            for cell_id in self.cell_offramp:
                # the offramp will be drawn as a line in the middle of the cell
                val = self.num_cells - 1 - cell_id
                offramp_ytick['tick'].append(val)
                offramp_ytick['label'].append('Offramp {0}'.format(cell_id))
                ax.plot([-0.5, self.num_steps - 0.5], [val, val], color='k', linewidth=1)

        # plot flow sensors
        flow_sensor_ytick = {}
        flow_sensor_ytick['tick'] = []
        flow_sensor_ytick['label'] = []
        if len(self.sensors['freeway']) != 0:
            for s_id in self.sensors['freeway']:
                # the flow sensor will be drawn as a line in the begining of the cells.
                cell_id = self.sensors[s_id]['cell']

                val = self.num_cells - 0.5 - cell_id
                flow_sensor_ytick['tick'].append(val)
                flow_sensor_ytick['label'].append('Flow Sensor {0}'.format(cell_id))
                ax.plot([-0.5, self.num_steps - 0.5], [val, val], color='k', linewidth=1)

        # set y ticks
        ax.set_yticks(onramp_ytick['tick'] + offramp_ytick['tick'] + flow_sensor_ytick['tick'])
        ax.set_yticklabels(onramp_ytick['label'] + offramp_ytick['label'] + flow_sensor_ytick['label'])

        ax.set_title('Estimated density profile')
        cax = fig.add_axes([0.95, 0.25, 0.01, 0.5])
        fig.colorbar(im, cax=cax, orientation='vertical')
        plt.draw()

    def plot_one_state_evolution(self, value, option='density'):
        """
        This function plots the estimated evolution of one state.
        :param value: cell_id
        :param option: either density of flow
        :return: one figure
        """

        fig = plt.figure(figsize=(16, 8), dpi=100)
        ax = fig.add_subplot(111)
        if option == 'density':
            ax.plot(np.squeeze(np.array(self.est_density[value, :])), 'r--', label='Estimated')
        elif option == 'flow':
            if value == 0:
                # plot qin
                index = self.x_index['qin']
            elif value == self.num_cells:
                # plot qout
                index = self.x_index['qout']
            elif value in self.x_index['onramp'].keys():
                # plot onramp flow
                index = self.x_index['onramp'][value]
            elif value in self.x_index['offramp'].keys():
                # plot offramp flow
                index = self.x_index['offramp'][value]
            else:
                raise Exception('KeyError: value must be the cell id')

            ax.plot(np.squeeze(np.array(self.est_state_all[index, :])), 'r--', label='Estimated')

        plt.title('{0} at {1}'.format(option, value))
        plt.xlabel('Time (step)')
        plt.ylabel('Value')

        plt.grid(True)
        plt.legend()

        plt.draw()

    def compare_one_state_evolution(self, value, true_state, option='density'):
        """
        This function plots the estimated evolution of one state specified by value
        :param value: cell_id
        :param true_state: array, the true state in the cell_id cell
        :param option: either density of flow
        :return: one figure
        """

        fig = plt.figure(figsize=(16, 8), dpi=100)
        ax = fig.add_subplot(111)
        if option == 'density':
            ax.plot(np.squeeze(np.array(self.est_density[value, :])), 'r--', label='Estimated')
        elif option == 'flow':
            if value == 0:
                # plot qin
                index = self.x_index['qin']
            elif value == self.num_cells:
                # plot qout
                index = self.x_index['qout']
            elif value in self.x_index['onramp'].keys():
                # plot onramp flow
                index = self.x_index['onramp'][value]
            elif value in self.x_index['offramp'].keys():
                # plot offramp flow
                index = self.x_index['offramp'][value]
            else:
                raise Exception('KeyError: value must be the cell id')

            ax.plot(np.squeeze(np.array(self.est_state_all[index, :])), 'r--', label='Estimated')
            ax.plot(np.squeeze(np.array(true_state)), 'b', label='True')

        plt.title('{0} at {1}'.format(option, value))
        plt.xlabel('Time (step)')
        plt.ylabel('Value')

        plt.grid(True)
        plt.legend()

        plt.draw()

    def __limit_value(self, value, v_range):
        """
        This function limits the value in the range. A saturation function.
        :param value: the value to confine
        :param v_range: the range [min, max]
        :return: the limited value
        """
        if np.isnan(value):
            print('Warning: trying to limit nan value in range {0}'.format(v_range))
            return value

        return np.min([v_range[1], np.max([value, v_range[0]])])

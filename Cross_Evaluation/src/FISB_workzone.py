import numpy as np
from Estimator import Estimator
from FISB import FISB

__author__ = 'Juan Carlos Martinez'


class FISB_workzone(Estimator, FISB):

    def __init__(self, time_grid=None, space_grid=None, sensors=None, loc_onramp=None, loc_offramp=None,
                 queue_threshold=17.88, omg_c=-15, omg_f=80, sigma_fac=15, tau_fac=15, v_crit=60, delta_v=20, lag=0,
                 time_diff_cutoff=300):
        """
        constructor
        :param time_grid: speed estimation time grid
        :param space_grid: speed estimation space grid
        :param sensors: sensor info as a dictionary
        :param loc_onramp: on ramp info
        :param loc_offramp: off ramp info
        :param queue_threshold: threshold to consider queue
        :param omg_c: wave propagation speed in congestion [km/h]
        :param omg_f: wave propagation speed in free flow [km/h]
        :param sigma_fac: factor of sensor spacing for smoothing band []
        :param tau_fac: factor of sensor time aggregation for smoothing band []
        :param v_crit: critical speed for traffic state shift [km/h]
        :param delta_v: width of the traffic state shift smoothing band [km/h]
        :param lag: allowed time lag for inverse speed estimation [s]
        :param time_diff_cutoff: time difference cutoff after which the sensor reading and the
                                 point being estimated are too far apart in time to be considered
                                 for the inverse speed estimation [s]
        """

        # save parameters
        self.num_cells = len(space_grid)-1
        self.sim_steps = len(time_grid)-1
        self.len_cell = space_grid[1] - space_grid[0]
        self.step_dur = time_grid[1] - time_grid[0]
        self.queue_threshold = queue_threshold
        self.time_grid = time_grid
        self.space_grid = space_grid
        self.speed_data = None

        # construct parent classes
        Estimator.__init__(self, time_grid, space_grid, loc_onramp, loc_offramp, sensors, queue_threshold)

        # save sorted sensor locations along main road
        sensors_space_grid = []
        for i in self.sensors.keys():
            if i != 'onramp' and i != 'offramp' and i != 'freeway':
                sensors_space_grid.append(self.sensors[i]['loc'])
        sensors_space_grid = sorted(sensors_space_grid)
        self.sensors_space_grid = sensors_space_grid

        FISB.__init__(self, omg_c, omg_f, sigma_fac, tau_fac, v_crit, delta_v, lag, time_diff_cutoff,
                      self.sensors_space_grid)

    def generate_estimates(self):
        """
        this function generates the estimates for speed, travel time and length of queue
        using the FISB filtering method. the measurement data needs to be set for the
        estimator class beforehand
        :return:
        """

        # set measurement data for FISB class
        if self.speed_data is not None:
            FISB.set_fisb_data(self, self.speed_data)
        else:
            raise Exception('Set the measurement data...')

        for t_idx in range(0, self.sim_steps):


            tt_at_time = 0
            for x_idx in range(0, self.num_cells):
                # speed map
                speed = 1/self.get_w(self.time_grid[t_idx], self.space_grid[x_idx])
                self.est_speed[x_idx, t_idx] = speed
                tt_at_time += self.len_cells/speed
            # travel time
            self.est_traveltime[t_idx] = tt_at_time
            # queue
            if not np.isnan(self.est_speed[:, t_idx]).any():
                index = (self.est_speed[:, t_idx] <= self.queue_threshold)
                if sum(index) <= 2:
                    # all uncongested
                    self.est_queue[t_idx] = 0
                else:
                    self.est_queue[t_idx] = self.len_cell * (self.num_cells - np.argmax(index))
            else:
                self.est_queue[t_idx] = float('nan')

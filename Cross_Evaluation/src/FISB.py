import math
import numpy as np

__author__ = 'Juan Carlos Martinez'

"""
this is an implementation (capable of on-line) of the filtered inverse-speed based algorithm used for speed map
reconstruction

methods:
    get_w(t,x): this function estimates the inverse speed at location (t,x)
    set_fisb_data(speed_data): this function saves the speed data required for all estimations
"""


class FISB:

    def __init__(self, omg_c, omg_f, sigma_fac, tau_fac, v_crit, delta_v, lag, time_diff_cutoff, sensors_space_grid):
        """
        constructor
        :param omg_c: wave propagation speed in congestion [km/h]
        :param omg_f: wave propagation speed in free flow [km/h]
        :param sigma_fac: factor of sensor spacing for smoothing band []
        :param tau_fac: factor of sensor time aggregation for smoothing band []
        :param v_crit: critical speed for traffic state shift [km/h]
        :param delta_v: width of the traffic state shift smoothing band [km/h]
        :param lag: allowed time lag for inverse speed estimation [s] (this is 0 for on-line implementation)
        :param time_diff_cutoff: time difference cutoff after which the sensor reading and the
                                 point being estimated are too far apart in time to be considered
                                 for the inverse speed estimation [s]
        :param sensors_space_grid: space grid containing the sensors locations
        :return:
        """

        # store constructor input in class with required unit conversions
        self.__omg_c = omg_c*1000/3600  # conversion factor to [m/s]
        self.__omg_f = omg_f*1000/3600  # conversion factor to [m/s]
        self.__sigma_fac = sigma_fac
        self.__tau_fac = tau_fac
        self.__v_crit = v_crit*1000/3600  # conversion factor to [m/s]
        self.__delta_v = delta_v*1000/3600  # conversion factor to [m/s]
        self.__lag = lag
        self.__time_diff_cutoff = time_diff_cutoff
        self.__sensors_space_grid = sensors_space_grid
        self.__num_sensors = len(self.__sensors_space_grid)
        self.__sigma = None
        self.__tau = None
        self.__sensors_time_grid = None
        self.__sensors_speed_data = None

    def set_fisb_data(self, speed_data):
        """
        this function saves the speed data in the FISB class
        :param speed_data: dictionary holding speed data
        :return:
        """
        self.__sensors_time_grid = np.array(speed_data['time'])
        self.__sigma = self.__sigma_fac*(self.__sensors_space_grid[-1] - self.__sensors_space_grid[-2])
        self.__tau = self.__tau_fac*(self.__sensors_time_grid[-1] - self.__sensors_time_grid[-2])
        self.__sensors_speed_data = np.array(speed_data['data'])

    def get_w(self, t, x):
        """
        this function returns the inverse speed estimate at location (t,x) by weighting the
        estimates in the congested and free flow states
        :param t: temporal location of the inverse speed estimate point
        :param x: spatial location of the inverse speed estimate point
        :param sigma: width for spatial smoothing [m]
        :return:
        """
        # set sigma based on the current configuration
        w_c = self.__get_w(t, x, 'cong')
        w_f = self.__get_w(t, x, 'free')
        gamma = 0.5*(1 + math.tanh((self.__v_crit - min(1/w_c, 1/w_f))/self.__delta_v))
        
        # set sigma back to None for the next configuration
        return gamma*w_c + (1-gamma)*w_f

    def __get_w(self, t, x, state):
        """
        this function returns the inverse speed estimate at location (t,x) given a state ('cong' for congestion,
        'free' for free flow)
        :param t: temporal location for the inverse speed estimate point
        :param x: spatial location for the inverse speed estimate point
        :param state: state of traffic, 'cong' for congested, 'free' for free-flowing
        :return w: inverse speed estimate at location (t,x) given traffic state
        """
        # w = self.__get_exponential_weighted_avg(t, x, omg)
        if state == 'cong':
            return self.__get_exponential_weighted_avg(t, x, self.__omg_c)
        elif state == 'free':
            return self.__get_exponential_weighted_avg(t, x, self.__omg_f)
        else:
            raise Exception('Error: Wrong traffic state...')

    def __get_exponential_weighted_avg(self, t, x, omg):
        """
        this function computes the exponential weighted average of the inverse speeds at the considered
        sensor locations in time and spate to estimate the inverse speed at location (t,x) given a wave propagation
        speed omg
        :param t: temporal location for inverse speed estimate point
        :param x: spatial location for inverse speed estimate point
        :param omg: wave propagation speed
        :return w: inverse speed estimate from weighted average
        """

        # truncate time and speed data for the ones that are of relevance
        time_idxs = np.intersect1d(np.where(self.__sensors_time_grid >= t - self.__time_diff_cutoff),
                                   np.where(self.__sensors_time_grid <= t + self.__lag))

        # set lists that will hold the weights and speed of each point used
        exponential_weights = []
        speeds = []

        # set weights and speeds
        for t_idx in time_idxs:
            for x_idx in range(0, self.__num_sensors):
                delta_t = self.__sensors_time_grid[t_idx] - t
                delta_x = self.__sensors_space_grid[x_idx] - x
                # add to lists only if there is available speed data
                if (not math.isinf(self.__sensors_speed_data[x_idx, t_idx]) and
                        not math.isnan(self.__sensors_speed_data[x_idx, t_idx])):
                    exponential_weights.append(
                        math.e**(-(abs(delta_x)/self.__sigma) - (abs(delta_t - delta_x/omg)/self.__tau)))
                    speeds.append(self.__sensors_speed_data[x_idx, t_idx])

        if not exponential_weights:
            return float('nan')
        else:
            # compute weighted exponential average
            exponential_weights = np.array(exponential_weights)
            return np.sum(np.true_divide(exponential_weights, np.array(speeds)))/np.sum(exponential_weights)



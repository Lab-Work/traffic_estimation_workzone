import bisect
from copy import deepcopy
import matplotlib.pyplot as plt
import numpy as np

__author__ = 'Yanning Li'
"""
This is the abstract class of estimator. It defines the input and output of the estimators. All later developed estimators
should use this as its super class
"""


class Estimator:
    """
    Abstract estimator class.
    """

    def __init__(self, time_grid=None, space_grid=None,
                 loc_onramp=None, loc_offramp=None,
                 sensors=None,
                 queue_threshold=17.88):
        """
        Constructor of the estimator. Any estimator should have at least the above inputs.
            Subclasses must initialize this first.
        - The estimator is designed to run on a discretized road with the finest grid from which we extract
            the true states.
            - This grid is later used for visualization, i.e., estimated constant density and speed for each cell.
            - This grid is used for computing a measure of the performance against the true states which are defined on
                the grid.
        - All sensors are assumed to be placed on the grid.
        - The aggregation interval of sensors = n * dt, where n is an integer and dt is the time grid resolution.
        :param time_grid: list of float seconds, num_steps+1
        :param space_grid: list of float m, num_cells+1 since 0 m is included
        :param sensors: keys: sensor_id
                [sensor_id]['loc'] : the location
                           ['cell'] : the cell
        :param loc_onramp: list of floats, meters
        :param loc_offramp: list of floats, meters
        :param sensors: dict, information for each sensor
                sensors[sensor_id]: dict with keys:
                    'loc', 'aggregation_sec','flow_std','speed_std','traveltime_std'
        :param queue_threshold: m/s, the threshold to determine congestion. If speed in a cell is below, then congested.
        :return:
        """

        # =========================================================================
        # The finest grid for visualization and evaluation.
        # convert to list to get the index easier
        if type(space_grid) is np.ndarray:
            self.space_grid = space_grid.tolist()
        elif type(space_grid) is list:
            self.space_grid = space_grid
        else:
            raise Exception('Error: space grid should be a list of floats')
        if type(time_grid) is np.ndarray:
            self.time_grid = time_grid.tolist()
        elif type(time_grid) is list:
            self.time_grid = time_grid
        else:
            raise Exception('Error: time grid should be a list of floats')

        # =========================================================================
        # Discretize the road into cell using the time-space grid. The algorithms may not necessarily need to run on the
        # grid, but the estimates are saved on this grid resolution.
        self.num_cells = len(self.space_grid) - 1
        self.num_steps = len(self.time_grid) - 1
        self.len_cells = self.space_grid[1] - self.space_grid[0]
        self.dur_steps = self.time_grid[1] - self.time_grid[0]

        self.loc_onramp = loc_onramp
        self.loc_offramp = loc_offramp
        self.cell_onramp, self.cell_offramp = self.__discretize_road()

        # =========================================================================
        # Keep a copy of the sensor network information, and place the sensors on the discretized road.
        self.sensors = self.__place_sensors(sensors)

        # =========================================================================
        # The queue threshold
        self.queue_threshold = queue_threshold

        # =========================================================================
        # Initialize the measurement data. We assume this is a batch estimator, i.e., get all the data from files.
        self.speed_data = None
        self.flow_data = None
        self.traveltime_data = None

        # =========================================================================
        # The output of each estimator should be the follows
        # Estimated density profile for entire time-space horizon on the finest grid
        self.est_speed = np.matrix(np.zeros((len(self.space_grid) - 1, len(self.time_grid) - 1))).astype(float)
        self.est_density = np.matrix(np.zeros((len(self.space_grid) - 1, len(self.time_grid) - 1))).astype(float)
        self.est_queue = np.array(np.zeros(len(self.time_grid) - 1).astype(float))
        self.est_traveltime = np.array(np.zeros(len(self.time_grid) - 1).astype(float))

    def __discretize_road(self):
        """
        This function converts the input road information to discretized cells
        :return: the cell ids of onramps and offrapms
        """
        # get the index of the cells with onramp or offramp
        if self.loc_onramp is not None:
            cell_onramp = [bisect.bisect(self.space_grid, i) - 1 for i in self.loc_onramp]
        else:
            cell_onramp = None
        if self.loc_offramp is not None:
            cell_offramp = [bisect.bisect(self.space_grid, i) - 1 for i in self.loc_offramp]
        else:
            cell_offramp = None

        return cell_onramp, cell_offramp

    # process and place the sensor network on grids
    def __place_sensors(self, sensors):
        """
        This function places sensors on the freeway, onramps, or offramps depending on their location.
        It also creates another key which maps the location to the cell
        :param sensors:
        :return:
        """
        if sensors is not None:
            # wireless sensor network
            wsn = deepcopy(sensors)

            freeway = []
            onramp = []
            offramp = []
            # process each sensor
            for s_id in wsn.keys():
                loc = wsn[s_id]['loc']
                # the location of flow and speed sensors are float values, the location of the BT sensor is a tuple
                if type(loc) is float or type(loc) is int:
                    if loc in self.space_grid:
                        # on the freeway
                        freeway.append(s_id)
                        wsn[s_id]['cell'] = self.space_grid.index(loc)
                    elif self.loc_onramp is not None and loc in self.loc_onramp:
                        # on the onramp
                        onramp.append(s_id)
                        wsn[s_id]['cell'] = bisect.bisect(self.space_grid, loc) - 1
                    elif self.loc_offramp is not None and loc in self.loc_offramp:
                        offramp.append(s_id)
                        wsn[s_id]['cell'] = bisect.bisect(self.space_grid, loc) - 1
                    else:
                        print('Location of sensor:{0}'.format(loc))
                        raise Exception('Error: flow sensors must locate on the grid or on/off ramps.')

                elif type(loc) is tuple:
                    loc1, loc2 = loc
                    wsn[s_id]['cell'] = (self.space_grid.index(loc1), self.space_grid.index(loc2))

                else:
                    raise Exception('Error: the location of the sensors must be int, float, or tuple of floats')

            # add the freeway, onramp, and offramp sensor information to sensor
            wsn['freeway'] = freeway
            wsn['onramp'] = onramp
            wsn['offramp'] = offramp

            return wsn

        else:
            return None

    # set the measurement data
    def set_meas_data(self, flow_data=None, speed_data=None, BT_data=None):
        """
        This function sets saves all the data in properties. It assumes the data for the entire time horizon is obtained
            prior of running the offline estimator.
        :param speed_data: Dict with keys: 'time', 'data'
                           Dict['time'] = 1d np.array of the timestamps of each speed data.
                           Dict['sensors'] = 1d np.array of sensor id, corresponds to the rows of data
                           Dict['data'] = 2d np.array of speed in m/s. DIM = len(loc_speed_sensors) x len(Dict['time'])
                                          Each row is the speed for each sensor in the same order of loc_speed_sensors
                                          Each column is a time slice.
        :param flow_data:  Dict with keys: 'time', 'data'
                           Dict['time'] = 1d np.array of the timestamps of each flow data.
                           Dict['sensors'] = 1d np.array of sensor id, corresponds to the rows of data
                           Dict['data'] = 2d np.array of flow in this interval. DIM = len(loc_flow_sensors) x len(Dict['time'])
                                          Each row is the flow for each sensor in the same order of loc_flow_sensors
                                          Each column is a time slice.
        :param BT_data:    Dict with keys: 'time', 'data'
                           Dict['time'] = 1d np.array of the timestamps of each BT data.
                           Dict['sensors'] = 1d np.array of sensor id, corresponds to the rows of data
                           Dict['data'] = 2d np.array of travel time in s. DIM = len(loc_BT_sensors) x len(Dict['time'])
                                          Each row is the travel time for each BT pair in the same order of loc_BT_sensors
                                          Each column is a time slice.
        :return:
        """
        self.speed_data = speed_data
        self.flow_data = flow_data
        self.traveltime_data = BT_data

    # Visualization
    def plot_est_speed(self, unit, limit):
        """
        This function plots the estimated speed profile over the entire time space domain in the specified unit
        :param unit: 'metric', 'imperial'; respectively 'm, s, m/s', and 'mile, hour, mph'
        :param limit: The limit of the colorbar in above units
        :return: A figure profile with x-axis being the time, and y-axis being the space. (flow direction upwards)
        """
        if unit == 'metric':
            # all internal values are in metric, so plot directly
            speed = np.flipud(self.est_speed)
            unit_str = 'm/s'
        elif unit == 'imperial':
            speed = np.flipud(self.__metric2imperial(self.est_speed, 'speed'))
            limit = self.__metric2imperial(np.array(limit), 'speed')
            unit_str = 'mph'
        else:
            raise Exception('Error: Unrecognized unit for plotting speed.')

        fig = plt.figure(figsize=(10, 10), dpi=100)
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        im = ax.imshow(speed, cmap=plt.get_cmap('jet'),
                       interpolation='nearest',
                       vmin=limit[0], vmax=limit[1])
        ax.autoscale(False)

        # ==============================================================
        # plot the sensors, and on/off ramps
        # plot onramps
        onramp_ytick = {'tick': [], 'label': []}
        if self.cell_onramp is not None:
            for cell_id in self.cell_onramp:
                # the onramp will be drawn as a line in the middle of the cell
                val = self.num_cells - 1 - cell_id
                onramp_ytick['tick'].append(val)
                onramp_ytick['label'].append('Onramp {0}'.format(cell_id))
                ax.plot([-0.5, self.num_steps - 0.5], [val, val], color='k', linewidth=1)

        # plot offramps
        offramp_ytick = {'tick': [], 'label': []}
        if self.cell_offramp is not None:
            for cell_id in self.cell_offramp:
                # the offramp will be drawn as a line in the middle of the cell
                val = self.num_cells - 1 - cell_id
                offramp_ytick['tick'].append(val)
                offramp_ytick['label'].append('Offramp {0}'.format(cell_id))
                ax.plot([-0.5, self.num_steps - 0.5], [val, val], color='k', linewidth=1)

        # plot flow sensors
        sensor_ytick = {'tick': [], 'label': []}
        if len(self.sensors['freeway']) != 0:
            for s_id in self.sensors['freeway']:
                # the flow sensor will be drawn as a line in the beginning of the cells.
                cell_id = self.sensors[s_id]['cell']

                val = self.num_cells - 0.5 - cell_id
                sensor_ytick['tick'].append(val)
                sensor_ytick['label'].append('{0}'.format(s_id))
                ax.plot([-0.5, self.num_steps - 0.5], [val, val], color='k', linewidth=1)

        # set y ticks
        ax.set_yticks(onramp_ytick['tick'] + offramp_ytick['tick'] + sensor_ytick['tick'])
        ax.set_yticklabels(onramp_ytick['label'] + offramp_ytick['label'] + sensor_ytick['label'])

        ax.set_title('Estimated speed profile ({0})'.format(unit_str))
        plt.xlabel('Time')
        plt.ylabel('Space, traffic direction $\Rightarrow$')
        cax = fig.add_axes([0.95, 0.25, 0.01, 0.5])
        fig.colorbar(im, cax=cax, orientation='vertical')
        plt.draw()

    # Visualization
    def plot_est_density(self, unit, limit):
        """
        This function plots the estimated density profile over the entire time space domain in the specified unit
        :param unit: 'metric', 'imperial'; respectively 'm, s, m/s', and 'mile, hour, mph'
        :param limit: The limit of the colorbar in above units
        :return: A figure profile with x-axis being the time, and y-axis being the space. (flow direction upwards)
        """
        if unit == 'metric':
            # all internal values are in metric, so plot directly
            density = np.flipud(self.est_density)
            unit_str = 'veh/m'
        elif unit == 'imperial':
            density = np.flipud(self.__metric2imperial(self.est_density, 'density'))
            limit = self.__metric2imperial(np.array(limit), 'density')
            unit_str = 'veh/mile'
        else:
            raise Exception('Error: Unrecognized unit for plotting density.')

        fig = plt.figure(figsize=(10, 10), dpi=100)
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        im = ax.imshow(density, cmap=plt.get_cmap('jet'),
                       interpolation='nearest',
                       vmin=limit[0], vmax=limit[1])
        ax.autoscale(False)

        # ==============================================================
        # plot the sensors, and on/off ramps
        # plot onramps
        onramp_ytick = {'tick': [], 'label': []}
        if self.cell_onramp is not None:
            for cell_id in self.cell_onramp:
                # the onramp will be drawn as a line in the middle of the cell
                val = self.num_cells - 1 - cell_id
                onramp_ytick['tick'].append(val)
                onramp_ytick['label'].append('Onramp {0}'.format(cell_id))
                ax.plot([-0.5, self.num_steps - 0.5], [val, val], color='k', linewidth=1)

        # plot offramps
        offramp_ytick = {'tick': [], 'label': []}
        if self.cell_offramp is not None:
            for cell_id in self.cell_offramp:
                # the offramp will be drawn as a line in the middle of the cell
                val = self.num_cells - 1 - cell_id
                offramp_ytick['tick'].append(val)
                offramp_ytick['label'].append('Offramp {0}'.format(cell_id))
                ax.plot([-0.5, self.num_steps - 0.5], [val, val], color='k', linewidth=1)

        # plot flow sensors
        sensor_ytick = {'tick': [], 'label': []}
        if len(self.sensors['freeway']) != 0:
            for s_id in self.sensors['freeway']:
                # the flow sensor will be drawn as a line in the beginning of the cells.
                cell_id = self.sensors[s_id]['cell']

                val = self.num_cells - 0.5 - cell_id
                sensor_ytick['tick'].append(val)
                sensor_ytick['label'].append('{0}'.format(s_id))
                ax.plot([-0.5, self.num_steps - 0.5], [val, val], color='k', linewidth=1)

        # set y ticks
        ax.set_yticks(onramp_ytick['tick'] + offramp_ytick['tick'] + sensor_ytick['tick'])
        ax.set_yticklabels(onramp_ytick['label'] + offramp_ytick['label'] + sensor_ytick['label'])

        ax.set_title('Estimated density profile ({0})'.format(unit_str))
        plt.xlabel('Time')
        plt.ylabel('Space, traffic direction $\Rightarrow$')
        cax = fig.add_axes([0.95, 0.25, 0.01, 0.5])
        fig.colorbar(im, cax=cax, orientation='vertical')
        plt.draw()

    # Visualization
    def plot_est_queue(self, unit):
        """
        This function plots the estimated queue length over the entire time horizon in the specified unit
        :param unit: 'metric', 'imperial'; respectively 'm, s, m/s', and 'mile, hour, mph'
        :return: A figure plot of queue length
        """
        if unit == 'metric':
            # all internal values are in metric, so plot directly
            queue = self.est_queue
            unit_str = 'm'
        elif unit == 'imperial':
            queue = self.__metric2imperial(self.est_queue, 'distance')
            unit_str = 'mile'
        else:
            raise Exception('Error: Unrecognized unit for plotting density.')

        fig = plt.figure(figsize=(10, 10), dpi=100)
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        ax.plot(self.time_grid[1:], queue, linewidth=1)
        plt.title('Estimated queue length')
        plt.xlabel('Time s')
        plt.ylabel('Length of queue ({0})'.format(unit_str))

        plt.draw()

    # Visualization
    def plot_est_traveltime(self, unit):
        """
        This function plots the estimated travel time over the entire time horizon in the specified unit
        :param unit: 'metric', 'imperial'; respectively 'm, s, m/s', and 'mile, hour, mph'
        :return: A figure plot of queue length
        """

        fig = plt.figure(figsize=(10, 10), dpi=100)
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        ax.plot(self.time_grid[1:], self.est_traveltime, linewidth=1)
        plt.title('Estimated queue length')
        plt.xlabel('Time s')
        plt.ylabel('Length of queue (s)')

        plt.draw()

    # utility functions
    @staticmethod
    def __metric2imperial(value=np.zeros((1, 1)), option='speed'):
        """
        A utility function which converts the metric (m, s, m/s) to imperial (mile, hour, m/h)
        :param value: float, np.array, or np.matrix. the to be converted value
        :param option: 'speed', 'density'
        :return: converted value
        """
        if type(value) is float or type(value) is np.ndarray or type(value) is np.matrix:
            if option == 'speed':
                return value * 3600.0 / 1609.34
            elif option == 'density':
                return value * 1609.34
            elif option == 'distance':
                return value / 1609.34
            else:
                raise Exception('Error: Unrecognized unit conversion option.')
        else:
            raise Exception('Error: Unrecognized value type for unit conversion.')

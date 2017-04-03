import os
import sys
from collections import OrderedDict

import numpy as np

__author__ = 'Yanning Li and Juan Carlos Martinez'
"""
This script generates the configuration file for estimation. The configuration file consists of blocks separated by a
blank line. Each block includes a specific sensor network configuration and algorithms to be used.
"""

# =========================================================================================== #
# ============================== Generate I80 Configuration ================================= #
# =========================================================================================== #

directory = os.path.dirname(os.path.realpath(__file__))
log_dir = '/data_fast/Yanning_workzone/'
f = open(log_dir+'I80_configurations_input.txt', 'w+')

# ===================================================
# Configure the I80 network
start_loc = 300
end_loc = 8300
cell_length = 200

# the order of main freeway sections
fwy_secs = [21216, 24671, 21217, 37745, 21419,
            37742, 15814, 17054, 15589, 21104,
            3412, 6831, 3399, 6840, 3400,
            6843, 3401, 6849, 40842, 40963,
            40962, 6799, 3248, 13019, 13008]

# the length of sections
len_secs = [78.0, 24.69, 73.26, 32.08, 686.73,
            34.11, 1407.16, 18.89, 991.51, 18.6,
            3330.64, 27.85, 443.20, 21.22, 58.37,
            19.19, 310.43, 25.32, 877.76, 20.22,
            537.52, 50.03, 683.021, 26.70, 637.068]
# the offset of freeway sections. that is the absolute distance of the section from the entrance of freeway
offset_secs = np.concatenate([[0], np.cumsum(np.array(len_secs))])


# ================================================================ #
# ======================= A few utility functions ================ #
# ================================================================ #
def __get_relative_loc(dist):
    """
    This function returns the relative location on the freeway using dist to the entrace
    :param dist: meters, to the entrance
    :return: sect, rel_dist(m)
    """
    sect = 0
    _rel_dist = 0.0

    for idx, offset in enumerate(offset_secs):

        if dist < offset:
            # found
            return sect, _rel_dist
        else:
            # keep updating
            sect = fwy_secs[idx]
            _rel_dist = dist - offset

    # if did not found
    raise Exception('Error: location {0} is out of the network'.format(dist))


# ================================================================ #
# ===================== The I80 BAU deployment =================== #
# ================================================================ #
if False:

    algorithms = ['linearFILL', 'fisb', 'enkfANupdated']

    I80 = OrderedDict()
    # location and types of sensors in the deployed system
    I80['RADAREB4'] = [900.0, 'radar']
    I80['RTMSEB5'] = [1700.0, 'rtms']
    I80['RADAREB6'] = [2300.0, 'radar']
    I80['RTMSEB7'] = [3300.0, 'rtms']
    I80['RADAREB8'] = [3900.0, 'radar']
    I80['RTMSEB9'] = [4900.0, 'rtms']
    I80['RADAREB10'] = [5500.0, 'radar']
    I80['RADAREB11'] = [6300.0, 'radar']
    I80['RTMSEB12'] = [7300.0, 'rtms']
    I80['RADAREB14'] = [8100.0, 'radar']

    # measurement error parameters
    flow_std = {'radar': 0.3, 'rtms': 0.15}
    speed_std = {'radar': 1.5, 'rtms': 5.0}
    num_ensembles = 200

    f.write('config_id:configI80\n')

    for sensor in I80.keys():
        abs_loc = I80[sensor][0]
        sensor_type = I80[sensor][1]
        sec, rel_dist = __get_relative_loc(abs_loc)

        line = 'sensor;id:{0};type:{1};section:{2};distance:{3};flow_std_specs:{4};speed_std_specs:{5};'.format(
                sensor, sensor_type, sec, rel_dist, flow_std[sensor_type], speed_std[sensor_type]) \
               + 'aggregation_sec:30\n'
        f.write(line)

    # each algorithm will run on the sensor grid
    for algorithm in algorithms:
        # write the algorithm line with parameters
        if algorithm == 'constNAN':
            f.write(
                'algorithm;id:constNAN;type:interpolation;interpolation_option:constant;' +
                'queue_threshold:17.88;missing_data:blank\n')
        elif algorithm == 'constNU':
            f.write(
                'algorithm;id:constNU;type:interpolation;interpolation_option:constant;' +
                'queue_threshold:17.88;missing_data:no_update\n')
        elif algorithm == 'constFILL':
            f.write(
                'algorithm;id:constFILL;type:interpolation;interpolation_option:constant;' +
                'queue_threshold:17.88;missing_data:fill\n')
        elif algorithm == 'linearNAN':
            f.write(
                'algorithm;id:linearNAN;type:interpolation;interpolation_option:linear;' +
                'queue_threshold:17.88;missing_data:blank\n')
        elif algorithm == 'linearNU':
            f.write(
                'algorithm;id:linearNU;type:interpolation;interpolation_option:linear;' +
                'queue_threshold:17.88;missing_data:no_update\n')
        elif algorithm == 'linearFILL':
            f.write(
                'algorithm;id:linearFILL;type:interpolation;interpolation_option:linear;' +
                'queue_threshold:17.88;missing_data:fill\n')

        elif algorithm == 'enkfANupdated':

            f.write('algorithm;id:enkfANupdated;type:enkf_AN;queue_threshold:17.88;num_ensembles:{0};'.format(
                num_ensembles) +
                    'vm:27.19;beta:0.6214;rhoc:0.0439;wc:-4.15;' +
                    'std_cell:0.01;std_oncell:0.015;std_offcell:0.012;std_qin:0.1;std_qout:0.2;' +
                    'init_rho:0;init_qin:0.5;init_qout:0.0\n')

        elif algorithm == 'fisb':

            f.write('algorithm;id:fisb;type:fisb;queue_threshold:17.88;' +
                    'omg_c:-14.95;omg_f:97.88;sigma_fac:0.75;tau_fac:0.75;v_crit:90.96;delta_v:20;lag:0;' +
                    'time_diff_cutoff:150\n')

        else:
            raise Exception('unrecognized algorithm')

    f.write('\n')


# ================================================================ #
# =================== Other configurations ======================= #
# ================================================================ #
if True:

    # ================================================
    # the complete set of algorithms to evaluate
    # algorithms = ['linearFILL', 'fisb', 'enkfANupdated']
    algorithms = ['linearFILL', 'fisb']
    # algorithms = ['enkfANupdated']

    # ================================================
    # The complete set of num_sensors to generate
    num_sensors_scenarios = [2, 3, 4, 5, 6, 7, 8, 9, 11, 14, 21, 41]

    # ================================================
    # the complete set of types of sensors to generate
    # sensor_types = ['IDEAL', 'RADAR', 'RADARx2', 'RADARx4', 'RADARx8',
    #                 'RTMS', 'RTMSx2', 'RTMSx4', 'RTMSx8', 'ICONE', 'ICONEx2']
    sensor_types = ['RADAR', 'RADARx2', 'RADARx4', 'RADARx8',
                    'RTMS', 'RTMSx2', 'RTMSx4', 'RTMSx8', 'ICONE', 'ICONEx2']

    # ================================================
    # Generate the configuration file
    # errors, in veh/s and m/s
    sensor_errors = OrderedDict()  # flow_std, speed_std, sensor class
    sensor_errors['IDEAL'] = [0.2, 5.0, 'vol_gt']
    sensor_errors['RADAR'] = [0.3, 5.0, 'radar']
    sensor_errors['RADARx2'] = [0.3, 5.0, 'radar']
    sensor_errors['RADARx4'] = [0.3, 5.0, 'radar']
    sensor_errors['RADARx8'] = [0.3, 5.0, 'radar']
    sensor_errors['RTMS'] = [0.2, 5.0, 'rtms']
    sensor_errors['RTMSx2'] = [0.2, 5.0, 'rtms']
    sensor_errors['RTMSx4'] = [0.2, 5.0, 'rtms']
    sensor_errors['RTMSx8'] = [0.2, 5.0, 'rtms']
    sensor_errors['ICONE'] = [0.3, 3.0, 'radar']
    sensor_errors['ICONEx2'] = [0.3, 3.0, 'radar']

    # get the spacing of each sensor configurations.
    num_cells = int((end_loc - start_loc) / cell_length)
    sensor_locations = OrderedDict()
    for num_sensors in num_sensors_scenarios:
        tmp_spacing = np.linspace(0, num_cells, num_sensors)

        # round the location of each sensor to its closest cell boundary
        sensor_locations[num_sensors] = np.array([round(i) for i in tmp_spacing]) * 200.0 + start_loc

        uneven_spacing = sensor_locations[num_sensors][1:] - sensor_locations[num_sensors][:-1]
        print('{0} sensors, avg {1}, {2}'.format(num_sensors,
                                                 np.mean(uneven_spacing), uneven_spacing))

    # for all different spacing
    for num_sensors in sensor_locations.keys():
        # for all different sensor types
        for sensor_type in sensor_types:
            config_name = 'configHOMO_{0}_{1}'.format(num_sensors, sensor_type)
            f.write('config_id:{0}\n'.format(config_name))

            # icone sensors uses 400 ensembles
            if 'ICONE' in sensor_type:
                num_ensembles = 400
                agg_interval = 60
            else:
                num_ensembles = 200
                agg_interval = 30

            for abs_loc in sensor_locations[num_sensors]:
                sec, rel_dist = __get_relative_loc(abs_loc)

                cus_line = 'sensor;id:{0}Loc{1}m;type:{2};section:{3};distance:{4};'.format(sensor_type,
                                                                                100 * int(round(abs_loc / 100.0)),
                                                                                sensor_errors[sensor_type][2],
                                                                                sec, rel_dist,) \
                           + 'flow_std_specs:{0};speed_std_specs:{1};aggregation_sec:{2}\n'.format(
                                    sensor_errors[sensor_type][0],
                                    sensor_errors[sensor_type][1],
                                    agg_interval)

                f.write(cus_line)

            # each algorithm will run on the sensor grid
            for algorithm in algorithms:
                # write the algorithm
                if algorithm == 'constNAN':
                    f.write(
                        'algorithm;id:constNAN;type:interpolation;interpolation_option:constant;' +
                        'queue_threshold:17.88;missing_data:blank\n')
                elif algorithm == 'constNU':
                    f.write(
                        'algorithm;id:constNU;type:interpolation;interpolation_option:constant;' +
                        'queue_threshold:17.88;missing_data:no_update\n')
                elif algorithm == 'constFILL':
                    f.write(
                        'algorithm;id:constFILL;type:interpolation;interpolation_option:constant;' +
                        'queue_threshold:17.88;missing_data:fill\n')
                elif algorithm == 'linearNAN':
                    f.write(
                        'algorithm;id:linearNAN;type:interpolation;interpolation_option:linear;' +
                        'queue_threshold:17.88;missing_data:blank\n')
                elif algorithm == 'linearNU':
                    f.write(
                        'algorithm;id:linearNU;type:interpolation;interpolation_option:linear;' +
                        'queue_threshold:17.88;missing_data:no_update\n')
                elif algorithm == 'linearFILL':
                    f.write(
                        'algorithm;id:linearFILL;type:interpolation;interpolation_option:linear;' +
                        'queue_threshold:17.88;missing_data:fill\n')

                elif algorithm == 'enkfANupdated':
                    # updated FD using RTMS sensor measurements
                    f.write('algorithm;id:enkfANupdated;type:enkf_AN;queue_threshold:17.88;num_ensembles:{0};'.format(
                        num_ensembles) +
                            'vm:27.19;beta:0.6214;rhoc:0.0439;wc:-4.15;' +
                            'std_cell:0.01;std_oncell:0.015;std_offcell:0.012;std_qin:0.1;std_qout:0.2;' +
                            'init_rho:0;init_qin:0.5;init_qout:0.0\n')

                elif algorithm == 'fisb':
                    f.write('algorithm;id:fisb;type:fisb;queue_threshold:17.88;' +
                            'omg_c:-14.95;omg_f:97.88;sigma_fac:0.75;tau_fac:0.75;v_crit:90.96;delta_v:20;lag:0;' +
                            'time_diff_cutoff:150\n')
                else:
                    raise Exception('unrecognized algorithm')

            f.write('\n')

f.close()

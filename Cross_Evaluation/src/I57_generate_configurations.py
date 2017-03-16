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
# ======================== Generate I57 Configuration ======================================= #
# =========================================================================================== #

# ===================================================
# Configure the I57 network
start_loc = 353.1
end_loc = 6753.1
cell_length = 200

# the order of main freeway sections
fwy_secs = [51413, 51427, 51426]
# the length of sections
len_secs = [6939.5, 13.6, 303.14]

# the offset of freeway sections. that is the absolute distance of the section from the entrance of freeway
offset_secs = np.concatenate([[0], np.cumsum(np.array(len_secs))])

update_file = True

if update_file is True:
    directory = os.path.dirname(os.path.realpath(__file__))
    if sys.platform == 'win32':
        f = open(directory + '\\..\\I57_configurations_input.txt', 'w+')
    elif sys.platform == 'darwin':
        f = open(directory + '/../I57_configurations_input.txt', 'w+')


# ================================================================ #
# ======================= A few utility functions ================ #
def __get_relative_loc(dist):
    """
    This function returns the relative location on the freeway using dist to the entrace
    :param dist: meters, to the entrance
    :return: sect, rel_dist(m)
    """
    sect = 0
    rel_dist = 0.0

    for i, offset in enumerate(offset_secs):

        if dist < offset:
            # found
            return sect, rel_dist
        else:
            # keep updating
            sect = fwy_secs[i]
            rel_dist = dist - offset

    # if did not found
    raise Exception('Error: location {0} is out of the network'.format(dist))


# ================================================================ #
# =============== Generate I57 deployment ============= #
# ================================================================ #
if False:
    I57 = OrderedDict()
    # The location of sensors mapped to the discretized grid.
    I57['RADARSB1'] = [353.1, 'radar']
    I57['RADARSB2'] = [1153.1, 'radar']
    I57['RADARSB3'] = [1953.1, 'radar']
    I57['RADARSB4'] = [2753.1, 'radar']
    I57['RADARSB5'] = [3553.1, 'radar']
    I57['RADARSB6'] = [4353.1, 'radar']
    I57['RTMSSB7'] = [5153.1, 'rtms']
    I57['RADARSB8'] = [5953.1, 'radar']
    I57['RADARSB9'] = [6753.1, 'radar']

    flow_std = {'radar': 0.3, 'rtms': 0.2}
    speed_std = {'radar': 5.0, 'rtms': 5.0}

    f.write('config_id:configI57\n')

    for sensor in I57.keys():

        abs_loc = I57[sensor][0]
        sensor_type = I57[sensor][1]
        if sensor_type == 'radar':
            sensor_model = 'RADAR'
        elif sensor_type == 'rtms':
            sensor_model = 'RTMS'

        sec, rel_dist = __get_relative_loc(abs_loc)
        cus_line = 'sensor;id:{0};type:{1};section:{2};distance:{3};flow_std_specs:{4};speed_std_specs:{5};aggregation_sec:30'.format(
            sensor, sensor_type, sec, rel_dist, flow_std[sensor_type], speed_std[sensor_type])

        # att_list = []
        # for key in sensor_paras[sensor_model].keys():
        #     entry = ':'.join( [key, str( sensor_paras[sensor_model][key])])
        #     att_list.append(entry)
        # att_line = ';'.join(att_list)
        # line = cus_line + att_line + '\n'
        f.write(cus_line + '\n')

    algorithms = ['enkfANupdated']
    num_ensembles = 200
    for algorithm in algorithms:
        # write the algorithm
        if algorithm == 'constNAN':
            f.write(
                'algorithm;id:constNAN;type:interpolation;interpolation_option:constant;queue_threshold:17.88;missing_data:blank\n')
        elif algorithm == 'constNU':
            f.write(
                'algorithm;id:constNU;type:interpolation;interpolation_option:constant;queue_threshold:17.88;missing_data:no_update\n')
        elif algorithm == 'constFILL':
            f.write(
                'algorithm;id:constFILL;type:interpolation;interpolation_option:constant;queue_threshold:17.88;missing_data:fill\n')
        elif algorithm == 'linearNAN':
            f.write(
                'algorithm;id:linearNAN;type:interpolation;interpolation_option:linear;queue_threshold:17.88;missing_data:blank\n')
        elif algorithm == 'linearNU':
            f.write(
                'algorithm;id:linearNU;type:interpolation;interpolation_option:linear;queue_threshold:17.88;missing_data:no_update\n')
        elif algorithm == 'linearFILL':
            f.write(
                'algorithm;id:linearFILL;type:interpolation;interpolation_option:linear;queue_threshold:17.88;missing_data:fill\n')
        elif 'enkf' in algorithm:
            f.write('algorithm;id:enkfANupdated;type:enkf_AN;queue_threshold:17.88;num_ensembles:{0};'.format(
                num_ensembles) +
                    'vm:23.56;beta:0.6214;rhoc:0.0578;wc:-4.88;' +
                    'std_cell:0.01;std_oncell:0.015;std_offcell:0.012;std_qin:0.1;std_qout:0.2;' +
                    'init_rho:0;init_qin:0.5;init_qout:0.0\n')
        elif algorithm == 'fisb':
            f.write('algorithm;id:fisb;type:fisb;queue_threshold:17.88;' +
                    'omg_c:-17.57;omg_f:84.83;sigma_fac:0.75;tau_fac:0.75;v_crit:76.94;delta_v:20;lag:0;' +
                    'time_diff_cutoff:150\n')
        else:
            raise Exception('unrecognized algorithm')
# ================================================================ #


# ================================================================ #
# ================= Exploring all configurations ================= #
# ================================================================ #
if False:

    # ================================================
    # the complete set of algorithms to evaluate
    algorithms = ['linearFILL', 'fisb', 'enkfANupdated']

    # ================================================
    # The complete set of num_sensors to generate
    num_sensors_scenarios = [2, 3, 4, 5, 6, 7, 9, 12, 17, 33]

    # ================================================
    # the complete set of types of sensors to generate
    sensor_types = ['IDEAL', 'RADAR', 'RADARx2', 'RADARx4', 'RADARx8',
                    'RTMS', 'RTMSx2', 'RTMSx4', 'RTMSx8', 'ICONE', 'ICONEx2']

    # errors, in veh/s and m/s
    sensor_errors = OrderedDict()
    sensor_errors['IDEAL'] = [0.2, 5.0, 'vol_gt']  # changed the variance
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

    # get the location of each sensor
    num_cells = int((end_loc - start_loc) / cell_length)
    sensor_locations = OrderedDict()
    for num_sensors in num_sensors_scenarios:
        tmp_spacing = np.linspace(0, num_cells, num_sensors)

        # round the location of each sensor to its closest cell boundary
        sensor_locations[num_sensors] = np.array([round(i) for i in tmp_spacing]) * 200.0 + start_loc

        uneven_spacing = sensor_locations[num_sensors][1:] - sensor_locations[num_sensors][:-1]
        print('{0} sensors, avg {1}, {2}'.format(num_sensors,
                                                 np.mean(uneven_spacing), uneven_spacing))

    # write the file
    if update_file is True:
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

                        f.write(
                            'algorithm;id:enkfANupdated;type:enkf_AN;queue_threshold:17.88;num_ensembles:{0};'.format(
                                num_ensembles) +
                            'vm:23.56;beta:0.6214;rhoc:0.0578;wc:-4.88;' +
                            'std_cell:0.01;std_oncell:0.015;std_offcell:0.012;std_qin:0.1;std_qout:0.2;' +
                            'init_rho:0;init_qin:0.5;init_qout:0.0\n')

                    elif algorithm == 'fisb':
                        f.write('algorithm;id:fisb;type:fisb;queue_threshold:17.88;' +
                                'omg_c:-17.57;omg_f:84.83;sigma_fac:0.75;tau_fac:0.75;v_crit:76.94;delta_v:20;lag:0;' +
                                'time_diff_cutoff:150\n')
                    else:
                        raise Exception('unrecognized algorithm')

                f.write('\n')


f.close()

from collections import OrderedDict

import numpy as np

from cross_evaluation import CrossEval
from sensor_paras import sensor_paras

__author__ = 'Yanning Li and Juan Carlos Martinez'
"""
This script generates the virtual sensor data.
"""

# =========================================================================================== #
# ================================== Configurations ========================================= #
# =========================================================================================== #
# set the folders and directories
workzone = 'I57'
log_dir = '../'
data_dir = '../'

# --------------------------------
# configure to select which road network configuration
small_net = False  # Small net is only used for generating virtual sensor data for FD calibration.
if small_net is True:
    fwy_secs = [701, 457, 456]
    len_secs = [1629.71, 2.44, 111.54]
    print('===========================================================================')
    print('============================== WARNING ====================================')
    print('===========================================================================')
    print('This is not I57 network. This is a small network just for FD calibration!')
    print('===========================================================================')
    print('===========================================================================')
    print('===========================================================================')
else:
    # the order of main freeway sections
    fwy_secs = [51413, 51427, 51426]
    # the length of sections
    len_secs = [6939.5, 13.6, 303.14]

# the offset of freeway sections. that is the absolute distance of the section from the entrance of freeway
offset_secs = np.concatenate([[0], np.cumsum(np.array(len_secs))])

# --------------------------------
# configure the replication
rep = 51642

# --------------------------------
# configure the sensors sensing location
start_loc = 353.1
end_loc = 6753.1
cell_length = 200
distance = np.arange(start_loc, end_loc, cell_length)
distance = np.concatenate( [distance, np.array([end_loc]) ] )
# round to 2 decimals
distance = [ round(i, 2) for i in distance ]

# for small net
# distance = [1400, 1600]

# --------------------------------
# configure sensors
# Here is the complete set of sensors to generate in I57 work zone (in total 13)
sensors_to_generate = ['IDEAL',
                       'RTMS', 'RTMSx2', 'RTMSx4', 'RTMSx8',
                       'RADAR', 'RADARx2', 'RADARx4', 'RADARx8',
                       'ICONE', 'ICONEx2']

update_file = True
generate_data = True

if update_file is True:
    f = open('../Virtual_sensor_data/{0}_rep{1}_to_generate.txt'.format(workzone, rep), 'w+')
    # f.write('\n')

# =========================================================================================== #
# ============================ Finished configuration ======================================= #
# =========================================================================================== #


# =====================================================
# A utility function
# =====================================================
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


# =========================================================================================== #
# ================================== BAU I57 deployment ===================================== #
# =========================================================================================== #
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

    flow_std = {'radar': 0.3, 'rtms': 0.15}
    speed_std = {'radar': 1.5, 'rtms': 5.0}

    for sensor in I57.keys():

        abs_loc = I57[sensor][0]
        sensor_type = I57[sensor][1]
        if sensor_type == 'radar':
            sensor_model = 'RADAR'
        elif sensor_type == 'rtms':
            sensor_model = 'RTMS'

        sec, rel_dist = __get_relative_loc(abs_loc)
        cus_line = 'id:{0};type:{1};section:{2};distance:{3};flow_std_specs:{4};speed_std_specs:{5};'.format(
            sensor,
            sensor_type, sec, rel_dist, flow_std[sensor_type], speed_std[sensor_type])

        att_list = []
        for key in sensor_paras[sensor_model].keys():
            entry = ':'.join([key, str(sensor_paras[sensor_model][key])])
            att_list.append(entry)
        att_line = ';'.join(att_list)
        line = cus_line + att_line + '\n'
        f.write(line)

# =========================================================================================== #
# ================================== IDEAL sensor =========================================== #
# =========================================================================================== #
# # generate IDEAL sensors
# # NOTE: flow_std_specs and speed_std_specs are in veh/s, and m/s, all others are in meters, and kph
if 'IDEAL' in sensors_to_generate and update_file is True:
    for i in range(0, len(distance)):

        sec, rel_dist = __get_relative_loc(distance[i])
        cus_line = 'id:IDEALLoc{0}m;type:vol_gt;section:{1};distance:{2};flow_std_specs:0.012;speed_std_specs:0.05;'.format(
            100 * int(round((distance[i]) / 100.0)), sec, rel_dist)

        att_list = []
        for key in sensor_paras['IDEAL'].keys():
            entry = ':'.join([key, str(sensor_paras['IDEAL'][key])])
            att_list.append(entry)
        att_line = ';'.join(att_list)
        line = cus_line + att_line + '\n'
        f.write(line)

# =========================================================================================== #
# ================================== RADAR sensor =========================================== #
# =========================================================================================== #
# generate RADAR sensor
# NOTE: flow_std_specs and speed_std_specs are in veh/s, and m/s, all others are in meters, and kph
if 'RADAR' in sensors_to_generate and update_file is True:
    for i in range(0, len(distance)):

        sec, rel_dist = __get_relative_loc(distance[i])
        cus_line = 'id:RADARLoc{0}m;type:radar;section:{1};distance:{2};flow_std_specs:0.3;speed_std_specs:1.5;'.format(
            100 * int(round((distance[i]) / 100.0)), sec, rel_dist)

        att_list = []
        for key in sensor_paras['RADAR'].keys():
            entry = ':'.join([key, str(sensor_paras['RADAR'][key])])
            att_list.append(entry)
        att_line = ';'.join(att_list)
        line = cus_line + att_line + '\n'
        f.write(line)

# ===============================
# generate RADARx2 sensor
# NOTE: flow_std_specs and speed_std_specs are in veh/s, and m/s, all others are in meters, and kph
if 'RADARx2' in sensors_to_generate and update_file is True:
    for i in range(0, len(distance)):

        sec, rel_dist = __get_relative_loc(distance[i])
        cus_line = 'id:RADARx2Loc{0}m;type:radar;section:{1};distance:{2};flow_std_specs:0.3;speed_std_specs:0.75;'.format(
            100 * int(round((distance[i]) / 100.0)), sec, rel_dist)

        att_list = []
        for key in sensor_paras['RADARx2'].keys():
            entry = ':'.join([key, str(sensor_paras['RADARx2'][key])])
            att_list.append(entry)
        att_line = ';'.join(att_list)
        line = cus_line + att_line + '\n'
        f.write(line)

# ===============================
# generate RADARx4 sensor
# NOTE: flow_std_specs and speed_std_specs are in veh/s, and m/s, all others are in meters, and kph
if 'RADARx4' in sensors_to_generate and update_file is True:
    for i in range(0, len(distance)):

        sec, rel_dist = __get_relative_loc(distance[i])
        cus_line = 'id:RADARx4Loc{0}m;type:radar;section:{1};distance:{2};flow_std_specs:0.3;speed_std_specs:0.5;'.format(
            100 * int(round((distance[i]) / 100.0)), sec, rel_dist)

        att_list = []
        for key in sensor_paras['RADARx4'].keys():
            entry = ':'.join([key, str(sensor_paras['RADARx4'][key])])
            att_list.append(entry)
        att_line = ';'.join(att_list)
        line = cus_line + att_line + '\n'
        f.write(line)

# ===============================
# generate RADARx8 sensor
# NOTE: flow_std_specs and speed_std_specs are in veh/s, and m/s, all others are in meters, and kph
if 'RADARx8' in sensors_to_generate and update_file is True:
    for i in range(0, len(distance)):

        sec, rel_dist = __get_relative_loc(distance[i])
        cus_line = 'id:RADARx8Loc{0}m;type:radar;section:{1};distance:{2};flow_std_specs:0.3;speed_std_specs:0.25;'.format(
            100 * int(round((distance[i]) / 100.0)), sec, rel_dist)

        att_list = []
        for key in sensor_paras['RADARx8'].keys():
            entry = ':'.join([key, str(sensor_paras['RADARx8'][key])])
            att_list.append(entry)
        att_line = ';'.join(att_list)
        line = cus_line + att_line + '\n'
        f.write(line)


# =========================================================================================== #
# =================================== RTMS sensor =========================================== #
# =========================================================================================== #
# RTMS I57
# NOTE: flow_std_specs and speed_std_specs are in veh/s, and m/s, all others are in meters, and kph
if 'RTMSI57' in sensors_to_generate and update_file is True:
    for i in range(0, len(distance)):

        sec, rel_dist = __get_relative_loc(distance[i])
        cus_line = 'id:RTMSI57Loc{0}m;type:rtms;section:{1};distance:{2};flow_std_specs:0.15;speed_std_specs:8.0;'.format(
            100 * int(round((distance[i]) / 100.0)), sec, rel_dist)

        att_list = []
        for key in sensor_paras['RTMSI57'].keys():
            entry = ':'.join([key, str(sensor_paras['RTMSI57'][key])])
            att_list.append(entry)
        att_line = ';'.join(att_list)
        line = cus_line + att_line + '\n'
        f.write(line)

# ===============================
# RTMS
# NOTE: flow_std_specs and speed_std_specs are in veh/s, and m/s, all others are in meters, and kph
if 'RTMS' in sensors_to_generate and update_file is True:
    for i in range(0, len(distance)):

        sec, rel_dist = __get_relative_loc(distance[i])
        cus_line = 'id:RTMSLoc{0}m;type:rtms;section:{1};distance:{2};flow_std_specs:0.15;speed_std_specs:5.0;'.format(
            100 * int(round((distance[i]) / 100.0)), sec, rel_dist)

        att_list = []
        for key in sensor_paras['RTMS'].keys():
            entry = ':'.join([key, str(sensor_paras['RTMS'][key])])
            att_list.append(entry)
        att_line = ';'.join(att_list)
        line = cus_line + att_line + '\n'
        f.write(line)

# ===============================
# RTMSx2
# NOTE: flow_std_specs and speed_std_specs are in veh/s, and m/s, all others are in meters, and kph
if 'RTMSx2' in sensors_to_generate and update_file is True:
    for i in range(0, len(distance)):

        sec, rel_dist = __get_relative_loc(distance[i])
        cus_line = 'id:RTMSx2Loc{0}m;type:rtms;section:{1};distance:{2};flow_std_specs:0.075;speed_std_specs:2.5;'.format(
            100 * int(round((distance[i]) / 100.0)), sec, rel_dist)

        att_list = []
        for key in sensor_paras['RTMSx2'].keys():
            entry = ':'.join([key, str(sensor_paras['RTMSx2'][key])])
            att_list.append(entry)
        att_line = ';'.join(att_list)
        line = cus_line + att_line + '\n'
        f.write(line)

# ===============================
# RTMSx4
# NOTE: flow_std_specs and speed_std_specs are in veh/s, and m/s, all others are in meters, and kph
if 'RTMSx4' in sensors_to_generate and update_file is True:
    for i in range(0, len(distance)):

        sec, rel_dist = __get_relative_loc(distance[i])
        cus_line = 'id:RTMSx4Loc{0}m;type:rtms;section:{1};distance:{2};flow_std_specs:0.0375;speed_std_specs:1.25;'.format(
            100 * int(round((distance[i]) / 100.0)), sec, rel_dist)

        att_list = []
        for key in sensor_paras['RTMSx4'].keys():
            entry = ':'.join([key, str(sensor_paras['RTMSx4'][key])])
            att_list.append(entry)
        att_line = ';'.join(att_list)
        line = cus_line + att_line + '\n'
        f.write(line)

# ===============================
# RTMSx8
# NOTE: flow_std_specs and speed_std_specs are in veh/s, and m/s, all others are in meters, and kph
if 'RTMSx8' in sensors_to_generate and update_file is True:
    for i in range(0, len(distance)):

        sec, rel_dist = __get_relative_loc(distance[i])
        cus_line = 'id:RTMSx8Loc{0}m;type:rtms;section:{1};distance:{2};flow_std_specs:0.02;speed_std_specs:0.6;'.format(
            100 * int(round((distance[i]) / 100.0)), sec, rel_dist)

        att_list = []
        for key in sensor_paras['RTMSx8'].keys():
            entry = ':'.join([key, str(sensor_paras['RTMSx8'][key])])
            att_list.append(entry)
        att_line = ';'.join(att_list)
        line = cus_line + att_line + '\n'
        f.write(line)

# =========================================================================================== #
# ==================================== ICONE sensor ========================================= #
# =========================================================================================== #
# ICONE
# NOTE: flow_std_specs and speed_std_specs are in veh/s, and m/s, all others are in meters, and kph
if 'ICONE' in sensors_to_generate and update_file is True:
    for i in range(0, len(distance)):

        sec, rel_dist = __get_relative_loc(distance[i])
        cus_line = 'id:ICONELoc{0}m;type:icone;section:{1};distance:{2};flow_std_specs:0.3;speed_std_specs:5.0;'.format(
            100 * int(round((distance[i]) / 100.0)), sec, rel_dist)

        att_list = []
        for key in sensor_paras['ICONE'].keys():
            entry = ':'.join([key, str(sensor_paras['ICONE'][key])])
            att_list.append(entry)
        att_line = ';'.join(att_list)
        line = cus_line + att_line + '\n'
        f.write(line)

# ===============================
# ICONEx2
# NOTE: flow_std_specs and speed_std_specs are in veh/s, and m/s, all others are in meters, and kph
if 'ICONEx2' in sensors_to_generate and update_file is True:
    for i in range(0, len(distance)):

        sec, rel_dist = __get_relative_loc(distance[i])
        cus_line = 'id:ICONEx2Loc{0}m;type:icone;section:{1};distance:{2};flow_std_specs:0.3;speed_std_specs:2.5;'.format(
            100 * int(round((distance[i]) / 100.0)), sec, rel_dist)

        att_list = []
        for key in sensor_paras['ICONEx2'].keys():
            entry = ':'.join([key, str(sensor_paras['ICONEx2'][key])])
            att_list.append(entry)
        att_line = ';'.join(att_list)
        line = cus_line + att_line + '\n'
        f.write(line)


if update_file is True:
    f.close()

# =========================================================================================== #
# ================================== Generate sensors ======================================= #
# =========================================================================================== #
# #
if generate_data is True:
    cross_eval = CrossEval(workzone=workzone, log_dir=log_dir, data_dir=data_dir,
                            grid_res=(5, 50), replications=[rep])

    # fetch all data needed for the evaluation
    cross_eval.just_fetch()

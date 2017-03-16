import numpy as np
from cross_evaluation import CrossEval
from sensor_paras import sensor_paras
from collections import OrderedDict

__author__ = 'Yanning Li and Juan Carlos Martinez'
"""
This script generates the virtual sensor data.
"""

# =========================================================================================== #
# ================================== Configurations I-80 ==================================== #
# =========================================================================================== #
# set the folders and directories
workzone = 'I80'
log_dir = '../'
data_dir = '../'

# --------------------------------
# configure to select which road network configuration
small_net = False   # Small net is only used for generating virtual sensor data for FD calibration.
if small_net is True:
    fwy_secs = [701, 457, 456]
    len_secs = [1629.71, 2.44, 111.54]

else:
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

# --------------------------------
# configure the replication
rep = 41368

# --------------------------------
# configure the sensors sensing locations
start_loc = 300
end_loc = 8300
cell_length = 200
distance = np.arange(start_loc, end_loc, cell_length)
distance = np.concatenate( [distance, np.array([end_loc]) ] )
distance = [round(i, 2) for i in distance]  # round to 2 decimals


# --------------------------------
# configure sensors
sensors_to_generate = ['IDEAL', 'RTMS', 'RTMSx2', 'RTMSx4', 'RTMSx8',
                       'RADAR', 'RADARx2', 'RADARx4', 'RADARx8',
                       'ICONE', 'ICONEx2']

update_file = True
generate_data = True

if update_file is True:
    f = open('../Virtual_sensor_data/{0}_rep{1}_to_generate.txt'.format(workzone, rep), 'w+')

# =========================================================================================== #
# ============================ Finished configuration ======================================= #
# =========================================================================================== #


# =====================================================
# A utility function
# =====================================================
def __get_relative_loc(dist):
    """
    This function returns the relative location on the freeway using dist to the entrance
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


# =========================================================================================== #
# ================================== I80 sensor ============================================= #
# =========================================================================================== #
# shift I80 sensors a bit to the closest space grid.
# sensor locations, I80[sensor_id] = [abs_loc, sensor type]
# The physical location of sensors are:
# """
# RADAREB4 abs loc:947.87
# RTMSEB5 abs loc:1606.47
# RADAREB6 abs loc:2310.07
# RTMSEB7 abs loc:3218.72
# RADAREB8 abs loc:4048.43
# RTMSEB9 abs loc:4878.03
# RADAREB10 abs loc:5616.03
# RADAREB11 abs loc:6323.43
# RTMSEB12 abs loc:7241.84
# RADAREB14 abs loc:8110.55
# """
# The locations the sensor measure are:
# """
# RADAREB4 abs loc:873.17
# RTMSEB5 abs loc:1606.47
# RADAREB6 abs loc:2235.37
# RTMSEB7 abs loc:3218.72
# RADAREB8 abs loc:3973.73
# RTMSEB9 abs loc:4878.03
# RADAREB10 abs loc:5541.33
# RADAREB11 abs loc:6248.73
# RTMSEB12 abs loc:7241.84
# RADAREB14 abs loc:8035.85
# """

if False:
    I80 = OrderedDict()
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

    flow_std = {'radar':0.3, 'rtms':0.15}
    speed_std = {'radar':1.5, 'rtms':5.0}

    for sensor in I80.keys():

        abs_loc = I80[sensor][0]
        sensor_type = I80[sensor][1]
        if sensor_type == 'radar':
            sensor_model = 'RADAR'
        elif sensor_type == 'rtms':
            sensor_model = 'RTMS'

        sec, rel_dist = __get_relative_loc(abs_loc)
        cus_line = 'id:{0};type:{1};section:{2};distance:{3};flow_std_specs:{4};speed_std_specs:{5};'.format(sensor,
            sensor_type, sec, rel_dist, flow_std[sensor_type], speed_std[sensor_type])

        att_list = []
        for key in sensor_paras[sensor_model].keys():
            entry = ':'.join( [key, str( sensor_paras[sensor_model][key] )] )
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
        cus_line = 'id:IDEALLoc{0}m;type:vol_gt;section:{1};distance:{2};flow_std_specs:0.012;speed_std_specs:0.05;'.\
            format(100 * int(round((distance[i]) / 100.0)), sec, rel_dist)

        att_list = []
        for key in sensor_paras['IDEAL'].keys():
            entry = ':'.join([key, str(sensor_paras['IDEAL'][key])])
            att_list.append(entry)
        att_line = ';'.join(att_list)
        line = cus_line + att_line + '\n'
        f.write(line)
# ===================== END SECTION =====================

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
        cus_line = 'id:RADARx2Loc{0}m;type:radar;section:{1};distance:{2};flow_std_specs:0.3;speed_std_specs:0.75;'.\
            format(100 * int(round((distance[i]) / 100.0)), sec, rel_dist)

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
        cus_line = 'id:RADARx4Loc{0}m;type:radar;section:{1};distance:{2};flow_std_specs:0.3;speed_std_specs:0.5;'.\
            format(100 * int(round((distance[i]) / 100.0)), sec, rel_dist)

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
        cus_line = 'id:RADARx8Loc{0}m;type:radar;section:{1};distance:{2};flow_std_specs:0.3;speed_std_specs:0.25;'.\
            format(100 * int(round((distance[i]) / 100.0)), sec, rel_dist)

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
        cus_line = 'id:RTMSx2Loc{0}m;type:rtms;section:{1};distance:{2};flow_std_specs:0.075;speed_std_specs:2.5;'.\
            format(100 * int(round((distance[i]) / 100.0)), sec, rel_dist)

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
        cus_line = 'id:RTMSx4Loc{0}m;type:rtms;section:{1};distance:{2};flow_std_specs:0.0375;speed_std_specs:1.25;'.\
            format(100 * int(round((distance[i]) / 100.0)), sec, rel_dist)

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
        cus_line = 'id:RTMSx8Loc{0}m;type:rtms;section:{1};distance:{2};flow_std_specs:0.02;speed_std_specs:0.6;'.\
            format(100 * int(round((distance[i]) / 100.0)), sec, rel_dist)

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
        cus_line = 'id:ICONEx2Loc{0}m;type:icone;section:{1};distance:{2};flow_std_specs:0.3;speed_std_specs:2.5;'.\
            format(100 * int(round((distance[i]) / 100.0)), sec, rel_dist)

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
# Generate the virual sensor data
if generate_data is True:
    cross_eval = CrossEval(workzone=workzone, log_dir=log_dir, data_dir=data_dir,
                            grid_res=(5, 200), replications=[rep])

    # fetch all data needed for the evaluation
    cross_eval.just_fetch()

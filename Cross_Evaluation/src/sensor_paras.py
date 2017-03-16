__author__ = 'Yanning Li'
"""
This file contains the parameters for sensors with different levels of accuracy. Use this when generating the configuration file.
Particularly, it contains the following types of sensors:

    - IDEAL: the ideal volume and flow sensor

    - RADARI57: the radar sensor with similar performance as observed in I57 work zone
    - RADAR: the radar with normal performance, e.g. in I80., assuming well calibrated to remove bias for single measurement
    - RADARx2: the radar sensor with 2 times improvement in performance: 2 times less Gaussian noise, and data missing
    - RADARxn: n times improvement

    - RTMSI57: the RTMS sensor with similar performance as observed in I57 and I80 work zone.
    - RTMS: the RTMS with normal performance, e.g. in I80, assuming well calibrated to remove bias for single measurement
    - RTMSx2: the RTMS sensor with 2 times improvement in performance: 2 times less Gaussian noise, and data missing
    - RTMSxn: n times imporovement


    - ICONE: the iCone sensor with normal performance as observed in other evaluation tests,
    - ICONEx2: two time performance

    - BT: the Bluetooth sensor with normal performance as observed in other evaluation tests
    - IDEALBT: the ideal BT sensor

Parameters:
    -- Volume Sensor Configuration -- #
    The following comments specify and explain the default parameters that need
    to be set for construction a 'volume' sensor

    # Deployment Parameters ##
    alpha                      : Angle (CW) between the sensor beam and the road cross-section [deg]
    theta                      : Angle (CW) of visibility on each side of sensor beam [deg]
    s                          : Road shoulder [m]
    L                          : Lane width [m]
    offset_dist                : offset distance to the upstream [m],
                                 e.g., sensor at x will measure traffic at x-offset_dist

    # General Noise Parameters ##
    noise_type                 : Noise type ['relative' or 'absolute']
    occlusion                  : Occlusion [True or False]
    p_occlusion_accept         : Percent range of a vehicle's crossing time that may overlap with another vehicle's
                                 crossing time without causing occlusion
    aggregation_sec            : Aggregation period [s]
    awake_sec                  : Detection phase [s]
    v_range                    : Speed range [kph]
    v_threshold                : Speed threshhold for free flow vs. congested flow [kph]
    p_missing_ff               : Percent missing range in free flow []
    p_missing_cf               : Percent missing range in congested flow []

    # Relative Noise Parameters ##
    v_bias_ff                  : Speed bias range in free flow [kph]
    v_bias_cf                  : Speed bias range in congested flow [kph]
    v_accuracy_p_ff              : Speed accuracy percent range in free flow [0, 100]
    v_accuracy_p_cf              : Speed accuracy percent range in congested flow [0, 100]

    # Absolute Noise Parameters ##
    v_noise_sigma_ff           : Speed noise std. range in free flow [kph]
    v_noise_sigma_cf           : Speed noise std. range in congested flow [kph]
    v_noise_mu_ff              : Speed noise mean range in free flow [kph]
    v_noise_mu_cf              : Speed noise mean range in congested flow [kph]

    # Added noise for counts, see the sensor model write up for details ##
    This is only for RTMS sensor
    c_bias_p_ff                  : count bias in free flow, [0, 100] %
    c_sigma_p_ff                 : count standard deviation in free flow, [0, 100] %
    c_bias_p_cf                  : count bias in free flow, [0, 100] %
    c_sigma_p_cf                 : count standard deviation in free flow, [0, 100] %



    Travel Time Sensor Configuration #
    The following comments specify and explain the default parameters that need
    to be set for construction a 'travel time' sensor

    # Deployment Parameters ##
    alpha                      : Angle (CW) between the sensor beam and the road cross-section [deg] (Set to 0 [Zero])
    theta                      : Angle (CW) of visibility on each side of sensor beam [deg] (Set to 0 [Zero])
    s                          : Road shoulder [m]
    L                          : Lane width [m]

    # General Noise Parameters ##
    aggregation_sec            : Aggregation period [s]
    veh_length                 : Average vehicle length [m]
    p_penetration              : Penetration percentage range [0, 100]
    t_noise_sigma              : Time recording noise sigma range [s]
    t_noise_mu                 : Time recording noise mean range [s]



"""

from collections import OrderedDict
sensor_paras = OrderedDict()

sensor_paras['veh_length'] = 4.5
sensor_paras['time_step'] = 0.2


# =========================================================================================== #
# ================================== IDEAL sensor =========================================== #
# =========================================================================================== #
sensor_paras['IDEAL'] = OrderedDict()

sensor_paras['IDEAL']['alpha'] = 0
sensor_paras['IDEAL']['theta'] = 0
sensor_paras['IDEAL']['s'] = 0.91   # 3 ft
sensor_paras['IDEAL']['L'] = 3.66   # 12 ft
sensor_paras['IDEAL']['offset_dist'] = 0    # m


sensor_paras['IDEAL']['v_range'] = [0,170]  #kph
sensor_paras['IDEAL']['aggregation_sec'] = 30
sensor_paras['IDEAL']['awake_sec'] = 30

sensor_paras['IDEAL']['occlusion'] = 'False'
sensor_paras['IDEAL']['p_occlusion_accept'] = [100,100]
sensor_paras['IDEAL']['v_threshold'] = 0
sensor_paras['IDEAL']['p_missing_ff'] = [0.0,0.0]
sensor_paras['IDEAL']['p_missing_cf'] = [0.0,0.0]


sensor_paras['IDEAL']['noise_type'] = 'absolute'
# relative
sensor_paras['IDEAL']['v_bias_ff']  = [0,0]
sensor_paras['IDEAL']['v_bias_cf'] = [0,0]
sensor_paras['IDEAL']['v_accuracy_p_ff'] = [0,0]
sensor_paras['IDEAL']['v_accuracy_p_cf'] = [0,0]
# absolute
sensor_paras['IDEAL']['v_noise_sigma_ff'] = [0,0]
sensor_paras['IDEAL']['v_noise_sigma_cf'] = [0,0]
sensor_paras['IDEAL']['v_noise_mu_ff'] = [0,0]
sensor_paras['IDEAL']['v_noise_mu_cf'] = [0,0]
# count relative error
sensor_paras['IDEAL']['c_bias_p_ff'] = [0, 0]
sensor_paras['IDEAL']['c_bias_p_cf'] = [0, 0]
sensor_paras['IDEAL']['c_sigma_p_ff'] = [0, 0]
sensor_paras['IDEAL']['c_sigma_p_cf'] = [0, 0]
# --------------------------------- #


# =========================================================================================== #
# ================================== RADAR sensor =========================================== #
# =========================================================================================== #
# ------------- RADARI57 ------------- #
sensor_paras['RADARI57'] = OrderedDict()

sensor_paras['RADARI57']['alpha'] = 84
sensor_paras['RADARI57']['theta'] = 2
sensor_paras['RADARI57']['s']= 0.91
sensor_paras['RADARI57']['L'] = 3.66    # m
sensor_paras['RADARI57']['offset_dist'] = 43.5   # m

sensor_paras['RADARI57']['v_range'] = [64.37, 159.33]   # kph [40, 99] mph
sensor_paras['RADARI57']['aggregation_sec'] = 30
sensor_paras['RADARI57']['awake_sec'] = 30

sensor_paras['RADARI57']['occlusion'] = 'True'
sensor_paras['RADARI57']['p_occlusion_accept'] = [30,30]
sensor_paras['RADARI57']['v_threshold'] = 64.37    # kph   40 mph
sensor_paras['RADARI57']['p_missing_ff'] = [0.0,0.0]
sensor_paras['RADARI57']['p_missing_cf'] = [100.0,100.0]

sensor_paras['RADARI57']['noise_type'] = 'absolute'
# relative, NOT used
sensor_paras['RADARI57']['v_bias_ff']  = [-1,-3]
sensor_paras['RADARI57']['v_bias_cf'] = [-3,-5]
sensor_paras['RADARI57']['v_accuracy_p_ff'] = [5,10]
sensor_paras['RADARI57']['v_accuracy_p_cf'] = [10,20]
# use absolute error
sensor_paras['RADARI57']['v_noise_mu_ff'] = [-4.83, -1.61]     # kph [-3, -1] mph
sensor_paras['RADARI57']['v_noise_mu_cf'] = [-4.83, -1.61]     # kph Not used
sensor_paras['RADARI57']['v_noise_sigma_ff'] = [0.8, 1.6]   # kph [0.5, 1] mph
sensor_paras['RADARI57']['v_noise_sigma_cf'] = [4.8, 8.0]    # kph [3. 5] mph
# count relative error
sensor_paras['RADARI57']['c_bias_p_ff'] = [0, 0]
sensor_paras['RADARI57']['c_bias_p_cf'] = [0, 0]
sensor_paras['RADARI57']['c_sigma_p_ff'] = [0, 0]
sensor_paras['RADARI57']['c_sigma_p_cf'] = [0, 0]
# --------------------------------- #

# ------------- RADAR ------------- #
sensor_paras['RADAR'] = OrderedDict()

sensor_paras['RADAR']['alpha'] = 84
sensor_paras['RADAR']['theta'] = 2
sensor_paras['RADAR']['s']= 0.91
sensor_paras['RADAR']['L'] = 3.66    # m
sensor_paras['RADAR']['offset_dist'] = 43.5    # m

sensor_paras['RADAR']['v_range'] = [8.05, 159.33]   # kph [5, 99] mph
sensor_paras['RADAR']['aggregation_sec'] = 30
sensor_paras['RADAR']['awake_sec'] = 30

sensor_paras['RADAR']['occlusion'] = 'True'
sensor_paras['RADAR']['p_occlusion_accept'] = [30,30]
sensor_paras['RADAR']['v_threshold'] = 64.37    # kph   40 mph
sensor_paras['RADAR']['p_missing_ff'] = [0.0,0.0]
sensor_paras['RADAR']['p_missing_cf'] = [13.0, 17.0]

sensor_paras['RADAR']['noise_type'] = 'absolute'
# relative, NOT used
sensor_paras['RADAR']['v_bias_ff']  = [0, 0]
sensor_paras['RADAR']['v_bias_cf'] = [0, 0]
sensor_paras['RADAR']['v_accuracy_p_ff'] = [0, 0]
sensor_paras['RADAR']['v_accuracy_p_cf'] = [0, 0]
# use absolute error
sensor_paras['RADAR']['v_noise_mu_ff'] = [0, 0]     # kph [-3, -1] mph
sensor_paras['RADAR']['v_noise_mu_cf'] = [0, 0]     # kph Not used
sensor_paras['RADAR']['v_noise_sigma_ff'] = [0.8, 1.6]   # kph [0.5, 1] mph
sensor_paras['RADAR']['v_noise_sigma_cf'] = [4.8, 8.0]    # kph [3, 5] mph
# count relative error
sensor_paras['RADAR']['c_bias_p_ff'] = [0, 0]
sensor_paras['RADAR']['c_bias_p_cf'] = [0, 0]
sensor_paras['RADAR']['c_sigma_p_ff'] = [0, 0]
sensor_paras['RADAR']['c_sigma_p_cf'] = [0, 0]
# --------------------------------- #

# ------------- RADARx2 ------------- #
sensor_paras['RADARx2'] = OrderedDict()

sensor_paras['RADARx2']['alpha'] = 84
sensor_paras['RADARx2']['theta'] = 2
sensor_paras['RADARx2']['s']= 0.91
sensor_paras['RADARx2']['L'] = 3.66    # m
sensor_paras['RADARx2']['offset_dist'] = 43.5    # m

sensor_paras['RADARx2']['v_range'] = [8.05, 159.33]   # kph [5, 99] mph
sensor_paras['RADARx2']['aggregation_sec'] = 30
sensor_paras['RADARx2']['awake_sec'] = 30

sensor_paras['RADARx2']['occlusion'] = 'True'
sensor_paras['RADARx2']['p_occlusion_accept'] = [30,30]
sensor_paras['RADARx2']['v_threshold'] = 64.37    # kph   40 mph
sensor_paras['RADARx2']['p_missing_ff'] = [0.0,0.0]
sensor_paras['RADARx2']['p_missing_cf'] = [6.5, 8.5]

sensor_paras['RADARx2']['noise_type'] = 'absolute'
# relative, NOT used
sensor_paras['RADARx2']['v_bias_ff']  = [0, 0]
sensor_paras['RADARx2']['v_bias_cf'] = [0, 0]
sensor_paras['RADARx2']['v_accuracy_p_ff'] = [0, 0]
sensor_paras['RADARx2']['v_accuracy_p_cf'] = [0, 0]
# use absolute error
sensor_paras['RADARx2']['v_noise_mu_ff'] = [0, 0]     # kph
sensor_paras['RADARx2']['v_noise_mu_cf'] = [0, 0]     # kph Not used
sensor_paras['RADARx2']['v_noise_sigma_ff'] = [0.4, 0.8]   # kph
sensor_paras['RADARx2']['v_noise_sigma_cf'] = [2.4, 4.0]    # kph
# count relative error
sensor_paras['RADARx2']['c_bias_p_ff'] = [0, 0]
sensor_paras['RADARx2']['c_bias_p_cf'] = [0, 0]
sensor_paras['RADARx2']['c_sigma_p_ff'] = [0, 0]
sensor_paras['RADARx2']['c_sigma_p_cf'] = [0, 0]
# --------------------------------- #


# ------------- RADARx4 ------------- #
sensor_paras['RADARx4'] = OrderedDict()

sensor_paras['RADARx4']['alpha'] = 84
sensor_paras['RADARx4']['theta'] = 2
sensor_paras['RADARx4']['s']= 0.91
sensor_paras['RADARx4']['L'] = 3.66    # m
sensor_paras['RADARx4']['offset_dist'] = 43.5    # m

sensor_paras['RADARx4']['v_range'] = [8.05, 159.33]   # kph [5, 99] mph
sensor_paras['RADARx4']['aggregation_sec'] = 30
sensor_paras['RADARx4']['awake_sec'] = 30

sensor_paras['RADARx4']['occlusion'] = 'True'
sensor_paras['RADARx4']['p_occlusion_accept'] = [30,30]
sensor_paras['RADARx4']['v_threshold'] = 64.37    # kph   40 mph
sensor_paras['RADARx4']['p_missing_ff'] = [0.0,0.0]
sensor_paras['RADARx4']['p_missing_cf'] = [3.0, 4.0]

sensor_paras['RADARx4']['noise_type'] = 'absolute'
# relative, NOT used
sensor_paras['RADARx4']['v_bias_ff']  = [0, 0]
sensor_paras['RADARx4']['v_bias_cf'] = [0, 0]
sensor_paras['RADARx4']['v_accuracy_p_ff'] = [0, 0]
sensor_paras['RADARx4']['v_accuracy_p_cf'] = [0, 0]
# use absolute error
sensor_paras['RADARx4']['v_noise_mu_ff'] = [0, 0]     # kph
sensor_paras['RADARx4']['v_noise_mu_cf'] = [0, 0]     # kph Not used
sensor_paras['RADARx4']['v_noise_sigma_ff'] = [0.2, 0.4]   # kph
sensor_paras['RADARx4']['v_noise_sigma_cf'] = [1.2, 2.0]    # kph
# count relative error
sensor_paras['RADARx4']['c_bias_p_ff'] = [0, 0]
sensor_paras['RADARx4']['c_bias_p_cf'] = [0, 0]
sensor_paras['RADARx4']['c_sigma_p_ff'] = [0, 0]
sensor_paras['RADARx4']['c_sigma_p_cf'] = [0, 0]
# --------------------------------- #


# ------------- RADARx8 ------------- #
sensor_paras['RADARx8'] = OrderedDict()

sensor_paras['RADARx8']['alpha'] = 84
sensor_paras['RADARx8']['theta'] = 2
sensor_paras['RADARx8']['s']= 0.91
sensor_paras['RADARx8']['L'] = 3.66    # m
sensor_paras['RADARx8']['offset_dist'] = 43.5    # m

sensor_paras['RADARx8']['v_range'] = [8.05, 159.33]   # kph [5, 99] mph
sensor_paras['RADARx8']['aggregation_sec'] = 30
sensor_paras['RADARx8']['awake_sec'] = 30

sensor_paras['RADARx8']['occlusion'] = 'True'
sensor_paras['RADARx8']['p_occlusion_accept'] = [30,30]
sensor_paras['RADARx8']['v_threshold'] = 64.37    # kph   40 mph
sensor_paras['RADARx8']['p_missing_ff'] = [0.0,0.0]
sensor_paras['RADARx8']['p_missing_cf'] = [1.5, 2]

sensor_paras['RADARx8']['noise_type'] = 'absolute'
# relative, NOT used
sensor_paras['RADARx8']['v_bias_ff']  = [0, 0]
sensor_paras['RADARx8']['v_bias_cf'] = [0, 0]
sensor_paras['RADARx8']['v_accuracy_p_ff'] = [0, 0]
sensor_paras['RADARx8']['v_accuracy_p_cf'] = [0, 0]
# use absolute error
sensor_paras['RADARx8']['v_noise_mu_ff'] = [0, 0]     # kph
sensor_paras['RADARx8']['v_noise_mu_cf'] = [0, 0]     # kph Not used
sensor_paras['RADARx8']['v_noise_sigma_ff'] = [0.1, 0.2]   # kph
sensor_paras['RADARx8']['v_noise_sigma_cf'] = [0.6, 1.0]    # kph
# count relative error
sensor_paras['RADARx8']['c_bias_p_ff'] = [0, 0]
sensor_paras['RADARx8']['c_bias_p_cf'] = [0, 0]
sensor_paras['RADARx8']['c_sigma_p_ff'] = [0, 0]
sensor_paras['RADARx8']['c_sigma_p_cf'] = [0, 0]
# --------------------------------- #


# =========================================================================================== #
# ================================== RTMS sensor ============================================ #
# =========================================================================================== #
# ------------- RTMSI57 -------------- #
sensor_paras['RTMSI57'] = OrderedDict()

# RTMS is mounted high and aims perpendicular to the road, hence there is no occlusion.
sensor_paras['RTMSI57']['alpha'] = 0
sensor_paras['RTMSI57']['theta'] = 0
sensor_paras['RTMSI57']['s']= 0.91
sensor_paras['RTMSI57']['L'] = 3.66
sensor_paras['RTMSI57']['offset_dist'] = 0.0    # m

sensor_paras['RTMSI57']['v_range'] = [0,160.93]    # kph   [0, 100] MPH
sensor_paras['RTMSI57']['awake_sec'] = 30
sensor_paras['RTMSI57']['aggregation_sec'] = 30

sensor_paras['RTMSI57']['occlusion'] = 'False'
sensor_paras['RTMSI57']['p_occlusion_accept'] = [100,100]
sensor_paras['RTMSI57']['v_threshold'] = 80.47     # ~50 mph
sensor_paras['RTMSI57']['p_missing_ff'] = [0.0, 0.0]     # %
sensor_paras['RTMSI57']['p_missing_cf'] = [2.0, 4.0]     # %

sensor_paras['RTMSI57']['noise_type'] = 'relative'
# relative
sensor_paras['RTMSI57']['v_bias_ff']  = [3.2, 8.0]
sensor_paras['RTMSI57']['v_bias_cf'] = [8.0, 12.0]    # kph [5, 7.5] mph
sensor_paras['RTMSI57']['v_accuracy_p_ff'] = [8, 12]
sensor_paras['RTMSI57']['v_accuracy_p_cf'] = [12, 18]
# absolute
sensor_paras['RTMSI57']['v_noise_mu_ff'] = [0, 0]
sensor_paras['RTMSI57']['v_noise_mu_cf'] = [0, 0]
sensor_paras['RTMSI57']['v_noise_sigma_ff'] = [0, 0]
sensor_paras['RTMSI57']['v_noise_sigma_cf'] = [0, 0]
# count relative error
sensor_paras['RTMSI57']['c_bias_p_ff'] = [0, 0]
sensor_paras['RTMSI57']['c_bias_p_cf'] = [0, 0]
sensor_paras['RTMSI57']['c_sigma_p_ff'] = [3, 7]
sensor_paras['RTMSI57']['c_sigma_p_cf'] = [7, 13]
# --------------------------------- #


# ------------- RTMS -------------- #
sensor_paras['RTMS'] = OrderedDict()

sensor_paras['RTMS']['alpha'] = 0
sensor_paras['RTMS']['theta'] = 0
sensor_paras['RTMS']['s']= 0.91
sensor_paras['RTMS']['L'] = 3.66
sensor_paras['RTMS']['offset_dist'] = 0.0    # m

sensor_paras['RTMS']['v_range'] = [0,160.93]    # kph   [0, 100] MPH
sensor_paras['RTMS']['awake_sec'] = 30
sensor_paras['RTMS']['aggregation_sec'] = 30

sensor_paras['RTMS']['occlusion'] = 'False'
sensor_paras['RTMS']['p_occlusion_accept'] = [100,100]
sensor_paras['RTMS']['v_threshold'] = 64.37    # kph   40 mph
sensor_paras['RTMS']['p_missing_ff'] = [0.0, 0.0]     # %
sensor_paras['RTMS']['p_missing_cf'] = [2.0, 4.0]     # %

sensor_paras['RTMS']['noise_type'] = 'relative'
# relative
sensor_paras['RTMS']['v_bias_ff']  = [0, 0]
sensor_paras['RTMS']['v_bias_cf'] = [0, 0]
sensor_paras['RTMS']['v_accuracy_p_ff'] = [8, 12]
sensor_paras['RTMS']['v_accuracy_p_cf'] = [12, 18]
# absolute
sensor_paras['RTMS']['v_noise_mu_ff'] = [0, 0]
sensor_paras['RTMS']['v_noise_mu_cf'] = [0, 0]
sensor_paras['RTMS']['v_noise_sigma_ff'] = [0, 0]
sensor_paras['RTMS']['v_noise_sigma_cf'] = [0, 0]
# count relative error
sensor_paras['RTMS']['c_bias_p_ff'] = [0, 0]
sensor_paras['RTMS']['c_bias_p_cf'] = [0, 0]
sensor_paras['RTMS']['c_sigma_p_ff'] = [3, 7]
sensor_paras['RTMS']['c_sigma_p_cf'] = [7, 13]
# --------------------------------- #


# ------------- RTMSx2 -------------- #
sensor_paras['RTMSx2'] = OrderedDict()

sensor_paras['RTMSx2']['alpha'] = 0
sensor_paras['RTMSx2']['theta'] = 0
sensor_paras['RTMSx2']['s']= 0.91
sensor_paras['RTMSx2']['L'] = 3.66
sensor_paras['RTMSx2']['offset_dist'] = 0.0    # m

sensor_paras['RTMSx2']['v_range'] = [0,160.93]    # kph   [0, 100] MPH
sensor_paras['RTMSx2']['awake_sec'] = 30
sensor_paras['RTMSx2']['aggregation_sec'] = 30

sensor_paras['RTMSx2']['occlusion'] = 'False'
sensor_paras['RTMSx2']['p_occlusion_accept'] = [100,100]
sensor_paras['RTMSx2']['v_threshold'] = 64.37    # kph   40 mph
sensor_paras['RTMSx2']['p_missing_ff'] = [0.0, 0.0]     # %
sensor_paras['RTMSx2']['p_missing_cf'] = [1.0, 2.0]     # %

sensor_paras['RTMSx2']['noise_type'] = 'relative'
# relative
sensor_paras['RTMSx2']['v_bias_ff']  = [0, 0]
sensor_paras['RTMSx2']['v_bias_cf'] = [0, 0]
sensor_paras['RTMSx2']['v_accuracy_p_ff'] = [4, 6]
sensor_paras['RTMSx2']['v_accuracy_p_cf'] = [6, 9]
# absolute
sensor_paras['RTMSx2']['v_noise_mu_ff'] = [0, 0]
sensor_paras['RTMSx2']['v_noise_mu_cf'] = [0, 0]
sensor_paras['RTMSx2']['v_noise_sigma_ff'] = [0, 0]
sensor_paras['RTMSx2']['v_noise_sigma_cf'] = [0, 0]
# count relative error
sensor_paras['RTMSx2']['c_bias_p_ff'] = [0, 0]
sensor_paras['RTMSx2']['c_bias_p_cf'] = [0, 0]
sensor_paras['RTMSx2']['c_sigma_p_ff'] = [1.5, 3.5]
sensor_paras['RTMSx2']['c_sigma_p_cf'] = [3.5, 6.5]
# --------------------------------- #



# ------------- RTMSx4 -------------- #
sensor_paras['RTMSx4'] = OrderedDict()

sensor_paras['RTMSx4']['alpha'] = 0
sensor_paras['RTMSx4']['theta'] = 0
sensor_paras['RTMSx4']['s']= 0.91
sensor_paras['RTMSx4']['L'] = 3.66
sensor_paras['RTMSx4']['offset_dist'] = 0.0    # m

sensor_paras['RTMSx4']['v_range'] = [0,160.93]    # kph   [0, 100] MPH
sensor_paras['RTMSx4']['awake_sec'] = 30
sensor_paras['RTMSx4']['aggregation_sec'] = 30

sensor_paras['RTMSx4']['occlusion'] = 'False'
sensor_paras['RTMSx4']['p_occlusion_accept'] = [100,100]
sensor_paras['RTMSx4']['v_threshold'] = 64.37    # kph   40 mph
sensor_paras['RTMSx4']['p_missing_ff'] = [0.0, 0.0]     # %
sensor_paras['RTMSx4']['p_missing_cf'] = [0.5, 1.0]     # %

sensor_paras['RTMSx4']['noise_type'] = 'relative'
# relative
sensor_paras['RTMSx4']['v_bias_ff']  = [0, 0]
sensor_paras['RTMSx4']['v_bias_cf'] = [0, 0]
sensor_paras['RTMSx4']['v_accuracy_p_ff'] = [2, 3]
sensor_paras['RTMSx4']['v_accuracy_p_cf'] = [3, 4.5]
# absolute
sensor_paras['RTMSx4']['v_noise_mu_ff'] = [0, 0]
sensor_paras['RTMSx4']['v_noise_mu_cf'] = [0, 0]
sensor_paras['RTMSx4']['v_noise_sigma_ff'] = [0, 0]
sensor_paras['RTMSx4']['v_noise_sigma_cf'] = [0, 0]
# count relative error
sensor_paras['RTMSx4']['c_bias_p_ff'] = [0, 0]
sensor_paras['RTMSx4']['c_bias_p_cf'] = [0, 0]
sensor_paras['RTMSx4']['c_sigma_p_ff'] = [0.75, 1.75]
sensor_paras['RTMSx4']['c_sigma_p_cf'] = [1.75, 3.25]
# --------------------------------- #

# ------------- RTMSx8 -------------- #
sensor_paras['RTMSx8'] = OrderedDict()

sensor_paras['RTMSx8']['alpha'] = 0
sensor_paras['RTMSx8']['theta'] = 0
sensor_paras['RTMSx8']['s']= 0.91
sensor_paras['RTMSx8']['L'] = 3.66
sensor_paras['RTMSx8']['offset_dist'] = 0.0    # m

sensor_paras['RTMSx8']['v_range'] = [0,160.93]    # kph   [0, 100] MPH
sensor_paras['RTMSx8']['awake_sec'] = 30
sensor_paras['RTMSx8']['aggregation_sec'] = 30

sensor_paras['RTMSx8']['occlusion'] = 'False'
sensor_paras['RTMSx8']['p_occlusion_accept'] = [100,100]
sensor_paras['RTMSx8']['v_threshold'] = 64.37    # kph   40 mph
sensor_paras['RTMSx8']['p_missing_ff'] = [0.0, 0.0]     # %
sensor_paras['RTMSx8']['p_missing_cf'] = [0.25, 0.5]     # %

sensor_paras['RTMSx8']['noise_type'] = 'relative'
# relative
sensor_paras['RTMSx8']['v_bias_ff']  = [0, 0]
sensor_paras['RTMSx8']['v_bias_cf'] = [0, 0]
sensor_paras['RTMSx8']['v_accuracy_p_ff'] = [1, 1.5]
sensor_paras['RTMSx8']['v_accuracy_p_cf'] = [1.5, 2.25]
# absolute
sensor_paras['RTMSx8']['v_noise_mu_ff'] = [0, 0]
sensor_paras['RTMSx8']['v_noise_mu_cf'] = [0, 0]
sensor_paras['RTMSx8']['v_noise_sigma_ff'] = [0, 0]
sensor_paras['RTMSx8']['v_noise_sigma_cf'] = [0, 0]
# count relative error
sensor_paras['RTMSx8']['c_bias_p_ff'] = [0, 0]
sensor_paras['RTMSx8']['c_bias_p_cf'] = [0, 0]
sensor_paras['RTMSx8']['c_sigma_p_ff'] = [0.375, 0.875]
sensor_paras['RTMSx8']['c_sigma_p_cf'] = [0.875, 1.625]
# --------------------------------- #


# =========================================================================================== #
# ==================================== ICONE sensor ========================================= #
# =========================================================================================== #
# ------------- ICONE ------------- #
sensor_paras['ICONE'] = OrderedDict()

# ICONEs should be pointed 30ft upstream to the close lane for every 1 ft to the side of road
sensor_paras['ICONE']['alpha'] = 84
sensor_paras['ICONE']['theta'] = 2
sensor_paras['ICONE']['s'] = 0.91
sensor_paras['ICONE']['L'] = 3.66
sensor_paras['ICONE']['offset_dist'] = 43.5    # m

sensor_paras['ICONE']['v_range'] = [8.05, 159.33]   # kph
sensor_paras['ICONE']['aggregation_sec'] = 60
sensor_paras['ICONE']['awake_sec'] = 30

sensor_paras['ICONE']['occlusion'] = 'True'
sensor_paras['ICONE']['p_occlusion_accept'] = [30,30]
sensor_paras['ICONE']['v_threshold'] = 64.37    # kph   40 mph
sensor_paras['ICONE']['p_missing_ff'] = [0.0,0.0]
sensor_paras['ICONE']['p_missing_cf'] = [13.0,17.0]

sensor_paras['ICONE']['noise_type'] = 'absolute'
# relative
sensor_paras['ICONE']['v_bias_ff']  = [0, 0]
sensor_paras['ICONE']['v_bias_cf'] = [0, 0]
sensor_paras['ICONE']['v_accuracy_p_ff'] = [0, 0]
sensor_paras['ICONE']['v_accuracy_p_cf'] = [0, 0]
# absolute
sensor_paras['ICONE']['v_noise_mu_ff'] = [0, 0]       # kph
sensor_paras['ICONE']['v_noise_mu_cf'] = [0, 0]
sensor_paras['ICONE']['v_noise_sigma_ff'] = [0.8, 1.6]    # kph
sensor_paras['ICONE']['v_noise_sigma_cf'] = [4.8, 8.0]       # kph  in literature was [40, 50] mph
# count relative error
sensor_paras['ICONE']['c_bias_p_ff'] = [0, 0]
sensor_paras['ICONE']['c_bias_p_cf'] = [0, 0]
sensor_paras['ICONE']['c_sigma_p_ff'] = [0, 0]
sensor_paras['ICONE']['c_sigma_p_cf'] = [0, 0]
# --------------------------------- #


# ------------- ICONEx2 ------------- #
sensor_paras['ICONEx2'] = OrderedDict()

sensor_paras['ICONEx2']['alpha'] = 84
sensor_paras['ICONEx2']['theta'] = 2
sensor_paras['ICONEx2']['s'] = 0.91
sensor_paras['ICONEx2']['L'] = 3.66
sensor_paras['ICONEx2']['offset_dist'] = 43.5    # m

sensor_paras['ICONEx2']['v_range'] = [8.05, 159.33]   # kph
sensor_paras['ICONEx2']['aggregation_sec'] = 60
sensor_paras['ICONEx2']['awake_sec'] = 30

sensor_paras['ICONEx2']['occlusion'] = 'True'
sensor_paras['ICONEx2']['p_occlusion_accept'] = [30,30]
sensor_paras['ICONEx2']['v_threshold'] = 64.37    # kph
sensor_paras['ICONEx2']['p_missing_ff'] = [0.0,0.0]
sensor_paras['ICONEx2']['p_missing_cf'] = [6.5,8.5]

sensor_paras['ICONEx2']['noise_type'] = 'absolute'
# relative
sensor_paras['ICONEx2']['v_bias_ff']  = [0, 0]
sensor_paras['ICONEx2']['v_bias_cf'] = [0, 0]
sensor_paras['ICONEx2']['v_accuracy_p_ff'] = [0, 0]
sensor_paras['ICONEx2']['v_accuracy_p_cf'] = [0, 0]
# absolute
sensor_paras['ICONEx2']['v_noise_mu_ff'] = [0.0, 0.0]       # kph
sensor_paras['ICONEx2']['v_noise_mu_cf'] = [0.0, 0.0]
sensor_paras['ICONEx2']['v_noise_sigma_ff'] = [0.4, 0.8]    # kph
sensor_paras['ICONEx2']['v_noise_sigma_cf'] = [2.4, 4.0]       # kph  in literature was [40, 50] mph
# count relative error
sensor_paras['ICONEx2']['c_bias_p_ff'] = [0, 0]
sensor_paras['ICONEx2']['c_bias_p_cf'] = [0, 0]
sensor_paras['ICONEx2']['c_sigma_p_ff'] = [0, 0]
sensor_paras['ICONEx2']['c_sigma_p_cf'] = [0, 0]
# --------------------------------- #



# =========================================================================================== #
# ==================================== BT sensor ============================================ #
# =========================================================================================== #
# --- Travel Time Ground Truth ---- #
sensor_paras['IDEALBT'] = OrderedDict()

sensor_paras['IDEALBT']['alpha'] = 0
sensor_paras['IDEALBT']['theta'] = 0
sensor_paras['IDEALBT']['s'] = 0.91
sensor_paras['IDEALBT']['L'] = 3.66

sensor_paras['IDEALBT']['aggregation_sec'] = 30
sensor_paras['IDEALBT']['p_penetration'] = [100,100]
sensor_paras['IDEALBT']['t_noise_sigma'] = [0,0]
sensor_paras['IDEALBT']['t_noise_mu'] = [0,0]
# --------------------------------- #


# ----------- BT ----------- #
sensor_paras['BT'] = OrderedDict()

sensor_paras['BT']['alpha'] = 0
sensor_paras['BT']['theta'] = 0
sensor_paras['BT']['s'] = 0.91
sensor_paras['BT']['L'] = 3.66

sensor_paras['BT']['aggregation_sec'] = 30
sensor_paras['BT']['p_penetration'] = [5,10]
sensor_paras['BT']['t_noise_sigma'] = [1,3]
sensor_paras['BT']['t_noise_mu'] = [0,0]
# --------------------------------- #

# ----------- BTx2 ----------- #
sensor_paras['BTx2'] = OrderedDict()

sensor_paras['BTx2']['alpha'] = 0
sensor_paras['BTx2']['theta'] = 0
sensor_paras['BTx2']['s'] = 0.91
sensor_paras['BTx2']['L'] = 3.66

sensor_paras['BTx2']['aggregation_sec'] = 30
sensor_paras['BTx2']['p_penetration'] = [10,20]
sensor_paras['BTx2']['t_noise_sigma'] = [0.5,1.5]
sensor_paras['BTx2']['t_noise_mu'] = [0,0]
# --------------------------------- #
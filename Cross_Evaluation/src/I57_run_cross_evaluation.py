import matplotlib.pyplot as plt
import numpy as np

from cross_evaluation import CrossEval

__author__ = 'Yanning Li and Juan Carlos Martinez'
"""
This script is used to run cross evaluation on I-57 work zone using the cross_evaluation class.
Uncomment proper sections to perform each analysis.
"""

# ========================================================================================= #
# configuration for I57 work zone
# ========================================================================================= #
# set the folders and directories
workzone = 'I57'
log_dir = '../'
data_dir = '../'
config_file = '../{0}_configurations_input.txt'.format(workzone)

replications = [51642]

# the estimation grid, (5,50) for linearFILL and fisb, (5,200) for enkf
estimation_grid = (5, 50)

# initialize the cross evaluation instance
cross_eval = CrossEval(workzone=workzone, log_dir=log_dir, data_dir=data_dir,
                        grid_res=estimation_grid, replications=replications)

# load the simulation configuration
cross_eval.load_config(config_file)

# ========================================================================================= #
# run estimators; it will generate the estimation results for the smart work zone configurations
# that has not been estimated before. No need to run if all already generated.
# ========================================================================================= #
if False:
    print('\n\nStart estimation...\n')
    cross_eval.fetch_data()
    cross_eval.run_estimators()

# ========================================================================================= #
# Visualize the measurements
# ========================================================================================= #
# Visualize the measurement in time space domain either by configuration
if False:
    cross_eval.plot_speed_meas_for_config(51642, 'configI57', unit='imperial', limit=[0, 80])

# Visualize the measurement in time space domain by a list of sensors
if False:
    cross_eval.plot_speed_meas_for_sensors(51642,
                                           ['IDEALLoc300m', 'IDEALLoc500m', 'IDEALLoc700m', 'IDEALLoc900m',
                                            'IDEALLoc1100m',
                                            'IDEALLoc1300m', 'IDEALLoc1500m', 'IDEALLoc1700m', 'IDEALLoc1900m',
                                            'IDEALLoc2100m',
                                            'IDEALLoc2300m', 'IDEALLoc2500m', 'IDEALLoc2700m', 'IDEALLoc2900m',
                                            'IDEALLoc3100m',
                                            'IDEALLoc3300m', 'IDEALLoc3500m', 'IDEALLoc3700m', 'IDEALLoc3900m',
                                            'IDEALLoc4100m',
                                            'IDEALLoc4300m', 'IDEALLoc4500m', 'IDEALLoc4700m', 'IDEALLoc4900m',
                                            'IDEALLoc5100m',
                                            'IDEALLoc5300m', 'IDEALLoc5500m', 'IDEALLoc5700m', 'IDEALLoc5900m',
                                            'IDEALLoc6100m',
                                            'IDEALLoc6300m', 'IDEALLoc6500m', 'IDEALLoc6700m', 'IDEALLoc6900m',
                                            'IDEALLoc7100m', ],
                                           unit='imperial', limit=[0, 80])

# Visualize the measurement in time series
if False:
    sensors_to_plot = ['RTMS']
    loc = [4000, 6000]
    sensor_ids = []
    for i in sensors_to_plot:
        for j in loc:
            sensor_ids.append('{0}Loc{1}m'.format(i, j))
    cross_eval.plot_speed_meas_time_series_sensors(51642, sensor_ids, unit='imperial', limit=(0,80))
    cross_eval.plot_flow_meas_time_series_sensors(51642, sensor_ids)

# ========================================================================================= #
# Calibrate fundamental diagram
# ========================================================================================= #
if False:
    cg_file_path_1 = '../Virtual_sensor_data/I57/rep1100/'
    cg_file_path_2 = '../Virtual_sensor_data/I57/rep1101/'
    cg_file_path_3 = '../Virtual_sensor_data/I57/rep1102/'
    ff_file_path_1 = '../Virtual_sensor_data/I57/rep1103/'
    ff_file_path_2 = '../Virtual_sensor_data/I57/rep1104/'
    sensor_files = []
    sensors_to_plot = ['IDEAL']
    loc = [1400, 1600]
    for i in sensors_to_plot:
        for j in loc:
            sensor_files.append('{0}{1}Loc{2}m.txt'.format(cg_file_path_1, i, j))
            sensor_files.append('{0}{1}Loc{2}m.txt'.format(cg_file_path_2, i, j))
            sensor_files.append('{0}{1}Loc{2}m.txt'.format(cg_file_path_3, i, j))
            sensor_files.append('{0}{1}Loc{2}m.txt'.format(ff_file_path_1, i, j))
            sensor_files.append('{0}{1}Loc{2}m.txt'.format(ff_file_path_2, i, j))

    agg_interval = 30
    cross_eval.calib_FD_using_sensor_data(sensor_files, agg_interval)

    # validate the calibration result using the FD true state data
    cross_eval.calib_FD_on_grid(grid=(5, 50), replications=[51642])

# ========================================================================================= #
# generate the true states from the trajectory file
# ========================================================================================= #
if False:
    grids = [(5, 50)]
    for grid in grids:
        cross_eval.generate_true_state(grid, 17.88)

# ========================================================================================= #
# compute the true queue length and plot true speed, density, queue, and travel time data
# ========================================================================================= #
# compute the true queue length
if False:
    grids = [(5, 50)]
    replications_to_compute = [51642]
    v_threshold = 17.88
    for rep in replications_to_compute:
        for grid in grids:
            cross_eval.compute_true_queue_for_grid(rep, grid, v_threshold)
            cross_eval.compute_trueinst_traveltime_for_grid(rep, grid, max_speed=27.2,
                                                            min_speed=1.34)

# plot the velocity, queue length and travel time.
if False:
    plot_speed = True
    plot_queue = True
    plot_tt = True

    grids = [(5, 50)]
    replications_to_plot = [51642]
    save_fig = True
    for grid in grids:
        if plot_speed:
            cross_eval.plot_true_speed_for_rep(grid, replications_to_plot[0], unit='imperial',
                                               limit=[0, 80], save_fig=save_fig, fig_size=(22, 8),
                                               fontsize=[50, 46, 42],
                                               title='I-57 True velocity field (mph)')
        # # cross_eval.plot_true_density_for_rep(grid, replications[0], unit='imperial', limit=[0,1000])
        # # cross_eval.plot_true_speed_field_partition(grid, replications[0], unit='imperial', limit=[0,80],
        # #                                           queue_buffer=800,save_fig=False,title='True speed 5s50m (mph)')
        if plot_queue:
            cross_eval.plot_true_queue_for_rep(grid, replications_to_plot[0], unit='imperial',
                                               fig_size=(22, 8), title='I-57 True queue length',
                                               fontsize=[50, 46, 42],
                                               limit=[0, 4], save_fig=save_fig)
        if plot_tt:
            cross_eval.plot_true_traveltime_for_rep(grid, replications_to_plot[0],
                                                    unit='metric', limit=[0, 8300],
                                                    title='I-57 True travel time',
                                                    fontsize=[50, 46, 42],
                                                    fig_size=(22, 8), save_fig=save_fig)
            # cross_eval.plot_measured_traveltime_for_rep(grid, replications[0], unit='metric', limit=[0, 8300])
            #
            # cross_eval.plot_true_speed_density_for_rep(grid, 41236, unit='metric', speed_limit=[0,40], density_limit=[0,0.12])
            # cross_eval.plot_trueinst_and_measured_traveltime_for_rep(grid, replications_to_compute, limit=[0,30],
            #                                                          save_fig=True, title=None)

# =====================================================================================================
# Compute the queue length and travel time for scenario
# =====================================================================================================
if False:
    # compute speed for scenarios
    num_sensors = [2, 3, 4, 5, 6, 7, 9, 12, 17, 33]
    types_sensors = ['IDEAL', 'RADAR', 'RADARx2', 'RADARx4', 'RADARx8',
                     'RTMS', 'RTMSx2', 'RTMSx4', 'RTMSx8', 'ICONE', 'ICONEx2']

    algs = ['enkfANupdated', 'fisb', 'linearFILL']

    v_threshold = 17.88
    replications_to_compute = [51642]

    for num in num_sensors:
        for ty in types_sensors:
            config_id = 'configHOMO_{0}_{1}'.format(num, ty)
            for alg in algs:

                for rep in replications_to_compute:
                    if 'enkf' in alg:
                        cross_eval.compute_queue_for_scenario(rep,
                                                              config_id, alg, v_threshold,
                                                              estimate_grid=(5, 200))
                        cross_eval.compute_traveltime_for_scenario(rep, config_id, alg,
                                                                   grid_res=(5, 200),
                                                                   max_speed=27.2, min_speed=1.34)
                    else:
                        cross_eval.compute_queue_for_scenario(rep,
                                                              config_id, alg, v_threshold,
                                                              estimate_grid=(5, 50))
                        cross_eval.compute_traveltime_for_scenario(rep, config_id, alg,
                                                                   grid_res=(5, 50),
                                                                   max_speed=27.2, min_speed=1.34)

# =====================================================================================================
# Visualize estimation results: speed, density, travel time, queue, and errors
# =====================================================================================================
# cross_eval.plot_speed_for_scenario( 51642, 'configI57', 'linearFILL',
#                                     unit='imperial', limit=[0,80], save_fig=True,
#                                     grid_res=(5,50),
#                                     fig_size=(22,8), fontsize=[36,34,32],
#                                     title='Velocity estimation using linearFILL (mph)')

if False:
    # plot speed for scenarios
    # num_sensors = [2, 3, 4, 5, 6, 7, 9, 12, 17, 33]
    # types_sensors = ['IDEAL', 'RADAR', 'RADARx2', 'RADARx4', 'RADARx8',
    #                 'RTMS', 'RTMSx2', 'RTMSx4', 'RTMSx8', 'ICONE', 'ICONEx2']
    num_sensors = [9]
    types_sensors = ['RTMS']
    algs = ['enkfANupdated', 'linearFILL', 'fisb']
    alg_labels = {'enkfANupdated': 'Kalman filter', 'linearFILL': 'spatial interpolation',
                  'fisb': 'spatio-temporal filtering'}
    save_fig = True
    # #
    for num in num_sensors:
        for ty in types_sensors:
            config_id = 'configHOMO_{0}_{1}'.format(num, ty)
            for alg in algs:
                cross_eval.plot_speed_for_scenario(51642, config_id, alg,
                                                   unit='imperial', limit=[0, 80], save_fig=save_fig,
                                                   grid_res=(5, 50),
                                                   fig_size=(22, 8), fontsize=[36, 34, 32],
                                                   title='Velocity estimation using {0} (mph)'.format(alg_labels[alg]))
                # cross_eval.plot_density_for_scenario( 51642,config_id, alg,
                #                                       unit='imperial', limit=[0,500], save_fig=save_fig)
                cross_eval.plot_queue_for_scenario(51642, config_id, alg,
                                                   unit='imperial', ylim=[0, 4], save_fig=save_fig,
                                                   true_grid=(5, 50), fig_size=(22, 8), fontsize=[36, 34, 32],
                                                   title='Queue estimation using {0}'.format(alg_labels[alg]))
                cross_eval.plot_traveltime_for_scenario(51642, config_id, alg, true_grid=(5, 50),
                                                        fig_size=(22, 8), fontsize=[36, 34, 32], ylim=(0, 60),
                                                        unit='imperial', save_fig=save_fig,
                                                        title='Travel time estimation using {0}'.format(
                                                            alg_labels[alg]))

                # cross_eval.plot_speed_est_error_for_scenario( replications[0], config_id, alg,
                #                                               unit='imperial', limit=[0,80], save_fig=save_fig)

# =====================================================================================================
# analyze the BAU case, compute the MAE errors
# =====================================================================================================
if False:
    updated = False
    computeMAE = True

    # first update the queue length and travel time estimation using standard methods
    if not updated:
        config_id = 'configI57'
        replications_BAU = [51642]
        v_threshold = 17.88
        algs = ['linearFILL', 'fisb', 'enkfANupdated']
        for rep in replications_BAU:
            for alg in algs:
                if 'enkf' in alg:
                    cross_eval.compute_queue_for_scenario(rep,
                                                          config_id, alg, v_threshold,
                                                          estimate_grid=(5, 200))
                    cross_eval.compute_traveltime_for_scenario(rep, config_id, alg,
                                                               grid_res=(5, 200),
                                                               max_speed=27.2, min_speed=1.34)
                else:
                    cross_eval.compute_queue_for_scenario(rep,
                                                          config_id, alg, v_threshold,
                                                          estimate_grid=(5, 50))
                    cross_eval.compute_traveltime_for_scenario(rep, config_id, alg,
                                                               grid_res=(5, 50),
                                                               max_speed=27.2, min_speed=1.34)
                print('Updated queue and travel time estimation for rep {0} and alg {1}'.format(rep, alg))

    # now compute the MAE errors for I80 BAU
    if computeMAE:

        computeMAE_speed = True
        computeMAE_speed_aq = True
        computeMAE_queue = True
        computeMAE_tt = True

        # configure the computation
        replications_BAU = [51642]
        config_id = 'configI57'
        queue_buffer = 800
        algs = ['linearFILL', 'fisb', 'enkfANupdated']
        norm = 'L1'

        # compute the average speed MAE
        if computeMAE_speed:
            errs = {}
            avg_err = {}
            for alg in algs:
                if alg not in errs.keys():
                    errs[alg] = []
                    avg_err[alg] = 0

                for rep in replications_BAU:
                    err = cross_eval.compute_speed_error_for_scenario_on_grid(rep, config_id, alg, (5, 50),
                                                                              norm=norm, area='all',
                                                                              queue_buffer=queue_buffer)
                    errs[alg].append(err)

                # compute the mean
                avg_err[alg] = np.mean(errs[alg])
                print('BAU velocity MAE with algorithm {0}: {1} mph; avg: {2} mph'.format(
                    alg,
                    np.array(errs[alg]) * 3600.0 / 1609,
                    avg_err[alg] * 3600.0 / 1609))

        # compute the average speed MAE around the queue
        if computeMAE_speed_aq:
            errs = {}
            avg_err = {}
            for alg in algs:
                if alg not in errs.keys():
                    errs[alg] = []
                    avg_err[alg] = 0

                for rep in replications_BAU:
                    err = cross_eval.compute_speed_error_for_scenario_on_grid(rep, config_id, alg, (5, 50),
                                                                              norm=norm, area='aroundqueue',
                                                                              queue_buffer=queue_buffer)
                    errs[alg].append(err)

                # compute the mean
                avg_err[alg] = np.mean(errs[alg])
                print('BAU velocity MAE around queue with algorithm {0}: {1} mph; avg: {2} mph'.format(
                    alg,
                    np.array(errs[alg]) * 3600.0 / 1609,
                    avg_err[alg] * 3600.0 / 1609))

        # compute the average speed MAE of the queue
        if computeMAE_queue:
            errs = {}
            avg_err = {}
            for alg in algs:
                if alg not in errs.keys():
                    errs[alg] = []
                    avg_err[alg] = 0

                for rep in replications_BAU:
                    err = cross_eval.compute_queuelength_error_for_scenario_on_grid(rep, config_id, alg, (5, 50),
                                                                                    norm=norm)
                    errs[alg].append(err)

                # compute the mean
                avg_err[alg] = np.mean(errs[alg])
                print('BAU queue MAE with algorithm {0}: {1} miles; avg: {2} miles'.format(
                    alg,
                    np.array(errs[alg]) / 1609.0,
                    avg_err[alg] / 1609.0))

        # compute the average speed MAE of the queue
        if computeMAE_tt:
            errs = {}
            avg_err = {}
            for alg in algs:
                if alg not in errs.keys():
                    errs[alg] = []
                    avg_err[alg] = 0

                for rep in replications_BAU:
                    err = cross_eval.compute_traveltime_error_for_scenario_on_grid(rep, config_id, alg,
                                                                                   (5, 50), norm='L1')
                    errs[alg].append(err)

                # compute the mean
                avg_err[alg] = np.mean(errs[alg])
                print('BAU travel time MAE with algorithm {0}: {1} min; avg: {2} min'.format(
                    alg,
                    np.array(errs[alg]) / 60.0,
                    avg_err[alg] / 60.0))

# =====================================================================================================
# Visualization and cross analysis of results
# =====================================================================================================
# ============================================
# # along number of sensors
if False:
    plot_speed = True
    plot_speed_aq = False
    plot_queue = False
    plot_tt = False

    replications_to_average = [51642]
    compare_grid = (5, 50)
    areas = ['all', 'aroundqueue']
    config_ids = []
    num_sensors = [2, 3, 5, 6, 7, 9, 12, 17, 33]  # 4 sensors configuration is removed from visualization
    alg_ids = ['linearFILL', 'fisb', 'enkfANupdated']
    norm = 'L1'
    x_axis = 'config'
    xlabel = 'Spacing of sensors (mile)'
    plot_cost_eff = False
    save_fig = True
    plot_style = 'line'

    sensor_type = 'ICONE'

    for num in num_sensors:
        config_ids.append('configHOMO_{0}_{1}'.format(num, sensor_type))
    # spacings.append( round((8200.0/num)/1609.0,2) )

    spacings = ['4', '2', '1', '6/8', '5/8',
                '4/8', '3/8', '2/8', '1/8']

    if plot_speed:
        cross_eval.compare_speed_error(replications_to_average, config_ids, alg_ids, estimation_grid,
                                       norm=norm, area='all', queue_buffer=800, ff_speed=27.2,
                                       x_axis=x_axis, xlabel=xlabel, plot_style=plot_style,
                                       xticklabels=spacings, unit='imperial', save_fig=save_fig,
                                       title='MAE of velocity', ylim=[4, 30], figsize=(22, 10), fontsize=[40, 36, 32],
                                       plot_cost_eff=plot_cost_eff, cost_eff_num_sensors=num_sensors,
                                       cost_eff_title='Improvement of velocity MAE per sensor',
                                       cost_eff_fontsize=[40, 36, 32],
                                       save_fig_name='I57_{0}_{1}_speed_all.pdf'.format(norm, sensor_type))

    # cross_eval.compare_speed_error(replications[0], config_ids, alg_ids, estimation_grid,
    #                                norm=norm, area='freeflow', queue_buffer=800,
    #                                x_axis=x_axis, xlabel=xlabel, plot_style=plot_style,
    #                                xticklabels=spacings, unit='imperial', save_fig=save_fig,
    #                                plot_cost_eff=plot_cost_eff, cost_eff_num_sensors=num_sensors,
    #                                title='MAE of speed in free flow',
    #                                cost_eff_title='Cost-effectiveness of free flow speed estimation',
    #                                save_fig_name='{0}_{1}_speed_free.pdf'.format(norm, sensor_type))
    # cross_eval.compare_speed_error(replications[0], config_ids, alg_ids, estimation_grid,
    #                                norm=norm, area='congflow', queue_buffer=800,
    #                                x_axis=x_axis, xlabel=xlabel, plot_style=plot_style,
    #                                xticklabels=spacings, unit='imperial', save_fig=save_fig,
    #                                plot_cost_eff=plot_cost_eff, cost_eff_num_sensors=num_sensors,
    #                                title='MAE of speed at congested flow',
    #                                cost_eff_title='Cost-effectiveness of congested flow speed estimation',
    #                                save_fig_name='{0}_{1}_speed_cong.pdf'.format(norm, sensor_type))

    if plot_speed_aq:
        cross_eval.compare_speed_error(replications_to_average, config_ids, alg_ids, estimation_grid,
                                       norm=norm, area='aroundqueue', queue_buffer=800, ylim=[5, 30],
                                       x_axis=x_axis, xlabel=xlabel, plot_style=plot_style,
                                       xticklabels=spacings, unit='imperial', save_fig=save_fig,
                                       title='MAE of velocity around the queue', figsize=(22, 10),
                                       fontsize=[40, 36, 32],
                                       plot_cost_eff=plot_cost_eff, cost_eff_num_sensors=num_sensors,
                                       cost_eff_fontsize=[40, 36, 32],
                                       cost_eff_title='Improvement of velocity around queue MAE per sensor',
                                       save_fig_name='I57_{0}_{1}_speed_aroundqueue.pdf'.format(norm, sensor_type))
    if plot_queue:
        cross_eval.compare_queuelength_error(replications_to_average, config_ids, alg_ids, compare_grid,
                                             norm=norm, x_axis=x_axis, xlabel=xlabel,
                                             xticklabels=spacings, unit='imperial', save_fig=save_fig,
                                             plot_style=plot_style, ylim=[0, 1], figsize=(22, 10),
                                             title='MAE of queue', fontsize=[40, 36, 32],
                                             plot_cost_eff=plot_cost_eff, cost_eff_num_sensors=num_sensors,
                                             cost_eff_fontsize=[40, 36, 32],
                                             cost_eff_title='Improvement of queue length MAE per sensor',
                                             save_fig_name='I57_{0}_{1}_queue.pdf'.format(norm, sensor_type))
    if plot_tt:
        cross_eval.compare_traveltime_error(replications_to_average, config_ids, alg_ids, estimation_grid,
                                            norm=norm, x_axis=x_axis, xlabel=xlabel, xticklabels=spacings,
                                            unit='metric', plot_style=plot_style, save_fig=save_fig,
                                            title='MAE of travel time', ylim=[0, 25], figsize=(22, 10),
                                            fontsize=[40, 36, 32],
                                            plot_cost_eff=plot_cost_eff, cost_eff_num_sensors=num_sensors,
                                            cost_eff_title='Improvement of travel time MAE per sensor',
                                            cost_eff_fontsize=[40, 36, 32],
                                            save_fig_name='I57_{0}_{1}_tt.pdf'.format(norm, sensor_type))

# ============================================
# along types of sensors
if False:
    plot_speed = True
    plot_speed_aq = True
    plot_queue = True
    plot_tt = True

    replications_to_average = [51642]
    config_ids = []
    types_sensors = ['ICONE', 'RADAR', 'RTMS', 'IDEAL']
    alg_ids = ['linearFILL', 'fisb', 'enkfANupdated']
    num_sensors = 9

    for sensor in types_sensors:
        config_ids.append('configHOMO_{0}_{1}'.format(num_sensors, sensor))

    norm = 'L1'
    x_axis = 'config'
    xlabel = 'Types of sensors'
    xticklabels = ['LER', 'RADAR', 'RTMS', 'IDEAL']
    plot_cost_eff = False
    save_fig = True
    plot_style = 'bar'

    if plot_speed:
        cross_eval.compare_speed_error(replications_to_average, config_ids, alg_ids, estimation_grid,
                                       norm=norm, area='all', queue_buffer=800, figsize=(11, 10),
                                       x_axis=x_axis, xlabel=xlabel, plot_style=plot_style,
                                       xticklabels=xticklabels, unit='imperial', save_fig=save_fig,
                                       plot_cost_eff=plot_cost_eff, cost_eff_num_sensors=num_sensors,
                                       title='MAE of velocity', ylim=[0, 30], ff_speed=27.2, fontsize=[40, 36, 36],
                                       cost_eff_title='Cost-effectiveness of speed estimation',
                                       save_fig_name='I57_{0}_{1}_speed_all.pdf'.format(norm, num_sensors))
    # cross_eval.compare_speed_error(replications, config_ids, alg_ids, estimation_grid,
    #                                norm=norm, area='freeflow', queue_buffer=800,
    #                                x_axis=x_axis, xlabel=xlabel, plot_style=plot_style,
    #                                xticklabels=xticklabels, unit='imperial', save_fig=save_fig,
    #                                plot_cost_eff=plot_cost_eff, cost_eff_num_sensors=num_sensors,
    #                                title='MAE of speed in free flow',
    #                                cost_eff_title='Cost-effectiveness of free flow speed estimation',
    #                                save_fig_name='{0}_{1}_speed_free.pdf'.format(norm, num_sensors))
    # cross_eval.compare_speed_error(replications[0], config_ids, alg_ids, estimation_grid,
    #                                norm=norm, area='congflow', queue_buffer=800,
    #                                x_axis=x_axis, xlabel=xlabel, plot_style=plot_style,
    #                                xticklabels=xticklabels, unit='imperial', save_fig=save_fig,
    #                                plot_cost_eff=plot_cost_eff, cost_eff_num_sensors=num_sensors,
    #                                title='MAE of speed at congested flow',
    #                                cost_eff_title='Cost-effectiveness of congested flow speed estimation',
    #                                save_fig_name='{0}_{1}_speed_cong.pdf'.format(norm, num_sensors))
    if plot_speed_aq:
        cross_eval.compare_speed_error(replications_to_average, config_ids, alg_ids, estimation_grid,
                                       norm=norm, area='aroundqueue', queue_buffer=800, figsize=(11, 10),
                                       x_axis=x_axis, xlabel=xlabel, plot_style=plot_style, ff_speed=27.2,
                                       xticklabels=xticklabels, unit='imperial', save_fig=save_fig,
                                       plot_cost_eff=plot_cost_eff, cost_eff_num_sensors=num_sensors,
                                       title='MAE of velocity around the queue', ylim=[0, 40],
                                       cost_eff_title='Cost-effectiveness of speed estimation around the queue',
                                       save_fig_name='I57_{0}_{1}_speed_aroundqueue.pdf'.format(norm, num_sensors))
    if plot_queue:
        cross_eval.compare_queuelength_error(replications_to_average, config_ids, alg_ids, estimation_grid,
                                             norm=norm, x_axis=x_axis, xlabel=xlabel, figsize=(11, 10),
                                             xticklabels=xticklabels, unit='imperial', save_fig=save_fig,
                                             plot_cost_eff=plot_cost_eff, cost_eff_num_sensors=num_sensors,
                                             plot_style=plot_style, ylim=[0, 0.5],
                                             title='MAE of queue', fontsize=[40, 36, 36],
                                             save_fig_name='I57_{0}_{1}_queue.pdf'.format(norm, num_sensors))
    if plot_tt:
        cross_eval.compare_traveltime_error(replications_to_average, config_ids, alg_ids, estimation_grid,
                                            norm=norm, x_axis=x_axis, xlabel=xlabel, xticklabels=xticklabels,
                                            unit='metric', plot_style=plot_style, save_fig=save_fig,
                                            plot_cost_eff=plot_cost_eff, cost_eff_num_sensors=num_sensors,
                                            title='MAE of travel time', figsize=(11, 10), fontsize=[40, 36, 36],
                                            save_fig_name='I57_{0}_{1}_tt.pdf'.format(norm, num_sensors))

# ============================================
# along accuracy of sensors
if False:
    plot_speed = True
    plot_speed_aq = True
    plot_queue = True
    plot_tt = True

    replications_to_average = [51642]

    config_ids = []
    # types_sensors = ['RTMS', 'RTMSx2', 'RTMSx4', 'RTMSx8']
    types_sensors = ['RADAR', 'RADARx2', 'RADARx4', 'RADARx8']
    alg_ids = ['linearFILL', 'fisb', 'enkfANupdated']
    num_sensors = 5

    for sensor in types_sensors:
        config_ids.append('configHOMO_{0}_{1}'.format(num_sensors, sensor))

    norm = 'L1'
    x_axis = 'config'
    xlabel = 'Accuracy of sensors'
    plot_cost_eff = False
    save_fig = True
    plot_style = 'bar'

    num_sensors = '5RADAR'

    if plot_speed:
        cross_eval.compare_speed_error(replications_to_average, config_ids, alg_ids, estimation_grid,
                                       norm=norm, area='all', queue_buffer=800, ylim=(0, 30),
                                       x_axis=x_axis, xlabel=xlabel, plot_style=plot_style,
                                       xticklabels=types_sensors, unit='imperial', save_fig=save_fig,
                                       plot_cost_eff=plot_cost_eff, cost_eff_num_sensors=num_sensors,
                                       title='MAE of velocity', figsize=(11, 10), fontsize=[38, 32, 32],
                                       cost_eff_title='Cost-effectiveness of speed estimation',
                                       save_fig_name='{0}_{1}_speed_all.pdf'.format(norm, num_sensors))
    # cross_eval.compare_speed_error(replications[0], config_ids, alg_ids, estimation_grid,
    #                                norm=norm, area='freeflow', queue_buffer=800,
    #                                x_axis=x_axis, xlabel=xlabel, plot_style=plot_style,
    #                                xticklabels=types_sensors, unit='imperial', save_fig=save_fig,
    #                                plot_cost_eff=plot_cost_eff, cost_eff_num_sensors=num_sensors,
    #                                title='MAE of speed in free flow',
    #                                cost_eff_title='Cost-effectiveness of free flow speed estimation',
    #                                save_fig_name='{0}_{1}_speed_free.pdf'.format(norm, num_sensors))
    # cross_eval.compare_speed_error(replications[0], config_ids, alg_ids, estimation_grid,
    #                                norm=norm, area='congflow', queue_buffer=800,
    #                                x_axis=x_axis, xlabel=xlabel, plot_style=plot_style,
    #                                xticklabels=types_sensors, unit='imperial', save_fig=save_fig,
    #                                plot_cost_eff=plot_cost_eff, cost_eff_num_sensors=num_sensors,
    #                                title='MAE of speed at congested flow',
    #                                cost_eff_title='Cost-effectiveness of congested flow speed estimation',
    #                                save_fig_name='{0}_{1}_speed_cong.pdf'.format(norm, num_sensors))
    if plot_speed_aq:
        cross_eval.compare_speed_error(replications_to_average, config_ids, alg_ids, estimation_grid,
                                       norm=norm, area='aroundqueue', queue_buffer=800, figsize=(11, 10),
                                       x_axis=x_axis, xlabel=xlabel, plot_style=plot_style, ylim=(0, 40),
                                       xticklabels=types_sensors, unit='imperial', save_fig=save_fig,
                                       plot_cost_eff=plot_cost_eff, cost_eff_num_sensors=num_sensors,
                                       title='MAE of velocity around the queue', fontsize=[38, 32, 32],
                                       cost_eff_title='Cost-effectiveness of speed estimation around the queue',
                                       save_fig_name='{0}_{1}_speed_aroundqueue.pdf'.format(norm, num_sensors))

    if plot_queue:
        cross_eval.compare_queuelength_error(replications_to_average, config_ids, alg_ids, estimation_grid,
                                             norm=norm, x_axis=x_axis, xlabel=xlabel,
                                             xticklabels=types_sensors, unit='imperial', save_fig=save_fig,
                                             plot_cost_eff=plot_cost_eff, cost_eff_num_sensors=num_sensors,
                                             plot_style=plot_style, figsize=(11, 10),
                                             title='MAE of queue', ylim=(0, 1.5), fontsize=[38, 32, 32],
                                             save_fig_name='{0}_{1}_queue.pdf'.format(norm, num_sensors))

    if plot_tt:
        cross_eval.compare_traveltime_error(replications_to_average, config_ids, alg_ids, estimation_grid,
                                            norm=norm, x_axis=x_axis, xlabel=xlabel, xticklabels=types_sensors,
                                            unit='metric', plot_style=plot_style, save_fig=save_fig, ylim=(0, 25),
                                            plot_cost_eff=plot_cost_eff, cost_eff_num_sensors=num_sensors,
                                            title='MAE of travel time', figsize=(11, 10), fontsize=[38, 32, 32],
                                            save_fig_name='{0}_{1}_tt.pdf'.format(norm, num_sensors))

# ============================================
# cost - accuracy plot
if False:
    plot_speed = True
    plot_speed_aq = True
    plot_queue = True
    plot_tt = True

    replications_to_average = [51642]
    compare_grid = (5, 50)

    num_sensors = [2, 3, 5, 6, 7, 9, 12, 17, 33]  # 4 sensors configuration is removed from visualization
    alg_ids = ['linearFILL', 'fisb', 'enkfANupdated']
    norm = 'L1'
    x_axis = 'config'
    xlabel = 'Monthly cost of system (thousand $)'
    save_fig = True
    plot_style = 'line'
    sensor_types = ['RTMS', 'RADAR']

    # cost of RTMS from VerMac: 600~1000
    # cost of RADAR from VerMac: $400~700.
    # normalize to RADAR
    unit_cost = {'RTMS': 1000, 'RADAR': 700}

    config_ids_prefix = []
    for num in num_sensors:
        config_ids_prefix.append('configHOMO_{0}'.format(num))

    if plot_queue:
        cross_eval.cost_compare_queue_error(replications_to_average, config_ids_prefix, alg_ids, estimation_grid,
                                            norm=norm, figsize=(22, 10), sensors=sensor_types, unit_costs=unit_cost,
                                            x_axis=x_axis, xlabel=xlabel, plot_style=plot_style, ylim=(0, 1),
                                            xticklabels=None, unit='imperial', save_fig=save_fig,
                                            title='Cost - accuarcy of queue estimation', fontsize=[38, 32, 32],
                                            save_fig_name='cost_queue.pdf')

    if plot_tt:
        cross_eval.cost_compare_traveltime_error(replications_to_average, config_ids_prefix, alg_ids, estimation_grid,
                                                 norm=norm, figsize=(22, 10), sensors=sensor_types,
                                                 unit_costs=unit_cost,
                                                 x_axis=x_axis, xlabel=xlabel, plot_style=plot_style, ylim=(0, 35),
                                                 xticklabels=None, unit='imperial', save_fig=save_fig,
                                                 title='Cost - accuarcy of travel time estimation',
                                                 fontsize=[38, 32, 32],
                                                 save_fig_name='cost_tt.pdf')

    if plot_speed:
        cross_eval.cost_compare_speed_error(replications_to_average, config_ids_prefix, alg_ids, estimation_grid,
                                            norm=norm, area='all', queue_buffer=800,
                                            figsize=(22, 10), sensors=sensor_types, unit_costs=unit_cost,
                                            x_axis=x_axis, xlabel=xlabel, plot_style=plot_style, ylim=(0, 34),
                                            xticklabels=None, unit='imperial', save_fig=save_fig,
                                            title='Cost - Accuracy of velocity', fontsize=[38, 32, 32],
                                            save_fig_name='cost_velocity.pdf')

    if plot_speed_aq:
        cross_eval.cost_compare_speed_error(replications_to_average, config_ids_prefix, alg_ids, estimation_grid,
                                            norm=norm, area='aroundqueue', queue_buffer=800,
                                            figsize=(22, 10), sensors=sensor_types, unit_costs=unit_cost,
                                            x_axis=x_axis, xlabel=xlabel, plot_style=plot_style, ylim=(0, 37),
                                            xticklabels=None, unit='imperial', save_fig=save_fig,
                                            title='Cost - Accuracy of velocity around queue', fontsize=[38, 32, 32],
                                            save_fig_name='cost_velocity_aroundqueue.pdf')

plt.show()

# ============================================
# along the number of ensembles
if False:
    algorithms = ['enkf50AN', 'enkf100AN', 'enkf150AN', 'enkf200AN', 'enkf300AN', 'enkf400AN', 'enkf500AN']
    config_ids = ['configHOMO_9_IDEAL', 'configHOMO_9_RTMS', 'configHOMO_9_RADAR', 'configHOMO_9_ICONE']

    cross_eval.compare_speed_error(replications[0], config_ids, algorithms, estimation_grid,
                                   norm='L1', area='all', queue_buffer=800,
                                   x_axis='alg', xlabel='Number of ensembles',
                                   xticklabels=[i.strip('enkf').strip('AN') for i in algorithms ],
                                   unit='imperial', save_fig=False, plot_cost_eff=False)

    cross_eval.compare_queuelength_error(replications[0], config_ids, algorithms, estimation_grid,
                                   norm='L1', x_axis='alg', xlabel='Number of ensembles',
                                   xticklabels=[i.strip('enkf').strip('AN') for i in algorithms ],
                                   unit='imperial', save_fig=False, plot_cost_eff=False)


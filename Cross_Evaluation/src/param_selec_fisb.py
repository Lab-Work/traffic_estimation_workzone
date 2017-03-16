import numpy as np
import os
import time
import matplotlib.pyplot as plt
from FISB_workzone import FISB_workzone
from cross_evaluation import Cross_Eval

"""
this script is used to obtain the optimal parameters sigma and tau for the FISB algorithm
these are set on the range [agg_sec/2, agg_sec*2] and [agg_space/2, agg_space*2] based on
the original author's recommendations

set all of the parameters under the USER INPUT section
true states for the required grid must have been generated beforehand
"""
__author__ = 'Juan Carlos Martinez'


# ================================================== #
#                      USER INPUT                    #
# ================================================== #

# directories
log_dir = '/Users/jotaceemeeme/Desktop/GitHub/IDOT_code/'
data_dir = '/Users/jotaceemeeme/Desktop/GitHub/IDOT_code/'

# sensor configuration id
config_id = 'configHOMO_9_IDEAL'
sensor_locs = ['100', '1100', '2100', '3100', '4100', '5300', '6300', '7300', '8300']
sensor_types = ['vol_gt', 'vol_gt', 'vol_gt', 'vol_gt', 'vol_gt', 'vol_gt', 'vol_gt', 'vol_gt', 'vol_gt']
sensor_sects = ['24671', '15814', '15814', '15589', '3412', '3412', '3412', '3401', '40842']
sensor_dists = ['22.0', '171.13', '1171.13', '745.08', '734.97', '1934.97', '2934.97', '34.5', '698.75']
sensor_flow_std_specs = ['0.1', '0.1', '0.1', '0.1', '0.1', '0.1', '0.1', '0.1', '0.1']
sensor_speed_std_specs = ['0.5', '0.5', '0.5', '0.5', '0.5', '0.5', '0.5', '0.5', '0.5']
sensor_agg_sec = ['30', '30', '30', '30', '30', '30', '30', '30', '30']

# work zone, replications and estimation grid
work_zone = 'I80'
replication = [41366]  # in square brackets [#]
grid = (5, 50)  # (time, space)

# algorithm parameters
queue_threshold = 17.88
omg_c = -15
omg_f = 80
v_crit = 60
delta_v = 20
lag = 0
time_diff_cutoff = 150

# tuning parameters
num_pars = 8  # number of parameters between [agg/2, agg*2] considered
queue_reg_len = 800  # [m] on each side of end-of-queue used for error

# ================================================== #
#           NO USER INPUT AFTER THIS LINE            #
# ================================================== #


def plot_speed_with_queue(speed_file, queue_file, queue_reg_len, grid):

    # extract data
    speed = np.flipud(np.loadtxt(speed_file, delimiter=',', dtype=float).T)
    queue = np.loadtxt(queue_file, delimiter=',', dtype=float)

    # set queue boundary points
    queue_pts = np.array([[(queue_end - queue_reg_len)//grid[1], queue_end//grid[1],
                           (queue_end + queue_reg_len)//grid[1]]
                          for queue_end in queue]).clip(min=0, max=speed.shape[1]-1).T

    # plot speed map
    fig = plt.figure(figsize=(15, 8), dpi=100)
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    im = ax.imshow(speed, cmap=plt.get_cmap('jet_r'), interpolation='nearest', aspect='auto')
    cbar = fig.colorbar(im)
    cbar.set_label('Speed [m/s]')
    ax.autoscale(False)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Distance [m]')
    ax.set_title('Speed Map [Queue Region: +/- {0}m]'.format(queue_reg_len))

    # plot queue boundaries
    for t in range(0, queue_pts.shape[1]):
        us_pt = queue_pts[0, t]
        q_pt = queue_pts[1, t]
        ds_pt = queue_pts[2, t]
        ax.plot([t - 0.5, t + 0.5], [us_pt - 0.5, us_pt - 0.5], color='k', linewidth=2)
        ax.plot([t - 0.5, t + 0.5], [q_pt - 0.5, q_pt - 0.5], color='k', linewidth=2)
        ax.plot([t - 0.5, t + 0.5], [ds_pt - 0.5, ds_pt - 0.5], color='k', linewidth=2)




# compute parameter lists
tau_facs = np.linspace(1/2, 2, num_pars)
sigma_facs = np.linspace(1/2, 2, num_pars)

# initialize cross evaluation class
cross_eval = Cross_Eval(workzone=work_zone, log_dir=log_dir, data_dir=data_dir, grid_res=grid, replications=replication)

# true states files
speed_true_state_file = data_dir + 'Estimation_results/{0}/rep{1}/truestate/truestate_{2}s{3}m_speed.txt'.\
    format(work_zone, replication[0], grid[0], grid[1])
queue_true_state_file = data_dir + 'Estimation_results/{0}/rep{1}/truestate/truestate_{2}s{3}m_queue.txt'.\
    format(work_zone, replication[0], grid[0], grid[1])

# try to open files
if not os.path.exists(speed_true_state_file) or not os.path.exists(queue_true_state_file):
    print('True states have not been generated...')
else:

    # write configuration file
    config_file_name = log_dir + 'fisb_opt_I80_configuration_inputs.txt'
    with open(config_file_name, 'w+') as config_file:

        # write configuration file for fisb optimization
        # iterate over every combination of parameters
        ct = 33  # config counter
        facs = []
        for sigma_fac in sigma_facs[4:]:
            for tau_fac in tau_facs:

                config_file.write('config_id:{0}_temp_{1}\n'.format(config_id, ct))
                # sensor parameters
                for idx in range(0, len(sensor_locs)):
                    line = 'sensor;id:IDEALLoc{0}m;'.format(sensor_locs[idx])
                    line += 'type:{0};'.format(sensor_types[idx])
                    line += 'section:{0};'.format(sensor_sects[idx])
                    line += 'distance:{0};'.format(sensor_dists[idx])
                    line += 'flow_std_specs:{0};'.format(sensor_flow_std_specs[idx])
                    line += 'speed_std_specs:{0};'.format(sensor_speed_std_specs[idx])
                    line += 'aggregation_sec:{0}'.format(sensor_agg_sec[idx])
                    line += '\n'
                    config_file.write(line)
                # algorithm parameters
                line = 'algorithm;id:fisb;type:fisb;'
                line += 'queue_threshold:{0};'.format(queue_threshold)
                line += 'omg_c:{0};'.format(omg_c)
                line += 'omg_f:{0};'.format(omg_f)
                line += 'sigma_fac:{0};'.format(sigma_fac)
                line += 'tau_fac:{0};'.format(tau_fac)
                line += 'v_crit:{0};'.format(v_crit)
                line += 'delta_v:{0};'.format(delta_v)
                line += 'lag:{0};'.format(lag)
                line += 'time_diff_cutoff:{0}'.format(time_diff_cutoff)
                line += '\n'
                config_file.write(line)
                config_file.write('\n')

                # update config counter
                facs.append([sigma_fac, tau_fac])
                ct += 1

    # close and load finished config file
    cross_eval.load_config(config_file_name)
    cross_eval.run_estimators()

    # metric for error is l1 error norm
    min_l1_error = np.inf
    opt_sigma_fac = None
    opt_tau_fac = None
    opt_temp_config_num = None

    # plot_speed_with_queue(speed_true_state_file, queue_true_state_file, queue_reg_len, grid)

    # true states from file
    speed_true_state = np.loadtxt(speed_true_state_file, delimiter=',', dtype=float)
    queue_true_state = np.loadtxt(queue_true_state_file, delimiter=',', dtype=float)

    # compute boundaries for end of queue region
    queue_reg_bounds = np.array([[queue_end - queue_reg_len, queue_end + queue_reg_len]
                                 for queue_end in queue_true_state]).clip(min=0, max=grid[1]*speed_true_state.shape[1])

    # matrix to show error computation
    l1_errors = np.zeros((len(sigma_facs), len(tau_facs)))
    l1_errors[:] = np.nan

    # keep track of file names for every configuration
    speed_file_names = []
    queue_file_names = []

    # find optimal parameters using all configurations generated
    for temp_config_num in range(1, ct):
        speed_file_name = (data_dir + 'Estimation_results/{0}/rep{1}/{2}s{3}m/'.format(work_zone, replication[0],
                                                                                       grid[0], grid[1]) + config_id +
                           '_temp_{0}_fisb_speed.txt'.format(temp_config_num))
        queue_file_name = (data_dir + 'Estimation_results/{0}/rep{1}/{2}s{3}m/'.format(work_zone, replication[0],
                                                                                       grid[0], grid[1]) + config_id +
                           '_temp_{0}_fisb_queue.txt'.format(temp_config_num))
        speed_file_names.append(speed_file_name)
        queue_file_names.append(queue_file_name)

        # load speed map from configuration
        speed = np.loadtxt(speed_file_name, delimiter=',', dtype=float)

        diff_mat = np.abs(speed_true_state - speed)

        error_mat = np.zeros(diff_mat.shape)
        error_mat[:] = np.nan

        for t_idx in range(0, queue_reg_bounds.shape[0]):
            queue_reg_idxs = np.arange(queue_reg_bounds[t_idx, 0]//grid[1],
                                       queue_reg_bounds[t_idx, 1]//grid[1], 1, dtype=int)
            queue_reg_idxs = queue_reg_bounds.shape[1] - queue_reg_idxs
            error_mat[t_idx, queue_reg_idxs] = diff_mat[t_idx, queue_reg_idxs]

        # update error
        l1_error = np.nanmean(error_mat)
        l1_errors[(temp_config_num-1) % len(sigma_facs), (temp_config_num-1)//len(sigma_facs)] = l1_error
        if l1_error < min_l1_error:
            min_l1_error = l1_error
            opt_sigma_fac = facs[temp_config_num - 1][0]
            opt_tau_fac = facs[temp_config_num - 1][1]
            opt_temp_config_num = temp_config_num

    fig, ax = plt.subplots()
    tau_diff = tau_facs[1] - tau_facs[0]
    sigma_diff = sigma_facs[1] - sigma_facs[0]
    im = ax.imshow(l1_errors, origin='lower', interpolation='nearest',
                   extent=[-tau_diff/2 + tau_facs[0], tau_facs[-1] + tau_diff/2,
                           -sigma_diff/2 + sigma_facs[0], sigma_facs[-1] + sigma_diff/2])
    ax.set_title('l1 error [grid:({0}s, {1}m), queue buffer: +/- {2}m]'.format(grid[0], grid[1], queue_reg_len))
    ax.set_ylabel('Sigma Factors')
    ax.set_yticks(sigma_facs)
    ax.set_xlabel('Tau Factors')
    ax.set_xticks(tau_facs)
    cbar = fig.colorbar(im)
    cbar.set_label('l1 error')
    plt.show()

    print('Optimal sigma factor: {0}'.format(opt_sigma_fac))
    print('Optimal tau factor: {0}'.format(opt_tau_fac))
    print('l1 error norm: {0}'.format(min_l1_error))
    print('Optimal temp configuration: {0}'.format(opt_temp_config_num))

    # plot_speed_with_queue(speed_file_names[0], queue_true_state_file, queue_reg_len, grid)
    # plot_speed_with_queue(speed_file_names[opt_temp_config_num - 1], queue_true_state_file, queue_reg_len, grid)
    # plot_speed_with_queue(speed_file_names[-1], queue_true_state_file, queue_reg_len, grid)

# os.remove(log_dir + 'temporary_I80_configuration_inputs.txt')

# cross_eval.plot_true_speed_field_partition(grid, replication[0])
#
plt.show()


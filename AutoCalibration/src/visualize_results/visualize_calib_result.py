__author__ = 'Yanning'
"""
This script visualizes the estimation results
"""


from post_calibration_class import *
import matplotlib.pyplot as plt

# ========================================================================== #
# Visualize the simulated speeds.
# ========================================================================== #
# file_path = '/Users/Yanning/Dropbox/DesktopServer/Workzone_autocalibration/I57_SB/' \
#             'Logs/thread2/'
# time_stamp = '20160708_224007'
#
# obj_ratio = [1.0, 0.0]
#
# compare_sol = ['valid', 'default', 'newseed']
#
# pc = Post_Calibration()
#
# for sol in compare_sol:
#     pc.read_sim_data(file_path, time_stamp + '_sol_' + sol + '.csv', sol)
#
# ## for I00_EB
# # pc.plot_calib_speed_data(compare_sol, ['EB1', 'EB2', 'EB3', 'EB4'])
# # pc.plot_calib_count_data(compare_sol, ['EB1', 'EB2', 'EB3', 'EB4'])
#
#
# ##  for I80_EB
# # pc.plot_calib_speed_data(compare_sol, ['EB4', 'EB5', 'EB6', 'EB7'])
# # pc.plot_calib_speed_data(compare_sol, ['EB8', 'EB9', 'EB10', 'EB12'])
# # pc.plot_calib_speed_data(compare_sol, ['EB12', 'EB14', 'EB15', 'EB16'])
#
#
# ## for I57_SB
# pc.plot_calib_speed_data(compare_sol, ['SB1', 'SB2', 'SB3', 'SB4'])
# pc.plot_calib_speed_data(compare_sol, ['SB3', 'SB4', 'SB5', 'SB6'])
# pc.plot_calib_speed_data(compare_sol, ['SB6', 'SB7', 'SB8', 'SB9'])
#
# ## visualize the optimization steps
# pc.read_opt_steps(file_path + time_stamp + '_opt_log.csv')
# pc.plot_opt_steps(obj_ratio)


# ========================================================================== #
# Visualize the sensitivity analysis result
# ========================================================================== #
# sensitivity_analysis_list = ['car_minDist', 'car_sensitivityFactor', 'car_reactionTime', 'car_maxAccel',
#                              'truck_minDist', 'truck_sensitivityFactor', 'truck_reactionTime', 'truck_maxAccel',
#                              'car_speedAcceptance', 'truck_speedAcceptance']
# pc.read_SA_steps(file_path, folder_name, time_stamp)
# pc.plot_SA_result(sensitivity_analysis_list[0:4])
# pc.plot_SA_result(sensitivity_analysis_list[4:8])
# pc.plot_SA_result(sensitivity_analysis_list[8:10])


# pc.plot_single_sensitivity_analysis_result(sensitivity_analysis_list[5])


# ========================================================================== #
# Visualize the time space diagram
# ========================================================================== #
# this is the time space timegram
# det_list = ['EB4', 'EB5', 'EB6', 'EB7', 'EB8', 'EB9', 'EB10', 'EB12', 'EB14', 'EB15', 'EB16']
# det_list = ['SB1', 'SB2', 'SB3', 'SB4', 'SB5', 'SB6', 'SB7', 'SB8', 'SB9']
# pc.plot_calib_data_time_space_diagram(['optimal', 'valid'], det_list, [], [])

# plot I57 FD
# pc.plot_I57_SB1_FD('valid')




# ========================================================================== #
# Analysis and visualization of the response surface
# ========================================================================== #
# result_file = '20160512_231055_visualization_step_log.csv'
# solution_file = '20160513_191955_opt_log.csv'
# mat_size = (46,31)
#
# # parameters for plot the high dimensional parallel coordinates
# paras_to_plot = ['car_maxAccel', 'car_reactionTime', 'car_minDist', 'car_sensitivityFactor',
#               'truck_maxAccel', 'truck_reactionTime', 'truck_minDist', 'truck_sensitivityFactor']
# para_bounds = {}
# para_bounds['car_maxAccel'] = [2.6, 3.4]
# para_bounds['car_reactionTime'] = [0.5, 1.2]
# para_bounds['car_minDist'] = [0.5, 1.5]
# para_bounds['car_sensitivityFactor'] = [0.9, 1.1]
# para_bounds['truck_maxAccel'] = [0.6, 1.8]
# para_bounds['truck_reactionTime'] = [0.5, 1.2]
# para_bounds['truck_minDist'] = [1.0, 2.5]
# para_bounds['truck_sensitivityFactor'] = [0.9, 1.1]
#
# # true parameters
# true_paras = {}
# true_paras['car_maxAccel'] = [3.2, 0.2, 2.6, 3.4]   # range: [2.6, 3.4], # 81
# true_paras['car_reactionTime'] = [0.7, 1.2, 1.6, 1] # range: [0.5, 1.2], # 8
# true_paras['car_minDist'] = [0.85, 0.3, 0.5, 1.5]   # range: [0.5, 1.5], # 101
# true_paras['car_sensitivityFactor'] = [1.05, 0, 1.05, 1.05]     # range: [0.9, 1.1], # 21
# true_paras['truck_maxAccel'] = [0.9, 0.5, 0.6, 1.8]     # range: [0.6, 1.8], # 121
# true_paras['truck_reactionTime'] = [0.9, 1.3, 1.7, 1]   # range: [0.5, 1.2], # 8
# true_paras['truck_minDist'] = [1.34, 0.5, 1.0, 2.5]     # range: [1.0, 2.5], # 151
# true_paras['truck_sensitivityFactor'] = [1.03, 0, 1.03, 1.03]   # range: [0.9, 1.1], # 21
#
# true_paras_val = [true_paras, 0]
#
# plot_simulation_mesh(result_file, solution_file, mat_size)


# ========================================================================== #
# Visualization of the calibration steps in a high dimensional parallel plot
# ========================================================================== #
# solution_file = '/Users/Yanning/Dropbox/DesktopServer/Workzone_autocalibration/I80_EB/' \
#                 'Logs/Calib2/Calib2_round3/20160706_094240_opt_log.csv'
# solution_file = '/Users/Yanning/Dropbox/DesktopServer/Workzone_autocalibration/I57_SB/' \
#                 'Logs/I57_800_opt_step.csv'
# parameters for plot the high dimensional parallel coordinates
# paras_to_plot = ['car_speedAcceptance','truck_maxAccel',
#                  'car_sensitivityFactor', 'truck_sensitivityFactor',
#                  'car_reactionTime', 'truck_reactionTime',
#                  'car_minHeadway', 'truck_minHeadway']
#
# para_bounds = {}
# para_bounds['car_speedAcceptance'] = [1.05, 1.3]
# para_bounds['truck_maxAccel'] = [0.6, 1.8]
# para_bounds['car_sensitivityFactor'] = [0.5, 1.5]
# para_bounds['truck_sensitivityFactor'] = [0.5, 1.5]
# para_bounds['car_reactionTime'] = [0.6, 1.6]
# para_bounds['truck_reactionTime'] = [0.6, 1.6]
# para_bounds['car_minHeadway'] = [1, 4.5]
# para_bounds['truck_minHeadway'] = [1, 4.5]


# true parameters
# true_paras_val = None
# #
# plot_high_dimension(solution_file, paras_to_plot, true_paras_val, para_bounds)


# ========================================================================== #
# find out the optimal value in opt step file with additional constraints
# original
# para_bounds = {}
# para_bounds['car_speedAcceptance'] = [0.85, 1.3]
# para_bounds['truck_maxAccel'] = [0.6, 1.8]
# para_bounds['car_sensitivityFactor'] = [0.5, 1.5]
# para_bounds['truck_sensitivityFactor'] = [0.5, 1.5]
#
# para_bounds['car_reactionTime'] = [0.6, 1.6]
# para_bounds['truck_reactionTime'] = [0.6, 1.6]
# para_bounds['car_minHeadway'] = [1, 4.5]
# para_bounds['truck_minHeadway'] = [1, 4.5]


"""
Form I80
The final calibrated values
=============== Calib 1 ==============
Total number of unique steps 1058/1060
Total number of feasible steps 294/1060
Optimal parameters with obj: 168.3689373
--- car_speedAcceptance:[1.09, 0.1, 0.85, 1.1]
--- truck_maxAccel:[0.82, 0.5, 0.6, 1.8]
--- car_sensitivityFactor:[0.5, 0.0, 0.5, 0.5]
--- truck_sensitivityFactor:[1.02, 0.0, 1.02, 1.02]
--- car_reactionTime:[0.6, 1.2, 1.6, 1.0]
--- truck_reactionTime:[0.8, 1.3, 1.7, 1.0]
--- car_minHeadway:[1.5, 0.0, 1.5, 1.5]
--- truck_minHeadway:[2.5, 0.0, 2.5, 2.5]

=============== Calib 2 ==============
Total number of unique steps 1322/1460
Total number of feasible steps 234/1460
Optimal parameters with obj: 178.581337235
--- car_speedAcceptance:[1.05, 0.1, 0.85, 1.1]
--- truck_maxAccel:[0.66, 0.5, 0.6, 1.8]
--- car_sensitivityFactor:[0.5, 0.0, 0.5, 0.5]
--- truck_sensitivityFactor:[0.61, 0.0, 0.61, 0.61]
--- car_reactionTime:[1.0, 1.2, 1.6, 1.0]
--- truck_reactionTime:[0.8, 1.3, 1.7, 1.0]
--- car_minHeadway:[1.0, 0.0, 1.0, 1.0]
--- truck_minHeadway:[2.0, 0.0, 2.0, 2.0]
"""


"""
For I57, optimal parameters
Total number of unique steps 791/794
Total number of feasible steps 135/794
Optimal parameters with obj: 253.74 (default: 377)
--- car_speedAcceptance:[1.1, 0.1, 0.85, 1.1]
--- truck_maxAccel:[0.85, 0.5, 0.6, 1.8]
--- car_sensitivityFactor:[0.61, 0.0, 0.61, 0.61]
--- truck_sensitivityFactor:[0.87, 0.0, 0.87, 0.87]
--- car_reactionTime:[0.6, 1.2, 1.6, 1.0]
--- truck_reactionTime:[0.6, 1.3, 1.7, 1.0]
--- car_minHeadway:[3.0, 0.0, 3.0, 3.0]
--- truck_minHeadway:[4.5, 0.0, 4.5, 4.5]
"""

para_bounds = {}
para_bounds['car_speedAcceptance'] = [0.85, 1.3]
para_bounds['truck_maxAccel'] = [0.6, 1.8]
para_bounds['car_sensitivityFactor'] = [0.5, 1.5]
para_bounds['truck_sensitivityFactor'] = [0.5, 1.5]

para_bounds['car_reactionTime'] = [0.6, 1.6]
para_bounds['truck_reactionTime'] = [0.6, 1.6]
para_bounds['car_minHeadway'] = [1.5, 2.0]
para_bounds['truck_minHeadway'] = [1.5, 2.0]

# print('=============== Calib 2 ==============')
# solution_file = '/Users/Yanning/Dropbox/DesktopServer/Workzone_autocalibration/I80_EB/' \
#                 'Logs/Calib1/Calib1_opt_log.csv'
solution_file = '/Users/Yanning/Dropbox/DesktopServer/Workzone_autocalibration/I57_SB/' \
                'Logs/I57_800_opt_step.csv'

find_feasible_optimum(solution_file, para_bounds)



# ========================================================================== #
# Visualizing the sorted opt steps to as the stopping criteria
# ========================================================================== #

# pc = Post_Calibration()
# obj_ratio = [1.0, 0.0]
#
# pc.read_opt_steps('I57_800_opt_step.csv')
# pc.plot_opt_steps(obj_ratio)
# pc.plot_sorted_opt_steps(obj_ratio)
#
#
# plt.show()








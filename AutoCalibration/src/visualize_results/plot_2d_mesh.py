_author__ = 'Yanning Li'
"""
This script is used to plot the 2d mesh for the visual_2d result.
"""

import numpy as np
import sys
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.mplot3d import Axes3D



def main(argv):

    result_file = '20160512_231055_visualization_step_log.csv'
    solution_file = '20160513_191955_opt_log.csv'
    mat_size = (46,31)

    # parameters for plot the high dimensional parallel coordinates
    paras_to_plot = ['car_maxAccel', 'car_reactionTime', 'car_minDist', 'car_sensitivityFactor',
                  'truck_maxAccel', 'truck_reactionTime', 'truck_minDist', 'truck_sensitivityFactor']
    para_bounds = {}
    para_bounds['car_maxAccel'] = [2.6, 3.4]
    para_bounds['car_reactionTime'] = [0.5, 1.2]
    para_bounds['car_minDist'] = [0.5, 1.5]
    para_bounds['car_sensitivityFactor'] = [0.9, 1.1]
    para_bounds['truck_maxAccel'] = [0.6, 1.8]
    para_bounds['truck_reactionTime'] = [0.5, 1.2]
    para_bounds['truck_minDist'] = [1.0, 2.5]
    para_bounds['truck_sensitivityFactor'] = [0.9, 1.1]

    # true parameters
    true_paras = {}
    true_paras['car_maxAccel'] = [3.2, 0.2, 2.6, 3.4]   # range: [2.6, 3.4], # 81
    true_paras['car_reactionTime'] = [0.7, 1.2, 1.6, 1] # range: [0.5, 1.2], # 8
    true_paras['car_minDist'] = [0.85, 0.3, 0.5, 1.5]   # range: [0.5, 1.5], # 101
    true_paras['car_sensitivityFactor'] = [1.05, 0, 1.05, 1.05]     # range: [0.9, 1.1], # 21
    true_paras['truck_maxAccel'] = [0.9, 0.5, 0.6, 1.8]     # range: [0.6, 1.8], # 121
    true_paras['truck_reactionTime'] = [0.9, 1.3, 1.7, 1]   # range: [0.5, 1.2], # 8
    true_paras['truck_minDist'] = [1.34, 0.5, 1.0, 2.5]     # range: [1.0, 2.5], # 151
    true_paras['truck_sensitivityFactor'] = [1.03, 0, 1.03, 1.03]   # range: [0.9, 1.1], # 21

    true_paras_val = [true_paras, 0]

    plot_high_dimension(solution_file, paras_to_plot, true_paras_val, para_bounds)
    # plot_simulation_mesh(result_file, solution_file, mat_size)





if __name__ == "__main__":
    sys.exit(main(sys.argv))












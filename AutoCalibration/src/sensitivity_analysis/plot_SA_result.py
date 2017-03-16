__author__ = 'Yanning Li'
"""
This script plots the sensitivity analysis result for all parameters under investigation.

"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict


def main(argv):

    SA_file = 'SA_12_factors.csv'

    sensitivity_analysis_list = ['car_speedAcceptance', 'car_maxAccel',
                                 'truck_speedAcceptance', 'truck_maxAccel',
                                 'car_minDist', 'car_sensitivityFactor',
                                 'truck_minDist', 'truck_sensitivityFactor',
                                 'car_reactionTime', 'truck_reactionTime',
                                 'car_cov', 'truck_cov']
    SA_result = read_SA_file(SA_file)

    subplots_SA(sensitivity_analysis_list[0:4], SA_result)
    subplots_SA(sensitivity_analysis_list[4:8], SA_result)
    subplots_SA(sensitivity_analysis_list[8:12], SA_result)


    # plot_SA_for_para(sensitivity_analysis_list[5], SA_result)

    plt.show()


def read_SA_file(file_name):
    """
    This runction reads the logged data for the SA
    :param file_name: str, the full file name of the SA_result file
    :return: SA_result. a dict. keys are the parameters, values are two lists: the key and response values
             SA_result['car_maxAccel'] = [ [key_values], [response_values]  ]
    """

    f = open(file_name, 'r')

    SA_result = OrderedDict()

    for line in f:

        val = []
        obj = []

        line = line.strip()
        items = line.split(',')

        name = items[0].split('_')
        key_name = name[0] + '_' + name[1]
        if key_name not in SA_result.keys():
            SA_result[key_name] = [[],[]]

        if name[2] == 'val':
            for i in range(1, len(items)):
                val.append(float(items[i]))
            SA_result[key_name][0] = val
        elif name[2] == 'obj':
            for i in range(1, len(items)):
                obj.append((float(items[i])))
            SA_result[key_name][1] = obj
        else:
            print 'unrecognized entry {0}'.format(items[0])

    f.close()

    return SA_result



# plot the sensitivity analysis result
def subplots_SA(para_name_list, SA_result):
    """
    This script plots the SA result in a window with four subfigures
    :param para_name_list: A list of strings. Each string is the parameter name to be plotted
    :param SA_result: The SA_result output from read_SA_file
    :return: figures
    """

    # speed
    # Four axes, returned as a 2-d array
    f, axarr = plt.subplots(2, 2, figsize = (15,10))

    if len(para_name_list) >= 1:
        para_name = para_name_list[0]
        x = SA_result[para_name][0]
        y = SA_result[para_name][1]

        axarr[0, 0].plot(x, y, '*b')
        # compute the mean
        [mean_x, mean_y] = compute_unique_mean(x, y)
        axarr[0, 0].plot(mean_x, mean_y, 'r', linewidth=2)

        axarr[0, 0].set_title('Parameter {0}'.format(para_name))
        axarr[0, 0].set_xlabel('{0} value'.format(para_name))
        axarr[0, 0].set_ylabel('RMSE of speed')
        axarr[0, 0].set_xlim([ np.min(x) - 0.1*(np.max(x)-np.min(x)) ,
                               np.max(x) + 0.1*(np.max(x)-np.min(x))])
        axarr[0, 0].grid(True)

    if len(para_name_list) >= 2:
        para_name = para_name_list[1]
        x = SA_result[para_name][0]
        y = SA_result[para_name][1]

        axarr[0, 1].plot(x, y, '*')
        # compute the mean
        [mean_x, mean_y] = compute_unique_mean(x, y)
        axarr[0, 1].plot(mean_x, mean_y, 'r', linewidth=2)

        axarr[0, 1].set_title('Parameter {0}'.format(para_name))
        axarr[0, 1].set_xlabel('{0} value'.format(para_name))
        axarr[0, 1].set_ylabel('RMSE of speed')
        axarr[0, 1].set_xlim([ np.min(x) - 0.1*(np.max(x)-np.min(x)) ,
                               np.max(x) + 0.1*(np.max(x)-np.min(x))])
        axarr[0, 1].grid(True)

    if len(para_name_list) >= 3:
        para_name = para_name_list[2]
        x = SA_result[para_name][0]
        y = SA_result[para_name][1]

        axarr[1, 0].plot(x, y, '*')
        # compute the mean
        [mean_x, mean_y] = compute_unique_mean(x, y)
        axarr[1, 0].plot(mean_x, mean_y, 'r', linewidth=2)

        axarr[1, 0].set_title('Parameter {0}'.format(para_name))
        axarr[1, 0].set_xlabel('{0} value'.format(para_name))
        axarr[1, 0].set_ylabel('RMSE of speed')
        axarr[1, 0].set_xlim([ np.min(x) - 0.1*(np.max(x)-np.min(x)) ,
                               np.max(x) + 0.1*(np.max(x)-np.min(x))])
        axarr[1, 0].grid(True)

    if len(para_name_list) == 4:
        para_name = para_name_list[3]
        x = SA_result[para_name][0]
        y = SA_result[para_name][1]

        axarr[1, 1].plot(x, y, '*')
        # compute the mean
        [mean_x, mean_y] = compute_unique_mean(x, y)
        axarr[1, 1].plot(mean_x, mean_y, 'r', linewidth=2)

        axarr[1, 1].set_title('Parameter {0}'.format(para_name))
        axarr[1, 1].set_xlabel('{0} value'.format(para_name))
        axarr[1, 1].set_ylabel('RMSE of speed')
        axarr[1, 1].set_xlim([ np.min(x) - 0.1*(np.max(x)-np.min(x)) ,
                               np.max(x) + 0.1*(np.max(x)-np.min(x))])
        axarr[1, 1].grid(True)

    if len(para_name_list) > 4:
        print 'Only support up to 4 parameters in one figure window'

    plt.draw()



# plot the sensitivity analysis result
def plot_SA_for_para(para_name, SA_result):
    """
    This script plots the result for para_name in a single figure
    :param para_name: str, the name of the parameter to be plotted
    :param SA_result: the output of read_SA_result
    :return:
    """

    # speed
    fig_window = plt.figure(figsize=(15,10))
    fig = fig_window.add_subplot(111)


    x = SA_result[para_name][0]
    y = SA_result[para_name][1]

    plt.plot(x, y, '*b')
    # compute the mean
    [mean_x, mean_y] = compute_unique_mean(x, y)
    plt.plot(mean_x, mean_y, 'r', linewidth=2)

    plt.title('Parameter Truck sensitivity factor', fontsize=24)
    plt.xlabel('Sensitivity factor value for trucks', fontsize=24)
    plt.ylabel('RMSE of speed', fontsize=24)
    # plt.xlim([ np.min(x) - 0.1*(np.max(x)-np.min(x)) , np.max(x) + 0.1*(np.max(x)-np.min(x))])
    plt.xlim([ np.min(x) - 0.1*(np.max(x)-np.min(x)),
               np.max(x) + 0.1*(np.max(x)-np.min(x)) ])
    plt.grid(True)
    fig.tick_params(axis='both', which='major', labelsize=20)


    plt.draw()





# find the unique values of parameter, and compute their mean objective value
# input: x: 1,1,1,1,1,2,2,2,2,2,....
#        y: 4,4,4,4,4,6,6,6,6,6,
# output: p = 1,2
#         o = 4,6
def compute_unique_mean(x,y):

    x_a = np.array(x)
    y_a = np.array(y)

    p_unique = np.unique(x_a)

    p_unique = np.sort(p_unique)

    para = []
    obj = []
    for p_value in p_unique:

        para.append(p_value)

        tmp_list = []
        for i in range(0, len(x_a)):
            if x_a[i] == p_value:
                tmp_list.append(y_a[i])

        obj.append(np.mean(tmp_list))

    return [para, obj]




if __name__ == "__main__":
    sys.exit(main(sys.argv))

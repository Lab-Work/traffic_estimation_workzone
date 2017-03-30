__author__ = 'Yanning Li'
"""
This is the post calibration class and static functions.
    - It visualizes the calibration results, e.g., the comparison of the simulated speed and the
    measured speed, using the default parameters, the optimal parameters, etc.
    - It visualizes the high-dimension parallel plot of the optimization steps.
    - It visualizes the sensitivity analysis results.
"""




import csv
import sys
import numpy as np
from matplotlib import cm
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D
from collections import OrderedDict
from copy import deepcopy



class Post_Calibration:

    def __init__(self):

        # calibration data is a dict; keys are the names
        # Two reserved names: valid, default
        # calib_data[sol_name] = [paras, sim_data]
        # paras[para_name] = [values]
        # sim_data[det_name] = [speed, count]
        self.calib_data = OrderedDict()

        # optimization steps
        self.opt_steps = None

        # save the sensitivity analysis result
        # sensitivity_reault[para_name] = [[val],[obj]]
        self.sensitivity_result = OrderedDict()


    # read optimization steps
    def read_opt_steps(self,  file_name):

        f = open(file_name, 'r')
        reader = csv.reader(f)

        tmp_list = list()
        for row in reader:
            tmp_list.append(row)

        # save in numpy array
        self.opt_steps = np.array(tmp_list)

        f.close()


    # read a data set from a file
    # if name is valid: then set valid_data;
    # if name is default: then set default data;
    # otherwise assign name
    def read_sim_data(self, file_path, file_name, name):

        f = open(file_path+file_name, 'r')

        # first line is the parameters
        line = f.readline()
        tup = line.split(',')
        # para = OrderedDict()
        # para['car'] = [tup[0], tup[1]]
        # para['truck'] = [tup[2], tup[3]]
        # TODO: changed file format; to fix reading paras
        para = None

        # read the rest data
        data = OrderedDict()
        while True:
            line = f.readline()
            # if EOF
            if not (bool(line)):
                break
            line = line.strip()
            items = line.split(',')

            if not items[0] in data.keys():
                # read speed;
                speed = items[1:]

                data[str(items[0])] = [speed, None]
            else:
                # save count data
                count = items[1:]
                data[str(items[0])][1] = count

        # name is the solution name
        # para is the set of parameters,
        # data[det_name] = [speed, count]
        self.calib_data[name] = (para, data)



    # plot the simulation steps along with the objective function value
    def plot_opt_steps(self, obj_ratio):

        f, ax = plt.subplots(figsize=(15, 8))

        obj_speed = (self.opt_steps[:,-2]).astype(np.float)
        obj_count = (self.opt_steps[:,-1]).astype(np.float)
        #print obj_speed
        #print obj_count
        obj_val = obj_ratio[0]*obj_speed + obj_ratio[1]* obj_count

        # iter_num = self.opt_steps[:,0]
        iter_num = np.arange(0, len(obj_val))

        ax.plot(iter_num, obj_speed, label = 'speed')
        # ax.plot(iter_num, obj_count, label = 'count')
        ax.plot(iter_num, obj_val, label = '{0}xS+{1}xC'.format(obj_ratio[0], obj_ratio[1]), linewidth=2)
        ax.grid(True)
        ax.set_title('Optimization steps', fontsize=24)
        ax.set_xlabel('Iteration step', fontsize=20)
        ax.set_ylabel('Objective function value', fontsize=20)
        ax.set_xlim([0, len(obj_val)])

        plt.legend( loc='upper right' )

        # plt.savefig('opt_steps.png')

        plt.draw()


    def plot_sorted_opt_steps(self, obj_ratio):
        """
        This function plots the sorted optimization objective values as a time series.
        :param obj_ratio: RMSE_speed*ratio[0] + RMSE_flow*ratio[1]
        :return:
        """

        f, ax = plt.subplots(figsize=(15, 8))

        obj_speed = (self.opt_steps[:,-2]).astype(np.float)
        obj_count = (self.opt_steps[:,-1]).astype(np.float)
        #print obj_speed
        #print obj_count
        obj_val = obj_ratio[0]*obj_speed + obj_ratio[1]* obj_count

        # the current optimal
        optimum = np.inf
        opt_val = []
        for val in obj_val:
            if val <= optimum:
                optimum = val
            opt_val.append(optimum)

        # sort the objective function value
        sorted_obj_val = sorted(obj_val, reverse=True)

        # iter_num = self.opt_steps[:,0]
        iter_num = np.arange(0, len(obj_val))

        # plot the unsorted values
        ax.plot(iter_num, obj_val, label = 'unsorted', linewidth=2)
        # ax.plot(iter_num, obj_count, label = 'count')
        ax.plot(iter_num, sorted_obj_val, label = 'sorted', linewidth=2)
        ax.plot(iter_num, opt_val, label='optimum', linewidth=2)
        ax.grid(True)
        ax.set_title('Optimization steps', fontsize=24)
        ax.set_xlabel('Iteration step', fontsize=20)
        ax.set_ylabel('Objective function value', fontsize=20)
        ax.set_xlim([0, len(obj_val)])

        plt.legend( loc='upper right' )

        # plt.savefig('opt_steps.png')

        plt.draw()


    # compare the data
    # input: compare_name_list is the name list for those solutions to be compared
    #        det_name_list is the detector names that to be plotted with up to four detectors (1,2,4)
    #                      if put [] there, will plot one detector in each window
    # We only plot the speed, since the count data is not reliable
    def plot_calib_speed_data(self, compare_name_list, det_name_list):

        # One figure window for speed with four subplots (one detector each)

        # speed
        # Four axes, returned as a 2-d array
        if len(det_name_list) == 4:
            f, axarr = plt.subplots(2, 2, figsize = (15,10))

            for name in compare_name_list:
                y = self.calib_data[name][1][det_name_list[0]][0]
                axarr[0, 0].plot(y, label = name)
            axarr[0, 0].set_title('Detector {0}'.format(det_name_list[0]), fontsize=25)
            axarr[0, 0].set_ylabel('Speed (mph)', fontsize=20)
            axarr[0, 0].grid(True)

            for name in compare_name_list:
                y = self.calib_data[name][1][det_name_list[1]][0]
                axarr[0, 1].plot(y, label = name)
            axarr[0, 1].set_title('Detector {0}'.format(det_name_list[1]), fontsize=25)
            axarr[0, 1].grid(True)

            for name in compare_name_list:
                y = self.calib_data[name][1][det_name_list[2]][0]
                axarr[1, 0].plot(y, label = name)
            axarr[1, 0].set_title('Detector {0}'.format(det_name_list[2]), fontsize=25)
            axarr[1, 0].grid(True)
            axarr[1, 0].set_ylabel('Speed (mph)', fontsize=20)
            axarr[1, 0].set_xlabel('Time (aggregated every 5 min)', fontsize=20)

            for name in compare_name_list:
                y = self.calib_data[name][1][det_name_list[3]][0]
                axarr[1, 1].plot(y, label = name)
            axarr[1, 1].set_title('Detector {0}'.format(det_name_list[3]), fontsize=25)
            axarr[1, 1].grid(True)
            axarr[1, 1].set_xlabel('Time (aggregated every 5 min)', fontsize=20)

            # Fine-tune figure; hide x ticks for top plots and y ticks for right plots
            # plt.setp([a.get_xticklabels() for a in axarr[0, :]], visible=False)
            # plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False)

            # plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
            plt.legend(loc='best')
            # plt.savefig('speed.png')

        if len(det_name_list) == 3:
            f, axarr = plt.subplots(2, 2, figsize = (15,10))

            for name in compare_name_list:
                y = self.calib_data[name][1][det_name_list[0]][0]
                axarr[0, 0].plot(y, label = name)
            axarr[0, 0].set_title('Detector {0}'.format(det_name_list[0]), fontsize=25)
            axarr[0, 0].set_ylabel('Speed (mph)', fontsize=20)
            axarr[0, 0].grid(True)

            # for name in compare_name_list:
            #     y = self.calib_data[name][1][det_name_list[1]][0]
            #     axarr[0, 1].plot(y, label = name)
            # axarr[0, 1].set_title('Detector {0}'.format(det_name_list[1]))
            # axarr[0, 1].grid(True)

            for name in compare_name_list:
                y = self.calib_data[name][1][det_name_list[1]][0]
                axarr[1, 0].plot(y, label = name)
            axarr[1, 0].set_title('Detector {0}'.format(det_name_list[1]), fontsize=25)
            axarr[1, 0].grid(True)
            axarr[1, 0].set_ylabel('Speed (mph)', fontsize=20)
            axarr[1, 0].set_xlabel('Time (aggregated every 5 min)', fontsize=20)

            for name in compare_name_list:
                y = self.calib_data[name][1][det_name_list[2]][0]
                axarr[1, 1].plot(y, label = name)
            axarr[1, 1].set_title('Detector {0}'.format(det_name_list[2]), fontsize=25)
            axarr[1, 1].grid(True)
            axarr[1, 1].set_xlabel('Time (aggregated every 5 min)', fontsize=20)

            # Fine-tune figure; hide x ticks for top plots and y ticks for right plots
            # plt.setp([a.get_xticklabels() for a in axarr[0, :]], visible=False)
            # plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False)

            # plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
            plt.legend( loc='best' )
            # plt.savefig('speed.png')

        elif len(det_name_list) == 2:
            f, axarr = plt.subplots(2, 1, figsize = (10,10))

            for name in compare_name_list:
                y = self.calib_data[name][1][det_name_list[0]][0]
                axarr[0].plot(y, label = name)
            axarr[0].set_title('Detector {0}'.format(det_name_list[0]), fontsize=25)
            axarr[0].set_ylabel('Speed (mph)', fontsize=20)
            axarr[0].grid(True)

            for name in compare_name_list:
                y = self.calib_data[name][1][det_name_list[1]][0]
                axarr[1].plot(y, label = name)
            axarr[1].set_title('Detector {0}'.format(det_name_list[1]), fontsize=25)
            axarr[1].set_ylabel('Speed (mph)', fontsize=20)
            axarr[1].set_xlabel('Time (aggregated every 5 min)', fontsize=20)
            axarr[1].grid(True)

            plt.legend( loc='best' )

        elif len(det_name_list) == 1:
            f, ax = plt.subplots(figsize = (10,10))

            for name in compare_name_list:
                y = self.calib_data[name][1][det_name_list[0]][0]
                ax.plot(y, label = name)
            ax.set_title('Detector {0}'.format(det_name_list[0]), fontsize=25)
            ax.set_ylabel('Speed (mph)', fontsize=20)
            ax.grid(True)

            plt.legend( loc='best' )

        elif len(det_name_list) == 0:
            pass
            # TODO: to be impolemented: if 0, then plot one detector in each window

        plt.draw()

        # count
        # Four axes, returned as a 2-d array
        # f1, axarr1 = plt.subplots(2, 2, figsize = (15,10))
        # for name in compare_name_list:
        #     y = self.calib_data[name][1]['EB1'][1]
        #     axarr1[0, 0].plot(y, label = name)
        # axarr1[0, 0].set_title('count EB1')
        # axarr1[0, 0].grid(True)
        # for name in compare_name_list:
        #     y = self.calib_data[name][1]['EB2'][1]
        #     axarr1[0, 1].plot(y, label = name)
        # axarr1[0, 1].set_title('count EB2')
        # axarr1[0, 1].grid(True)
        # for name in compare_name_list:
        #     y = self.calib_data[name][1]['EB3'][1]
        #     axarr1[1, 0].plot(y, label = name)
        # axarr1[1, 0].set_title('count EB3')
        # axarr1[1, 0].grid(True)
        # axarr1[1, 0].set_xlabel('Time')
        # for name in compare_name_list:
        #     y = self.calib_data[name][1]['EB4'][1]
        #     axarr1[1, 1].plot(y, label = name)
        # axarr1[1, 1].set_title('count EB4')
        # axarr1[1, 1].grid(True)
        # axarr1[1, 1].set_xlabel('Time')
        # # Fine-tune figure; hide x ticks for top plots and y ticks for right plots
        # plt.setp([a.get_xticklabels() for a in axarr1[0, :]], visible=False)
        # plt.setp([a.get_yticklabels() for a in axarr1[:, 1]], visible=False)
        #
        # plt.legend( loc='upper right' )
        # plt.savefig('count.png')




    # compare the data
    # input: compare_name_list is the name list for those solutions to be compared
    #        det_name_list is the detector names that to be plotted with up to four detectors (1,2,4)
    #                      if put [] there, will plot one detector in each window
    # We only plot the speed, since the count data is not reliable
    def plot_calib_count_data(self, compare_name_list, det_name_list):

        # One figure window for speed with four subplots (one detector each)

        # speed
        # Four axes, returned as a 2-d array
        if len(det_name_list) == 4:
            f, axarr = plt.subplots(2, 2, figsize = (15,10))

            for name in compare_name_list:
                y = self.calib_data[name][1][det_name_list[0]][1]
                axarr[0, 0].plot(y, label = name)
            axarr[0, 0].set_title('Detector {0}'.format(det_name_list[0]))
            axarr[0, 0].set_ylabel('Count (/5min)')
            axarr[0, 0].grid(True)

            for name in compare_name_list:
                y = self.calib_data[name][1][det_name_list[1]][1]
                axarr[0, 1].plot(y, label = name)
            axarr[0, 1].set_title('Detector {0}'.format(det_name_list[1]))
            axarr[0, 1].grid(True)

            for name in compare_name_list:
                y = self.calib_data[name][1][det_name_list[2]][1]
                axarr[1, 0].plot(y, label = name)
            axarr[1, 0].set_title('Detector {0}'.format(det_name_list[2]))
            axarr[1, 0].grid(True)
            axarr[1, 0].set_ylabel('Count (/5min)')
            axarr[1, 0].set_xlabel('Time (aggregated every 5 min)')

            for name in compare_name_list:
                y = self.calib_data[name][1][det_name_list[3]][1]
                axarr[1, 1].plot(y, label = name)
            axarr[1, 1].set_title('Detector {0}'.format(det_name_list[3]))
            axarr[1, 1].grid(True)
            axarr[1, 1].set_xlabel('Time (aggregated every 5 min)')

            # Fine-tune figure; hide x ticks for top plots and y ticks for right plots
            # plt.setp([a.get_xticklabels() for a in axarr[0, :]], visible=False)
            # plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False)

            # plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
            plt.legend( loc='upper right' )
            # plt.savefig('speed.png')

        if len(det_name_list) == 3:
            f, axarr = plt.subplots(2, 2, figsize = (15,10))

            for name in compare_name_list:
                y = self.calib_data[name][1][det_name_list[0]][1]
                axarr[0, 0].plot(y, label = name)
            axarr[0, 0].set_title('Detector {0}'.format(det_name_list[0]))
            axarr[0, 0].set_ylabel('Count (/5min)')
            axarr[0, 0].grid(True)

            # for name in compare_name_list:
            #     y = self.calib_data[name][1][det_name_list[1]][0]
            #     axarr[0, 1].plot(y, label = name)
            # axarr[0, 1].set_title('Detector {0}'.format(det_name_list[1]))
            # axarr[0, 1].grid(True)

            for name in compare_name_list:
                y = self.calib_data[name][1][det_name_list[1]][1]
                axarr[1, 0].plot(y, label = name)
            axarr[1, 0].set_title('Detector {0}'.format(det_name_list[1]))
            axarr[1, 0].grid(True)
            axarr[1, 0].set_ylabel('Count (/5min)')
            axarr[1, 0].set_xlabel('Time (aggregated every 5 min)')

            for name in compare_name_list:
                y = self.calib_data[name][1][det_name_list[2]][1]
                axarr[1, 1].plot(y, label = name)
            axarr[1, 1].set_title('Detector {0}'.format(det_name_list[2]))
            axarr[1, 1].grid(True)
            axarr[1, 1].set_xlabel('Time (aggregated every 5 min)')

            # Fine-tune figure; hide x ticks for top plots and y ticks for right plots
            # plt.setp([a.get_xticklabels() for a in axarr[0, :]], visible=False)
            # plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False)

            # plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
            plt.legend( loc='upper right' )
            # plt.savefig('speed.png')

        elif len(det_name_list) == 2:
            f, axarr = plt.subplots(2, 1, figsize = (10,10))

            for name in compare_name_list:
                y = self.calib_data[name][1][det_name_list[0]][1]
                axarr[0].plot(y, label = name)
            axarr[0].set_title('Detector {0}'.format(det_name_list[0]))
            axarr[0].set_ylabel('Speed (mph)')
            axarr[0].grid(True)

            for name in compare_name_list:
                y = self.calib_data[name][1][det_name_list[1]][1]
                axarr[1].plot(y, label = name)
            axarr[1].set_title('Detector {0}'.format(det_name_list[1]))
            axarr[1].set_ylabel('Speed (mph)')
            axarr[1].set_xlabel('Time (aggregated every 5 min)')
            axarr[1].grid(True)

            plt.legend( loc='upper right' )

        elif len(det_name_list) == 1:
            f, ax = plt.subplots(figsize = (10,10))

            for name in compare_name_list:
                y = self.calib_data[name][1][det_name_list[0]][1]
                ax.plot(y, label = name)
            ax.set_title('Detector {0}'.format(det_name_list[0]))
            ax.set_ylabel('Speed (mph)')
            ax.grid(True)

            plt.legend( loc='upper right' )

        elif len(det_name_list) == 0:
            pass
            # TODO: to be impolemented: if 0, then plot one detector in each window

        plt.draw()

        # count
        # Four axes, returned as a 2-d array
        # f1, axarr1 = plt.subplots(2, 2, figsize = (15,10))
        # for name in compare_name_list:
        #     y = self.calib_data[name][1]['EB1'][1]
        #     axarr1[0, 0].plot(y, label = name)
        # axarr1[0, 0].set_title('count EB1')
        # axarr1[0, 0].grid(True)
        # for name in compare_name_list:
        #     y = self.calib_data[name][1]['EB2'][1]
        #     axarr1[0, 1].plot(y, label = name)
        # axarr1[0, 1].set_title('count EB2')
        # axarr1[0, 1].grid(True)
        # for name in compare_name_list:
        #     y = self.calib_data[name][1]['EB3'][1]
        #     axarr1[1, 0].plot(y, label = name)
        # axarr1[1, 0].set_title('count EB3')
        # axarr1[1, 0].grid(True)
        # axarr1[1, 0].set_xlabel('Time')
        # for name in compare_name_list:
        #     y = self.calib_data[name][1]['EB4'][1]
        #     axarr1[1, 1].plot(y, label = name)
        # axarr1[1, 1].set_title('count EB4')
        # axarr1[1, 1].grid(True)
        # axarr1[1, 1].set_xlabel('Time')
        # # Fine-tune figure; hide x ticks for top plots and y ticks for right plots
        # plt.setp([a.get_xticklabels() for a in axarr1[0, :]], visible=False)
        # plt.setp([a.get_yticklabels() for a in axarr1[:, 1]], visible=False)
        #
        # plt.legend( loc='upper right' )
        # plt.savefig('count.png')




    # plot time space diagram for the result
    def plot_calib_data_time_space_diagram(self, compare_name_list, det_name_list, det_position, speed_axis_limit):

        # Plot the time-space diagram of the sensors
        # compare_name_list should be less than 4, put in one window
        # if len(compare_name_list) <= 4:
        #     f, axarr = plt.subplots(2, 2, figsize = (15,10))


        fig, axarr = plt.subplots(len(compare_name_list), 1, figsize=(10,10), dpi=100)

        subplot_index = 0
        im = []
        for name in compare_name_list:
            t_x_data = []
            for det in det_name_list:
                # append the speed data []
                t_x_data.append( self.calib_data[name][1][det][0] )

            t_x = np.array(t_x_data).astype(np.float)

            if name == 'valid':
                # axarr[subplot_index].set_title('Traffic condition from sensor measurement')
                axarr[subplot_index].set_ylabel('measurement', fontsize=18)
            elif name == 'optimal' or name == 'user':
                # axarr[subplot_index].set_title('Traffic condition from calibrated parameters')
                axarr[subplot_index].set_ylabel('calibrated', fontsize=18)
            elif name == 'default':
                # axarr[subplot_index].set_title('Traffic condition from default parameters')
                axarr[subplot_index].set_ylabel('default', fontsize=18)
            else:
                axarr[subplot_index].set_ylabel('{0}'.format(name), fontsize=18)


            if subplot_index == len(compare_name_list)-1:
                axarr[subplot_index].set_xlabel('time steps (5 min)', fontsize=18)

            im.append( axarr[subplot_index].imshow(t_x[::-1],cmap=plt.get_cmap('jet_r'),
                        interpolation='nearest',
                        vmin=0, vmax=80) )

            axarr[subplot_index].set_yticks(np.arange(0, len(det_name_list), 1.0))
            axarr[subplot_index].set_yticklabels(det_name_list[::-1])

            subplot_index += 1


        cax = fig.add_axes([0.9, 0.1, 0.02, 0.8])
        fig.colorbar(im[len(compare_name_list)-1], cax=cax)



        # fig = plt.figure(figsize=(16,8), dpi=100)
        # # plot the first one
        # name  = compare_name_list[1]
        # t_x_data = []
        # for det in det_name_list:
        #     # append the speed data []
        #     t_x_data.append( self.calib_data[name][1][det][0] )
        #
        #
        # t_x = np.array(t_x_data).astype(np.float)
        #
        # im = plt.imshow(t_x.T,cmap=plt.get_cmap('jet'),
        #                 interpolation='nearest',
        #                 vmin=0, vmax=80)
        #
        # # cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
        # # fig.colorbar(im[2], cax=cax)
        # # fig.colorbar()



        # fig = plt.figure(figsize=(16,8), dpi=100)
        # # plot the first one
        # name  = compare_name_list[2]
        # t_x_data = []
        # for det in det_name_list:
        #     # append the speed data []
        #     t_x_data.append( self.calib_data[name][1][det][0] )
        #
        #
        # t_x = np.array(t_x_data).astype(np.float)
        #
        # im = plt.imshow(t_x.T,cmap=plt.get_cmap('jet'),
        #                 interpolation='nearest',
        #                 vmin=0, vmax=80)

        plt.show()

        # else:
        #     print 'plot_calib_data_time_space_diagram does not support plotting more than four solutions'




    # plot the fundamental diagram SB1
    def plot_I57_SB1_FD(self, sol_name):

        SB1_speed =  np.array(self.calib_data[sol_name][1]['SB1'][0]).astype(float)
        SB1_flow =  np.array(self.calib_data[sol_name][1]['SB1'][1]).astype(float)*12.0

        SB1_density = SB1_flow/SB1_speed

        # speed
        fig_window = plt.figure(figsize=(15,10))
        fig = fig_window.add_subplot(111)


        x = SB1_density
        y = SB1_flow

        plt.plot(x, y, '*b')
        plt.title('Fundamental diagram for I57', fontsize=24)
        plt.xlabel('Traffic density (veh/mile)', fontsize=24)
        plt.ylabel('Traffic flow (veh/hr)', fontsize=24)
        # plt.xlim([ np.min(x) - 0.1*(np.max(x)-np.min(x)) , np.max(x) + 0.1*(np.max(x)-np.min(x))])
        plt.xlim([0, 14])
        plt.ylim([0, 650])
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



def plot_simulation_2d_mesh(result_file, solution_file, mat_size):
    """
    This function visualizes the 2d response surface to two parameters which is saved in result_file.
    It can also overlay the optimization steps as a directed arrow which is saved in solution_files
    :param result_file: a N by 3 matrix with header. First two columns are the parameter values, the last column is
        the objective funciton value.
    :param solution_file: the standard solution files. Each row comes with the parameters and then the objective value
    :param mat_size: a tuple, the levels of values of the two parameters in result_file
    :return:
    """

    # =======================================================
    # read data file
    # =======================================================
    x_grid = []
    y_grid = []
    f_val = []
    with open(result_file, 'r') as f:

        header = f.readline()

        header_items = header.strip().split(',')

        for line in f:

            items = line.strip().split(',')

            x_grid.append(float(items[0]))
            y_grid.append(float(items[1]))
            f_val.append(float(items[2]))

    x_mat = np.array(x_grid).reshape(mat_size)
    y_mat = np.array(y_grid).reshape(mat_size)
    f_mat = np.array(f_val).reshape(mat_size)

    # =======================================================
    # plot the solution steps
    # =======================================================
    sol_paras, sol_val = load_previous_opt_solutions(solution_file)

    solutions = OrderedDict()
    solutions[header_items[0]] = []
    solutions[header_items[1]] = []
    solutions['obj'] = []

    for sol_id in range(0,len(sol_paras)):

        para = sol_paras[sol_id]
        val = sol_val[sol_id]

        for key in para.keys():

            if key == header_items[0]:
                solutions[header_items[0]].append( para[key][0] )
            elif key == header_items[1]:
                solutions[header_items[1]].append( para[key][0] )

        solutions['obj'].append(val)

    # =======================================================
    # plot mesh surface
    # =======================================================
    fig = plt.figure(figsize=(13,10))
    ax = fig.gca(projection='3d')

    surf = ax.plot_surface(x_mat, y_mat, f_mat,
                           rstride=1, cstride=1, cmap=cm.jet,
                           linewidth=0, antialiased=False)

    # plot contour
    cset = ax.contour(x_mat, y_mat, f_mat, zdir='z', offset =-10, cmap=cm.jet)
    cset = ax.contour(x_mat, y_mat, f_mat, zdir='x', offset =0.6, cmap=cm.jet)
    cset = ax.contour(x_mat, y_mat, f_mat, zdir='y', offset =0.6, cmap=cm.jet)

    ax.set_xlim( 0.6, 1.3 )
    ax.set_ylim( 0.6, 1.1 )
    ax.set_zlim( -10 , np.max(f_val))

    # set axis
    ax.set_xlabel('{0}'.format(header_items[0]))
    ax.set_ylabel('{0}'.format(header_items[1]))

    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    fig.colorbar(surf, shrink=0.5, aspect=15)

    plt.draw()

    # =======================================================
    # plot contour
    # =======================================================
    levels = [0.01, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 90]
    fig = plt.figure(figsize=(13,10))
    cs = plt.contourf(x_mat, y_mat, f_mat, levels,
                 cmap=cm.jet)

    x = np.array(solutions[header_items[0]])
    y = np.array(solutions[header_items[1]])
    plt.quiver(x[:-1], y[:-1], x[1:]-x[:-1], y[1:]-y[:-1], scale_units='xy', angles='xy', scale=1)

    plt.xlabel('{0}'.format(header_items[0]))
    plt.ylabel('{0}'.format(header_items[1]))
    cbar = plt.colorbar(cs)

    plt.show()


def load_previous_opt_solutions(opt_step_file):
    """
    To save the optimization cost, if a value point has been evaluated in the simulation, it will be
    logged in txt file.
    A new optimization will start by loading previously explored value points. If a point has been
    previously evaluated, then it will NOT be re-simulated to save time.
    :param opt_step_file:
    :return: [paras_list, obj_list]
    """
    if opt_step_file is None:
        return None

    f = open(opt_step_file, 'r')

    paras_list = []
    obj_value_list = []

    for line in f:

        paras = OrderedDict()

        line = line.strip()
        items = line.split(';')

        # parse each line. first item is the counter, skip
        for i in range(1, len(items)):
            entry = items[i]

            key = entry.split(':')[0]

            if key == 'RMSE_avg':
                obj_val = float(entry.split(':')[1])
            elif key == 'RMSEs':
                continue
            else:
                # the parameters
                paras[key] = [float(v) for v in entry.split(':')[1].split(',')]

        paras_list.append(paras)
        obj_value_list.append(obj_val)

    f.close()

    return [paras_list, obj_value_list]



def plot_high_dimension(opt_step_file, list_paras, true_para_sol, para_bounds):
    """
    This function plots the high dimensional solution in a parallel coordinates.
    :param opt_step_file: the file name
    :param list_paras: list of strings, each string is the key in paras to plot
    :param true_para_sol: tuple, (true_para, true_sol)
    :param para_bounds: dict; key is para.keys(), value is a [min, max]
    :return: a plot
    """

    # ==========================================================================
    # Read previous solutions
    # ==========================================================================
    step_paras, step_sol = load_previous_opt_solutions(opt_step_file)

    # ==========================================================================
    # normalize and prepare the para data using the bounds
    # list of normalized paras, each entry is a set of normalized parameters
    norm_para = []
    for each_step_para in step_paras:

        para_val = []
        for para_name in list_paras:

            norm_v = (each_step_para[para_name][0]-para_bounds[para_name][0])\
                     /(para_bounds[para_name][1] - para_bounds[para_name][0])
            para_val.append( norm_v )


        norm_para.append(para_val)

    # normalize the true value
    if true_para_sol is not None:
        true_para, true_sol = true_para_sol
        norm_true = []
        for para_name in list_paras:
            norm_v = (true_para[para_name][0] - para_bounds[para_name][0])\
                     /(para_bounds[para_name][1] - para_bounds[para_name][0])
            norm_true.append( norm_v )

    # print(norm_true)

    # ==========================================================================
    # Plot the parallel coordinates
    # ==========================================================================
    dim = len(list_paras)
    x = np.arange(0, dim)

    fig = plt.figure(figsize=(30,8))
    ax = fig.add_subplot(111)

    # normalize colormap
    step_sol = np.array(step_sol)
    cNorm = colors.Normalize(vmin=np.min(step_sol[:,0]) , vmax=np.max(step_sol[:,0]))
    scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cm.jet)
    scalarMap.set_array([])

    # Z = [[0,0],[0,0]]
    # levels = range( 0 , 125, 5 )
    # CS3 = plt.contourf( Z, levels, cmap=cm.jet )
    # plt.clf()

    for idx in range(0,len(norm_para)):

        # step_sol[idx] = [RMS_speed, RMS_count]
        color_val = scalarMap.to_rgba( step_sol[idx][0] )
        ax.plot( x, norm_para[idx], color=color_val )

    if true_para_sol is not None:
        plt.title('True paras: {0}'.format( ','.join('{0:.2f}'.format(i) for i in norm_true )) )
    else:
        plt.title('Optimization steps')

    ax.set_xticklabels(list_paras)
    ax.grid()
    plt.colorbar(scalarMap)
    plt.show()






def find_feasible_optimum(opt_step_file, para_bounds):
    """
    This function finds the optimal parameter within the feasible bounds.
    :param opt_step_file: the file name
    :param para_bounds: dict; key is para.keys(), value is a [min, max]
    :return: a plot
    """

    # ==========================================================================
    # Read previous solutions
    # ==========================================================================
    step_paras, step_sol = load_previous_opt_solutions(opt_step_file)

    # optimal so far
    optimal_para = None
    optimal_val = np.inf

    unique_paras = []
    unique_vals = []
    num_feasible = 0
    # remove duplicates, and find the minimum
    for para_set, val in zip(step_paras, step_sol):

        if para_set in unique_paras:
            continue
        else:
            unique_paras.append(para_set)
            unique_vals.append(val)

            # check if the value is feasible in the bounds
            feasible_flag = True
            for para_name in para_set:
                # print('paraset[key]:{0}'.format(para_set[para_name]))
                if para_set[para_name][0] < para_bounds[para_name][0] or \
                    para_set[para_name][0] > para_bounds[para_name][1]:
                    feasible_flag = False

            if feasible_flag:
                num_feasible += 1
            else:
                continue

            # only consider the speed RMSE
            if val < optimal_val:
                optimal_para = para_set
                optimal_val = val
            elif val == optimal_val:
                print('Warning: two sets of parameters achieving same objective values')
                print('Another set of paras with obj:{0}'.format(val))
                for para_name in para_set:
                    print('--- {0}:{1}'.format(para_name, para_set[para_name]))

    print('Total number of unique steps {0}/{1}'.format(len(unique_paras), len(step_paras)))
    print('Total number of feasible steps {0}/{1}'.format(num_feasible, len(step_paras)))
    print('Optimal parameters with obj: {0}'.format(optimal_val))
    for para_name in optimal_para:
        print('--- {0}:{1}'.format(para_name, optimal_para[para_name]))





def __isNumber(s):
    """
    This function tests if string s is a number
    :param s: string
    :return:
    """
    try:
        float(s)
        return True
    except ValueError:
        return False









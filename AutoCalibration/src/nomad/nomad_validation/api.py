__author__ = 'Yanning Li'
"""
This script is used to validate NOMAD by several nonlinear optimization functions.
"""

# It will later be used for communicating with AIMSUN by writing parameters to files and reading objective value from AIMSUN file.

# run aimsun.exe -script thispythonfile.py and at the same time, run the other side of the code

import numpy as np
import os
from os.path import exists
import time
import sys
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



# this function reads the parameters
def read_paras(parafile='C:\paras.txt'):

    while not exists(parafile):
        print 'waiting for paras...'
        time.sleep(0.1)   # sleep 1 second

    if exists('C:\paras.txt'):
        f = open('C:\paras.txt', 'r')
        line = f.readline()

        tup = line.split(',')

        x1 = float(tup[0])
        x2 = float(tup[1])

        f.close()

    # delete the file once read
    os.remove(parafile)

    return (x1, x2)

# this function writes the objective function value to file simval.txt
def write_simval(simval=0, simvalfile='C:\simval.txt'):

    f = open(simvalfile,'w')
    f.write(str.format("{0:.16f}", simval))
    f.close()

    return 0



def compute_objective_peaks(paras):
    """
    Peaks function is a typical MATLAB 2d function. 
    optimal solution is (0.2283, -1.6256) = -6.5511
    :param paras: (x1,x2) -3 <= x1, x2 <= 3, 
    :return: val
    """
    x1 = paras[0]
    x2 = paras[1]

    val = 3*pow(1-x1, 2)*np.exp(-pow(x1, 2) - pow(x2+1, 2)) - \
          10*(x1/5 - pow(x1, 3) - pow(x2, 5))*np.exp(-pow(x1,2)-pow(x2, 2)) - \
          1/3*np.exp(-pow(x1+1, 2) - pow(x2, 2))

    return val


def compute_objective_eggholder(paras):
    """
    Eggholder is a 2D function which looks a lot like the parameter calibration problem for AIMSUN.
    Search Eggholder function in Wiki.
    Optimal solution is f(512, 404.2319) = -959.6407
    :param paras: (x,y), -512 <= x, y <= 512
    :return:
    """
    x, y = paras

    val = -(y+47)*np.sin( np.sqrt( np.abs( x/2 + (y+47) ) ) ) \
        - x*np.sin( np.sqrt( np.abs( x-(y+47) ) ) )

    return val



def plot_egghold_function():
    """
    This function plots the egg hold function in 2d mesh
    :return:
    """
    x = np.linspace(-512, 512, 200)
    y = np.linspace(-512, 512, 200)

    X, Y = np.meshgrid(x, y)
    Z = -(Y+47)*np.sin( np.sqrt( np.abs( X/2 + (Y+47) ) ) ) \
        - X*np.sin( np.sqrt( np.abs( X-(Y+47) ) ) )

    print 'computed Z'

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.jet,
                       linewidth=0, antialiased=False)
    ax.set_zlim(-1000, 1000)

    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.show()



# Here is the main function
def main(argv):

    # plot_egghold_function()
    # paras = read_paras(parafile='C:\paras.txt')

    # val = compute_objective(paras)

    # print 'paras: {0} and obj value:{1}\n'.format(paras, val)

    # write_simval(val, simvalfile='C:\simval.txt')

    print 0

if __name__ == '__main__':
    main(sys.argv[1:])















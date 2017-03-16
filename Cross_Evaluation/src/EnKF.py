import numpy as np
from copy import deepcopy

__author__ = 'Yanning Li'
"""
This is the abstract Ensemble Kalman filter class. It contains the necessary components for the iteration of
the EnKF. For each specific problem, the user need to create a subclass and overwrite functions for the state evolution equations
and the observation equation.
    - The equations are based on paper "An iterative ensemble Kalman Filter" by Rolf J. Lorentzen and Geir Nevdal, 2011
    - This approach uses an augmented state to account for the nonlinear observation equations.
    - To reduce the computational cost, this implementation uses the efficient version of the augmented state EnKF scheme.
"""


class EnKF:
    """
    This is the abstract class for EnKF. The user need to create a subclass for each specific problem with the following
    functions over written:
    - f_model()
    - h_obs()

    This class only gives the single step iteration, hence can be used to process streaming data

    """

    def __init__(self, dim_state=0, num_ensembles=0):
        """
        The constructor of a nonlinear EnKF. As a start, it should know the dimension of the states and observations,
            and the number of ensembles that will be used.
        :param dim_state: int, the dimension of the system states
        :param num_ensembles: int, the number of ensembles that will be used
        :return:
        """

        self.dim_state = dim_state
        self.num_ensembles = num_ensembles

        # The system state variables. See the PDF for the meaning of the variables
        # the state in forecast step. N is the dimension of the state and q is the number of ensembles
        # N x q matrix, state denoted by x, forecast denoted by _f
        self.X_f = np.matrix(np.zeros((self.dim_state,self.num_ensembles),float))

        # N x q matrix the state in analysis step. N is the dimension of the state and q is the number of ensembles
        # state analysis, denoted by x_a
        self.X_a = np.matrix(np.zeros((self.dim_state,self.num_ensembles),float))

        # the estimated state
        self.state_est = np.matrix( np.zeros( (self.dim_state, 1), float ) )

    def set_initial_ensembles(self, X_a):
        """
        This function generates the initial ensembles which are the disturbed initial states.
        :param X_a: np.matrix, dim_state x num_ensembles, initial ensembles.
        :return:
        """

        if X_a is not None:
            # if not None, then update. Otherwise simply set as 0
            tmp_X_a = deepcopy(X_a)
            if type(tmp_X_a) is not np.matrix:
                tmp_X_a = np.matrix(tmp_X_a)

            if tmp_X_a.shape != ( self.dim_state, self.num_ensembles ):
                raise Exception('Error: The dimension of the initial ensembles must be dim_state x num_ensembles.')

            self.X_a = tmp_X_a

    # @profile
    def update_estimate(self, y_obs, cov_noise, cur_time):
        """
        This is the function that should be called to update a new estimate x_a^k, given y_k,o
        :param y_obs: np.matrix, num_obs x 1, the observed measurement
        :param cov_noise: np.matrix, the covariance of noise of measurement.
        :param cur_time: seconds, current time of the observation
        :return:
        """

        if y_obs is None or y_obs.size == 0:
            # no observation, simply set x_a = x_f

            for j in range(0, self.num_ensembles):
                # propagate the state for each ensemble based on the model
                self.X_f[:,j] = self.f_model( self.X_a[:,j], j )

            # no analysis step, simply copy X_f to X_a (NOTE: must use deep copy)
            self.X_a = deepcopy(self.X_f)

            return self.X_a.mean(1)

        else:
            # if there is observation
            # get the dimension of the observation
            dim_obs = y_obs.size

            Y_f = np.matrix( np.zeros((dim_obs, self.num_ensembles), float) )
            Y_obs_ensembles = np.matrix( np.zeros((dim_obs, self.num_ensembles), float) )

            for j in range(0, self.num_ensembles):

                # print 'Status: Updating ensemble {0}'.format(j)
                # propagate the state for each ensemble
                self.X_f[:, j] = self.f_model( self.X_a[:, j], j )
                # propagate the forecast observation
                Y_f[:, j] = self.h_obs( cur_time, j )
                # Y_f[:, j] = self.h_obs( self.X_f[:,j] )   # The standard observation equation

                # generate the observation ensemble
                Y_obs_ensembles[:, j] = np.matrix(y_obs).reshape(dim_obs,1) + np.matrix(np.random.multivariate_normal( np.zeros(dim_obs),
                                                                                 cov_noise)).reshape(dim_obs,1)

            # TODO: a quick hack: all observation should be positive
            # Y_obs_ensembles[Y_obs_ensembles < 0] = 0

            K = self.__compute_Kalman_gain(self.X_f, Y_f, cov_noise)

            # update the analysis state
            self.X_a = self.X_f + K*(Y_obs_ensembles - Y_f)

            # return the estimate
            self.state_est = self.X_a.mean(1)
            return self.state_est

    @staticmethod
    def __compute_Kalman_gain(X_f, Y_f, R):
        """
        This function computes the Kalman gain
        :param X_f: np.matrix, dim_state x num_ensembles, the forecast state
        :param Y_f: np.matrix, dim_obs x num_ensembles, the forecast measurement
        :param R: np.matrix, dim_obs x dim_obs, the noise covariance
        :return: K, np.matrix, dim_state x dim_obs
        """

        X_f = np.matrix(X_f)
        Y_f = np.matrix(Y_f)
        R = np.matrix(R)

        # get number of ensembles
        J = X_f.shape[1]

        # compute the innovation sequence
        A = X_f - X_f.mean(1)
        B = Y_f - Y_f.mean(1)

        # compute the Kalman gain
        K = A*B.T*np.linalg.inv( B*B.T/(J-1.0) + R )\
                 /(J-1.0)

        return K

    def f_model(self, x_a, e_id):
        """
        This function propagates the state of all ensembles. This is the state propagation equation F(x)
        This function need to be overwritten for each specific problem.
        :param x_a: np.matrix, dim_state x 1, the analysis state in the previous step
        :param e_id: the ensemble id
        :return: x_f, np.matrix, dim_state x 1, the forecast state in the next step
        """
        raise Exception('Error: The system state evolution function f_model not defined.')

    def h_obs(self, cur_time, e_id):
        """
        This function computes the forecast observation from the forecast state H(x)
        This function need to be overwritten for each specific problem.
        :param cur_time: the current time step
        :param e_id: the ensemble id
        :return: y_f, np.matrix, dim_obs x 1, the forecast observation in the next step
        """
        raise Exception('Error: The system observation function h_obs not defined.')















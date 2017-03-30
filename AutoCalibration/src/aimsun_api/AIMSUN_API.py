import csv
import os
import random
import sys
import time
from collections import OrderedDict
from datetime import datetime
from os.path import exists

SITEPACKAGES = "C:\\Python27\\Lib\\site-packages"
sys.path.append(SITEPACKAGES)
import numpy as np

from PyANGBasic import *
from PyANGKernel import *
from PyANGConsole import *
from PyANGAimsun import *

__author__ = 'Yanning Li, Juan Carlos Martinez Mori'
"""
This script contains functions requried to automate AIMSUN.
"""


# ====================================================================================================
# Uncomment the following are just to remove the error messages

# def QString():
#     pass
#
#
# def QVariant():
#     pass
#
#
# def GKScheduleDemandItem():
#     pass
#
#
# GKTimeDuration = None
# GKVehicle = None
# GKSystem = None
# GKSection = None
# GKExperiment = None
# GKVehicleReactionTimes = None
#
#
# def GAimsunSimulator():
#     pass
#
#
# GKSimulationTask = None
# GKReplication = None
# GGetramModule = None
# GKTimeSerieIndex = None
# GK = None
# GKColumnIds = None
# Qt = None
# QTime = None
# GKGenericExperiment = None

# ====================================================================================================
# Define the columns for the files to be read
# pars.txt
INDEX_PARS_STATE = 0
INDEX_PARS_VEH = 1
INDEX_PARS_FROMTIME = 2
INDEX_PARS_DUR = 3
# flows.txt
INDEX_FLOWS_STATE = 0
INDEX_FLOWS_SECID = 1
INDEX_FLOWS_INFLOW = 2
# turns.txt
INDEX_TURNS_STATE = 0
INDEX_TURNS_FROM = 1
INDEX_TURNS_TO = 2
INDEX_TURNS_PERC = 3

KMH2MPH = 0.62137

_debug = False
_show_detector_data = False


class AimsunApi:
    """
    This class contains the methods that directly interact with AIMSUN
    """

    def __init__(self, cmd_logger):
        """
        Initialize the api with the command window logger
        :param cmd_logger: logger for porting cmd output to text file
        :return: an object
        """
        self.cmd_logger = cmd_logger
        self.model = None

    # ====================================================================================================
    # --------------------------- Functions for reading demand or data from files ------------------------
    def read_demand_from_file(self, demand_name, pars_file, flows_file, turns_file,
                              main_entrance_id, main_entrance_discount_ratio):
        """
        Read the demand from files to create traffic states.
        :param demand_name: create the demand for this specific one
        :param pars_file: the file for pars
        :param flows_file: the file for flows
        :param turns_file: the file for turns
        :param main_entrance_id: the main entrance id
        :param main_entrance_discount_ratio: the discount ratio, e.g., 1.0 no change. 1.1, 10% higher
        :return: aimsun demand object
        """

        # Set up a Traffic Demand and add it to the model
        if demand_name is not None:
            demand = self.model.getCatalog().findByName(QString(demand_name))
        else:
            demand = None

        if demand is None or demand.isA(QString("GKTrafficDemand")) is False:
            # create new demand
            demand = GKSystem.getSystem().newObject("GKTrafficDemand", self.model, -1, True)
            demand.setName(QString(demand_name))
            self.cmd_logger.print_cmd('Demand {0} not found. Creating new one...'.format(demand_name))

            # create under demand folder
            folder = self.__getDemandFolder()
            folder.append(demand)

        # clear schedule in demand to prevent overlapping
        demand.removeSchedule()

        # Read pars.txt into a dictionary
        pars_dict = self.__readParsFile(pars_file)
        flows_dict = self.__readFlowsFile(flows_file)
        turns_dict = self.__readTurnsFile(turns_file)

        for state_name in pars_dict.keys():
            state = self.__createState(state_name, pars_dict, flows_dict, turns_dict, main_entrance_id,
                                       main_entrance_discount_ratio)

            if _debug:
                self.cmd_logger.print_cmd(
                    'state from {0}, duration:{1}\n'.format(state.getFrom().toString(), state.getDuration().toString()))

                self.cmd_logger.print_cmd(
                    'state entrance flow {0}'.format(state.getEntranceFlow(self.model.getCatalog().find(int(330)), None)))

                self.cmd_logger.print_cmd('state turn 330->(340, 341): ({0},{1})'.format(
                    state.getTurningPercentage(self.model.getCatalog().find(int(330)),
                                               self.model.getCatalog().find(int(340)), None),
                    state.getTurningPercentage(self.model.getCatalog().find(int(330)),
                                               self.model.getCatalog().find(int(341)), None)))

            schedule = self.__createScheduleItem(state)

            demand.addToSchedule(schedule)

        # append the demand to the demands folder
        # folder = __getDemandFolder(model)
        # folder.append(demand)

        return demand

    def load_demand_from_ang(self, demand_name):
        """
        Load demand from the ang file with demand name
        :param demand_name: str, demand name
        :return:
        """
        # find the demand from the model
        demand = self.model.getCatalog().findByName(QString(demand_name))
        if demand is None or not demand.isA(QString("GKTrafficDemand")):
            self.cmd_logger.print_cmd('Error: no traffic demand named {0}\n'.format(demand_name))
            return None

        self.cmd_logger.print_cmd('Loaded demand named {0}\n'.format(demand_name))

        return demand

    def read_validation_data(self, start_time_str, end_time_str, name_strs, valid_file_path):
        """
        This code is to read validation data from csv file
        data format should be as follows: Make SURE there is Header
        Date/Time (UTC-06:00) Central Time (US & Canada),RadarSpeed (MPH),RadarVehiclesCount in the interval
        4/26/2015 15:55,,
        4/26/2015 16:00,72.23863137,51
        4/26/2015 16:05,71.23257753,91
        Note: to read data for a simulation from 03:00 to 04:00 (AIMSUN aggragate data at 03:05 + 00:05...)
        need to read 12 lines with start_time and end_time 03:00~to 04:00 (Carlos shifted the time 5 min earlier)

        :param start_time_str: "4/26/2015 16:00", string
        :param end_time_str:"4/26/2015 16:00", string; if start and end both None, the read all rows
        :param name_strs: a list of detector names
        :param valid_file_path: '', then current folder same as the script
        :return: dict: validation_data[detector_name] = np.array([[speed (mph)],[count (in interval)]])
        """
        dict_valid = OrderedDict()

        # if *_time_str is None, then assume the entire file is the validation data
        if start_time_str is None or end_time_str is None:

            for name in name_strs:

                # first list speed, second list count
                dict_valid[name] = [[], []]

                file_name = valid_file_path + name + '.csv'
                f_handle = open(file_name, 'r')

                data_set = csv.reader(f_handle, delimiter=',')
                # skip header
                next(data_set, None)

                for row in data_set:
                    dict_valid[name][0].append(float(row[1]))
                    dict_valid[name][1].append(float(row[2]))

                f_handle.close()

        # Otherwise on read data in a specific time period
        else:
            dt_format = '%m/%d/%Y %H:%M'
            t_start = datetime.strptime(start_time_str, dt_format)
            t_end = datetime.strptime(end_time_str, dt_format)

            for name in name_strs:
                # first list speed, second list count
                dict_valid[name] = [[], []]

                file_name = valid_file_path + name + '.csv'
                f_handle = open(file_name, 'r')

                data_set = csv.reader(f_handle, delimiter=',')
                # skip header
                next(data_set, None)

                for row in data_set:

                    cur_dt = datetime.strptime(row[0], dt_format)

                    if t_start <= cur_dt < t_end:
                        # no need to save the time
                        dict_valid[name][0].append(float(row[1]))
                        dict_valid[name][1].append(float(row[2]))

                f_handle.close()

        return dict_valid

    # ---------------------- Functions for setting up the simulation in AIMSUN ------------------------
    def set_model(self, model):
        """
        This function sets the aimsun GK Model
        :param model: GK Model
        :return: set in property
        """
        self.model = model

    def setup_scenario(self, scenario_name, demand):
        """
        Set up the scenario and connect scenario with demand
        :param scenario_name: string, the scenario name, must exist in ang file, or create a new one
        :param demand: demand object, The demand for this scenario
        :return:
        """
        self.cmd_logger.print_cmd('\nSetting up scenario...')

        scenario = self.model.getCatalog().findByName(QString(scenario_name))
        if scenario is None or not scenario.isA(QString("GKScenario")):
            scenario = GKSystem.getSystem().newObject("GKScenario", self.model, -1, True)
            scenario.setName(QString(scenario_name))
            self.cmd_logger.print_cmd(
                'Error: no traffic scenario named {0}. Creating new one...\n'.format(scenario_name))

        scenario.setDemand(demand)

        # append the state to the state folder
        # folder = __getScenariosFolder(model)
        # folder.append(scenario)

        # set parameters here
        # parameters are set in the ScenarioInput data class
        paras = scenario.getInputData()

        # set the detection and statistical intervals as 5 min
        det_interval = "00:05:00"
        paras.setDetectionInterval(GKTimeDuration.fromString(QString(det_interval)))
        paras.setStatisticalInterval(GKTimeDuration.fromString(QString(det_interval)))

        self.cmd_logger.print_cmd('---- Detection interval is set as : {0}\n'.format(det_interval))

        return scenario

    def setup_experiment(self, experiment_name, scenario):
        """
        This function sets up the experiment under scenario
        :param experiment_name: Name of the experiment
        :param scenario: the scenario name
        :return: the GK experiment object
        """
        self.cmd_logger.print_cmd('\nSetting up experiment...\n')

        experiment = self.model.getCatalog().findByName(QString(experiment_name))
        if experiment is None or not experiment.isA(QString("GKExperiment")):
            experiment = GKSystem.getSystem().newObject("GKExperiment", self.model, -1, True)
            self.cmd_logger.print_cmd(
                'ERROR: No traffic experiment named {0}. Creating a new one...\n'.format(experiment_name))

            # attach the new experiment to folder
            folder = self.__getScenariosFolder()
            folder.append(experiment)

        experiment.setScenario(scenario)

        return experiment

    def setup_replication(self, experiment, num_rep, seed_list):
        """
        This function sets up the replications
        :param experiment: experiment
        :param num_rep: number of replications
        :param seed_list: the seed for each replication
        :return: the ave_result replication
        """
        self.cmd_logger.print_cmd('\nSetting up replications...')

        if experiment != None and experiment.isA("GKExperiment") \
                and experiment.getSimulatorEngine() == GKExperiment.eMicro:

            # add replications here
            replication_list = experiment.getReplications()

            # ===============================
            # create new replications
            if len(replication_list) == 0:
                # create replications
                for i in range(0, num_rep):
                    replication = GKSystem.getSystem().newObject("GKReplication", self.model, -1, True)
                    replication.setExperiment(experiment)
                    replication_list.append(replication)

                    if seed_list is not None:
                        replication.setRandomSeed(seed_list[i])
                    self.cmd_logger.print_cmd('---- Created replication {0} with seed {1}'.format(replication.getId(),
                                                                                                  replication.getRandomSeed()))
            else:
                # show replcations:
                self.cmd_logger.print_cmd('---- Reloading {0} replications: {1} \n'.format(len(replication_list),
                                                                                           [replication.getId() for
                                                                                            replication in
                                                                                            replication_list]))

            # create the average experiment result
            avg_result = GKSystem.getSystem().newObject("GKExperimentResult", self.model)
            avg_result.setName('average_result')
            self.cmd_logger.print_cmd('Created new average replication: {0}'.format(avg_result.getName()))
            # print_cmd('Total number of replications is: {0}',format(len(experiment.getReplications()))

            # set the experiment of this result object
            avg_result.setExperiment(experiment)
            # add replcations to the average
            for replication in replication_list:
                avg_result.addReplication(replication)
                self.cmd_logger.print_cmd(
                    '---- Added replication {0} to {1}'.format(replication.getId(), avg_result.getName()))

            # compute the average; add to the experiment.
            experiment.addReplication(avg_result)

            return avg_result

    def set_random_seed(self, avg_result):
        """
        This function reset the seed for avg_result replications
        :param avg_result: the avg_result replication
        :param cmd_logger: logger for cmd output
        :return:
        """
        replication_list = avg_result.getReplications()

        self.cmd_logger.print_cmd('\nResetting seeds for replications:')

        i = 0
        for replication in replication_list:
            replication.setRandomSeed(random.randint(0, 10000))
            self.cmd_logger.print_cmd('----Reset replication {0} with seed {1}'.format(replication.getId(),
                                                                                       replication.getRandomSeed()))
            i += 1

    def create_simulator(self):
        """
        Create a simulator for model
        :return: simulator for the model
        """

        simulator = GAimsunSimulator()
        simulator.setModel(self.model)

        return simulator

    def set_new_paras(self, experiment, paras):
        """
        This function sets up the parameters for experiment
        :param experiment: the experiment
        :param paras: a dict, parameters
        :return:
        """
        # Unique car id is 53
        car_type = self.model.getCatalog().find(53)
        if not car_type.isA(QString("GKVehicle")):
            self.cmd_logger.print_cmd('Error: Car type is not correct: demand can not be correctly loaded\n')
            self.cmd_logger.print_cmd(type(car_type))
            # debug_catalog(model)
            return None

        # unique truck id is 56
        truck_type = self.model.getCatalog().find(56)
        if not truck_type.isA(QString("GKVehicle")):
            self.cmd_logger.print_cmd('Error: Truck type is not correct: demand can not be correctly loaded\n')
            return None

        # note the parameter should be in the same format as the preset parameters
        # even if we are only calibrating the mean
        for key in paras.keys():

            para_name = key.split('_')
            para_value = paras[key]

            if para_name[0] == 'car':
                # car paras
                # calibrate speedAcceptanceMean and speedAcceptanceDev for freeflow
                if para_name[1] == 'speedAcceptance':
                    car_type.setDataValueByID(GKVehicle.speedAcceptanceMean, QVariant(para_value[0]))
                    car_type.setDataValueByID(GKVehicle.speedAcceptanceDev, QVariant(para_value[1]))
                    car_type.setDataValueByID(GKVehicle.speedAcceptanceMin, QVariant(para_value[2]))
                    car_type.setDataValueByID(GKVehicle.speedAcceptanceMax, QVariant(para_value[3]))
                    self.cmd_logger.print_cmd('---- set {0}: {1}'.format(key, para_value))

                elif para_name[1] == 'minHeadway':
                    car_type.setDataValueByID(GKVehicle.minimunHeadwayMean, QVariant(para_value[0]))
                    car_type.setDataValueByID(GKVehicle.minimunHeadwayDev, QVariant(para_value[1]))
                    car_type.setDataValueByID(GKVehicle.minimunHeadwayMin, QVariant(para_value[2]))
                    car_type.setDataValueByID(GKVehicle.minimunHeadwayMax, QVariant(para_value[3]))
                    self.cmd_logger.print_cmd('---- set {0}: {1}'.format(key, para_value))

                # calibrate maxAccelMean for congflow
                elif para_name[1] == 'maxAccel':
                    car_type.setDataValueByID(GKVehicle.maxAccelMean, QVariant(para_value[0]))
                    car_type.setDataValueByID(GKVehicle.maxAccelDev, QVariant(para_value[1]))
                    car_type.setDataValueByID(GKVehicle.maxAccelMin, QVariant(para_value[2]))
                    car_type.setDataValueByID(GKVehicle.maxAccelMax, QVariant(para_value[3]))
                    self.cmd_logger.print_cmd('---- set {0}: {1}'.format(key, para_value))

                # calibrate reaction_time for congflow
                elif para_name[1] == 'reactionTime':
                    # [reaction_time, reaction_stop, reaction_light, reaction_prob]
                    car_react = GKVehicleReactionTimes(para_value[0], para_value[1],
                                                       para_value[2], para_value[3])

                    car_type.setVariableReactionTimes([car_react])
                    experiment.setVariableReactionTimesMicro(car_type, [car_react])
                    self.cmd_logger.print_cmd('---- set {0}: {1}'.format(key, para_value))

                # calibrate minDistMean for congflow
                elif para_name[1] == 'minDist':
                    car_type.setDataValueByID(GKVehicle.minDistMean, QVariant(para_value[0]))
                    car_type.setDataValueByID(GKVehicle.minDistDev, QVariant(para_value[1]))
                    car_type.setDataValueByID(GKVehicle.minDistMin, QVariant(para_value[2]))
                    car_type.setDataValueByID(GKVehicle.minDistMax, QVariant(para_value[3]))
                    self.cmd_logger.print_cmd('---- set {0}: {1}'.format(key, para_value))

                # calibrate sensitivityFactorMean = min = max for congflow
                elif para_name[1] == 'sensitivityFactor':
                    car_type.setDataValueByID(GKVehicle.sensitivityFactorMean, QVariant(para_value[0]))
                    car_type.setDataValueByID(GKVehicle.sensitivityFactorDev, QVariant(para_value[1]))
                    car_type.setDataValueByID(GKVehicle.sensitivityFactorMin, QVariant(para_value[2]))
                    car_type.setDataValueByID(GKVehicle.sensitivityFactorMax, QVariant(para_value[3]))
                    self.cmd_logger.print_cmd('---- set {0}: {1}'.format(key, para_value))

                # This should already been preset as unrealistically high such that speed can be fully
                # controlled by the speed acceptance
                elif para_name[1] == 'maxSpeed':
                    car_type.setDataValueByID(GKVehicle.maxSpeedMean, QVariant(para_value[0]))
                    car_type.setDataValueByID(GKVehicle.maxSpeedDev, QVariant(para_value[1]))
                    car_type.setDataValueByID(GKVehicle.maxSpeedMin, QVariant(para_value[2]))
                    car_type.setDataValueByID(GKVehicle.maxSpeedMax, QVariant(para_value[3]))
                    self.cmd_logger.print_cmd('---- set {0}: {1}'.format(key, para_value))

                else:
                    self.cmd_logger.print_cmd(
                        '\n---- ERROR: Could not recognize preset parameter entry {0}: {1}\n'.format(key, para_value))

            elif para_name[0] == 'truck':
                # truck paras
                if para_name[1] == 'speedAcceptance':
                    truck_type.setDataValueByID(GKVehicle.speedAcceptanceMean, QVariant(para_value[0]))
                    truck_type.setDataValueByID(GKVehicle.speedAcceptanceDev, QVariant(para_value[1]))
                    truck_type.setDataValueByID(GKVehicle.speedAcceptanceMin, QVariant(para_value[2]))
                    truck_type.setDataValueByID(GKVehicle.speedAcceptanceMax, QVariant(para_value[3]))
                    self.cmd_logger.print_cmd('---- set {0}: {1}'.format(key, para_value))

                elif para_name[1] == 'minHeadway':
                    truck_type.setDataValueByID(GKVehicle.minimunHeadwayMean, QVariant(para_value[0]))
                    truck_type.setDataValueByID(GKVehicle.minimunHeadwayDev, QVariant(para_value[1]))
                    truck_type.setDataValueByID(GKVehicle.minimunHeadwayMin, QVariant(para_value[2]))
                    truck_type.setDataValueByID(GKVehicle.minimunHeadwayMax, QVariant(para_value[3]))
                    self.cmd_logger.print_cmd('---- set {0}: {1}'.format(key, para_value))

                elif para_name[1] == 'maxAccel':
                    truck_type.setDataValueByID(GKVehicle.maxAccelMean, QVariant(para_value[0]))
                    truck_type.setDataValueByID(GKVehicle.maxAccelDev, QVariant(para_value[1]))
                    truck_type.setDataValueByID(GKVehicle.maxAccelMin, QVariant(para_value[2]))
                    truck_type.setDataValueByID(GKVehicle.maxAccelMax, QVariant(para_value[3]))
                    self.cmd_logger.print_cmd('---- set {0}: {1}'.format(key, para_value))

                elif para_name[1] == 'reactionTime':
                    # [reaction_time, reaction_stop, reaction_light, reaction_prob]
                    truck_react = GKVehicleReactionTimes(para_value[0], para_value[1],
                                                         para_value[2], para_value[3])

                    truck_type.setVariableReactionTimes([truck_react])
                    experiment.setVariableReactionTimesMicro(truck_type, [truck_react])
                    self.cmd_logger.print_cmd('---- set {0}: {1}'.format(key, para_value))

                elif para_name[1] == 'minDist':
                    truck_type.setDataValueByID(GKVehicle.minDistMean, QVariant(para_value[0]))
                    truck_type.setDataValueByID(GKVehicle.minDistDev, QVariant(para_value[1]))
                    truck_type.setDataValueByID(GKVehicle.minDistMin, QVariant(para_value[2]))
                    truck_type.setDataValueByID(GKVehicle.minDistMax, QVariant(para_value[3]))
                    self.cmd_logger.print_cmd('---- set {0}: {1}'.format(key, para_value))

                elif para_name[1] == 'sensitivityFactor':
                    truck_type.setDataValueByID(GKVehicle.sensitivityFactorMean, QVariant(para_value[0]))
                    truck_type.setDataValueByID(GKVehicle.sensitivityFactorDev, QVariant(para_value[1]))
                    truck_type.setDataValueByID(GKVehicle.sensitivityFactorMin, QVariant(para_value[2]))
                    truck_type.setDataValueByID(GKVehicle.sensitivityFactorMax, QVariant(para_value[3]))
                    self.cmd_logger.print_cmd('---- set {0}: {1}'.format(key, para_value))

                elif para_name[1] == 'maxSpeed':
                    truck_type.setDataValueByID(GKVehicle.maxSpeedMean, QVariant(para_value[0]))
                    truck_type.setDataValueByID(GKVehicle.maxSpeedDev, QVariant(para_value[1]))
                    truck_type.setDataValueByID(GKVehicle.maxSpeedMin, QVariant(para_value[2]))
                    truck_type.setDataValueByID(GKVehicle.maxSpeedMax, QVariant(para_value[3]))
                    self.cmd_logger.print_cmd('---- set {0}: {1}'.format(key, para_value))

                else:
                    self.cmd_logger.print_cmd(
                        '\n---- ERROR: Could not recognize preset parameter entry {0}: {1}\n'.format(key, para_value))

            else:
                raise Exception('ERROR: Could not recognize preset parameter entry {0}\n'.format(key))

    # ---------------------- Functions for running a simulation -------------------------------
    def simulate_rep_from_paras(self, experiment, paras,
                                simulator, avg_result, plugin,
                                valid_data, det_weight,
                                name):
        """
        This is the function that used for simulate and return the simulation result.
        It can be used in two ways:
        - Generate validation data if valid_data is None
        - Simulate and output the objective function value
        :param model: GK Model
        :param experiment: the GK Experiment
        :param paras: parameters to be simulated
        :param simulator: the simulator
        :param avg_result: the avg_result (a list of replications to simulation)
        :param plugin: the AIMSUN plugin for computing the average
        :param valid_data: validation data; if None, then will simply simulate and generate the validation
                            data. The obj_value in this case will be [0, 0]
        :param det_weight: the weight for detectors for computing the obj
        :param name: a string used for printing the status
        :param cmd_logger: logger for cmd output
        :return: [obj_value, sim_data];
                obj_value: [RMSE_speed, RMSE_count]
                sim_data: dict; [det_name] = [[speeds], [counts]]
        """
        if valid_data is not None:
            self.cmd_logger.print_cmd('\n------------Simulating with {0} paras-------------------'.format(name))
        else:
            self.cmd_logger.print_cmd('\n------------Generating validation data------------------')

        # set new parameters
        self.set_new_paras(experiment, paras)

        self.__simulate_experiment(simulator, avg_result)

        sim_data = self.__get_averaged_detector_data(avg_result, plugin)

        if valid_data is not None:
            obj_value = self.__evaluate_obj_val_RMSE(sim_data, valid_data, det_weight)
        else:
            obj_value = [0, [0]]

        return [obj_value, sim_data]

    # -------------------------- Utility functions for interfacing AIMSUN -------------------------------
    def __readParsFile(self, parsFilePath):
        """
        Read parameter files to create traffic states
        :param parsFilePath: the file name with full path.
            Each row of file: state name, car type and name, state start time, state duration
               e.g: cong_state1_car,53 Car,15:30:00,00:05:00
        :return: pars_dict: [state name] = [vehicle type, state start time, duration], all strings
        """
        pars_dict = OrderedDict()
        pars_file = open(parsFilePath, 'r')
        while True:
            line = pars_file.readline()
            if not (bool(line)):
                break
            line = line.strip()
            items = line.split(',')
            tmp_key = items[INDEX_PARS_STATE].strip()
            if tmp_key not in pars_dict.keys():
                pars_dict[tmp_key] = list()
            pars_dict[tmp_key].append(
                (items[INDEX_PARS_VEH].strip(), items[INDEX_PARS_FROMTIME].strip(), items[INDEX_PARS_DUR].strip()))
        pars_file.close()
        return pars_dict

    def __readFlowsFile(self, flowsFilePath):
        """
        Read flows file to create traffic states
        :param flowsFilePath: the file name with full path.
            Each row of file: state name, section id, inflow (veh/hr)
                       e.g: cong_state3_car,30109,50.0
        :return: flows_dict: [state name] = [section_id, inflow], all strings
        """
        flows_dict = OrderedDict()
        flows_file = open(flowsFilePath, 'r')
        while True:
            line = flows_file.readline()
            if not (bool(line)):
                break
            line = line.strip()
            items = line.split(',')
            tmp_key = items[INDEX_FLOWS_STATE].strip()
            if tmp_key not in flows_dict.keys():
                flows_dict[tmp_key] = list()
            flows_dict[tmp_key].append((items[INDEX_FLOWS_SECID].strip(), items[INDEX_FLOWS_INFLOW].strip()))
        flows_file.close()
        return flows_dict

    def __readTurnsFile(self, turnsFilePath):
        """
        Read turns file to create traffic states
        :param turnsFilePath: the file name with full path.
            Each row of file: state name, from section_id, to section_id, percent
               e.g: cong_state1_car,21217,23587,2.0
        :return: turns_dict: [state name] = [from_sec_id, to_sec_id, percent]
        """
        turns_dict = OrderedDict()
        turns_file = open(turnsFilePath, 'r')
        while True:
            line = turns_file.readline()
            if not (bool(line)):
                break
            line = line.strip()
            items = line.split(',')
            tmp_key = items[INDEX_TURNS_STATE].strip()
            if tmp_key not in turns_dict.keys():
                turns_dict[tmp_key] = list()
            turns_dict[tmp_key].append(
                (items[INDEX_TURNS_FROM].strip(), items[INDEX_TURNS_TO].strip(), items[INDEX_TURNS_PERC].strip()))
        turns_file.close()
        # print '\n\nturns_dict:{0}\n\n'.format(turns_dict)
        return turns_dict

    def __createState(self, state_name, pars_dict, flows_dict, turns_dict, main_entrance_id,
                      main_entrance_discount_ratio):
        """
        This function creates the traffic states from the dict read from files
        :param state_name: the statename to be created
        :param pars_dict: the dicts created by files
        :param flows_dict: the dicts created by files
        :param turns_dict: the dicts created by files
        :param main_entrance_id: the inflow to be discounted. Set as None if do not want any discount
        :param main_entrance_discount_ratio: set as 1 if not discount
        :return: return the State object
        """
        # create new state and set parameters
        state = GKSystem.getSystem().newObject("GKTrafficState", self.model)
        state.setName(state_name)
        it = QTime.fromString((pars_dict[state_name])[0][1], Qt.ISODate)
        duration = (pars_dict[state_name])[0][2]
        state.setInterval(it, GKTimeDuration.fromString(duration))

        if _debug is True:
            self.cmd_logger.print_cmd('AIMSUNFUNCTIONS: state from {0} duration {1}'.format(it.toString(), duration))

        # set the vehicle for this state
        vehicleString = str((pars_dict[state_name])[0][0]).split()
        vehId = int(vehicleString[0])  # make sure this is correct
        vehName = vehicleString[1]
        vehicle = state.getModel().getCatalog().find(vehId)
        if vehicle is None:
            # is wont work since the name is not unique, UserClass has object called car
            vehicle = state.getModel().getCatalog().findByName(vehName)
        state.setVehicle(vehicle)

        # set the inflow of the state
        for entrance in range(0, len(flows_dict[state_name])):
            # print (flows_dict[state_name])[entrance][0]
            fromSection = self.__findSection((flows_dict[state_name])[entrance][0])

            # discount the main entrance flow
            if fromSection.getId() == main_entrance_id:
                state.setEntranceFlow(fromSection, None,
                                      float((flows_dict[state_name])[entrance][1]) * main_entrance_discount_ratio)
            else:
                state.setEntranceFlow(fromSection, None,
                                      float((flows_dict[state_name])[entrance][1]))

        # set the turn percentage of the state
        if state_name in turns_dict.keys():
            for turn in range(0, len(turns_dict[state_name])):
                # print_cmd('For state {0}, has turns: {1}'.format(state_name, turns_dict[state_name]))
                fromSection = self.__findSection((turns_dict[state_name])[turn][0])
                toSection = self.__findSection((turns_dict[state_name])[turn][1])
                state.setTurningPercentage(fromSection, toSection, None, float((turns_dict[state_name])[turn][2]))
        else:
            self.cmd_logger.print_cmd('No Turn information for state {0}'.format(state_name))

        # for testing the aimsun automation
        if _debug:
            self.cmd_logger.print_cmdprint_cmd('AIMSUNFUNCTION: state turn 330->(340, 341): ({0},{1})'.format(
                state.getTurningPercentage(self.model.getCatalog().find(int(330)),
                                           self.model.getCatalog().find(int(340)), None),
                state.getTurningPercentage(self.model.getCatalog().find(int(330)),
                                           self.model.getCatalog().find(int(341)), None)))

        # append the state to the state folder
        folder = self.__getStateFolder()
        folder.append(state)

        return state

    def __findSection(self, entry):
        """
        This function returns a section object
        :param model: the GKModel
        :param entry: string, e.g., '330'
        :return: return GKSection object of the entry id
        """
        section = self.model.getCatalog().find(int(entry))
        if section.isA(QString("GKSection")) is False:
            section = None
        return section

    def __getStateFolder(self):
        """
        This function returns the folder object for the traffic states
        :param model: the GK Model
        :return:
        """
        folderName = "GKModel::trafficStates"
        folder = self.model.getCreateRootFolder().findFolder(folderName)
        if folder is None:
            folder = GKSystem.getSystem().createFolder(
                self.model.getCreateRootFolder(), folderName)

        # print_cmd('__getStateFolder: type: {0}'.format(type(folder))
        # print_cmd('__getStateFolder: name: {0}'.format(folder.getName())
        return folder

    def __getDemandFolder(self):
        """
        This function returns the demand folder
        :param model: GK model
        :return: Folder object
        """
        folderName = "GKModel::trafficDemand"
        folder = self.model.getCreateRootFolder().findFolder(folderName)
        if folder is None:
            folder = GKSystem.getSystem().createFolder(
                self.model.getCreateRootFolder(), folderName)
        return folder

    def __getScenariosFolder(self):
        """
        This function returns the folder object for the scenarios
        :param model: GK Model
        :return: return the folder
        """
        folderName = "GKModel::top::scenarios"
        folder = self.model.getCreateRootFolder().findFolder(folderName)
        if folder is None:
            folder = GKSystem.getSystem().createFolder(
                self.model.getCreateRootFolder(), folderName)
        return folder

    def __createScheduleItem(self, state):
        """
        Before the states can be added to the demand, each state must be added to a schedule item before being added to demand.
        NOTE: states contains the flow, turn percentage, and vehicle type. However, to add the state to the demand, we need to
        assign the state to a scheduleDemandItem. The simulator uses the fromTime and duration of the scheduleDemandItem, NOT
        the fromTime and duration set in the state!
        :param state: the state object
        :return: GKScheduleDemandItem object: associated with the state
        """
        schedule = GKScheduleDemandItem()
        schedule.setTrafficDemandItem(state)

        if _debug:
            self.cmd_logger.print_cmd(
                'schedule.state.duration {0}'.format(schedule.getTrafficDemandItem().getDuration().toString()))

        hr = state.getFrom().hour()
        minute = state.getFrom().minute()
        sec = state.getFrom().second()
        schedule.setFrom(3600 * hr + 60 * minute + sec)

        if _debug:
            self.cmd_logger.print_cmd('state.getDuration().toSeconds(): {0}'.format(state.getDuration().toSeconds()[0]))
            self.cmd_logger.print_cmd(type(state.getDuration().toSeconds()[0]))

        schedule.setDuration(state.getDuration().toSeconds()[0])

        return schedule

    def __simulate_experiment(self, simulator, avg_result):
        """
        This function simulates the replications in an experiment.
        :param simulator: the created simulator from create_simulator
        :param avg_result: the replication from setup_replication
        :return:
        """
        self.cmd_logger.print_cmd('\nReset replications...')

        # first reset replications
        avg_result.resetReplications()

        # add replications to simulator
        replication_list = avg_result.getReplications()

        for replication in replication_list:
            simulation_task = GKSimulationTask(replication, GKReplication.eBatch, "", "", True)  # Other approach
            simulator.addSimulationTask(simulation_task)
            self.cmd_logger.print_cmd('Added replication {0} to simulator with status {1}. '.format(replication.getId(),
                                                                                                    replication.getSimulationStatus()))
            self.cmd_logger.print_cmd(
                'pending {0}; done {1}; discarded {2}; loaded {3}'.format(GKGenericExperiment.ePending,
                                                                          GKGenericExperiment.eDone,
                                                                          GKGenericExperiment.eDiscarded,
                                                                          GKGenericExperiment.eLoaded))
        # simulate model
        if not simulator.isBusy():
            self.cmd_logger.print_cmd('Simulating...\n')
            sim_status = simulator.simulate()
        else:
            self.cmd_logger.print_cmd('Simulator is busy\n')

        # make sure correctly simulated
        if sim_status is True:
            self.cmd_logger.print_cmd('Simulation finished\n')
        else:
            self.cmd_logger.print_cmd('ERROR: Simulation failed\n')

    def __read_detector_data(self, data_origin):
        """
        This function reads the detector data from data_origin (a average_result type)
        :param model: GK Model
        :param data_origin: a list of replications:
                [a replication or the average_result(which is a subtype of GKReplication)]
        :return: a dict, avg_data[detector_name] = [[speed mph],[count/5min]]
        """
        avg_data = OrderedDict()

        det_type = self.model.getType("GKDetector")
        # read data for each replication and then the average
        for replication in data_origin:

            self.cmd_logger.print_cmd('\nReading Replication data: {0}'.format(replication.getName()))

            # get the column id
            speedColumn = det_type.getColumn(GK.BuildContents(GKColumnIds.eSpeed, replication, None))
            countColumn = det_type.getColumn(GK.BuildContents(GKColumnIds.eCount, replication, None))

            # read each detector
            for det in self.model.getCatalog().getObjectsByType(det_type).itervalues():

                det_name = str(det.getName())
                # print_cmd('----Reading Detector {0}...'.format(det_name))

                # add to dict
                # flow, speed
                avg_data[det_name] = [[], []]

                speedData = det.getDataValueTS(speedColumn)
                countData = det.getDataValueTS(countColumn)

                if countData.size() == 0 or speedData.size() == 0 or countData.size() != speedData.size():
                    self.cmd_logger.print_cmd('ERROR: Detector {0} has no data available'.format(det_name))
                else:
                    # print_cmd('----size of data is: {0}'.format(countData.size()))
                    # The speed data returned from AIMSUN is in km/h; 1 km/h = 0.62137 mph
                    for interval in range(countData.size()):
                        avg_data[det_name][0].append(speedData.getValue(GKTimeSerieIndex(interval))[0] * KMH2MPH)
                        avg_data[det_name][1].append(countData.getValue(GKTimeSerieIndex(interval))[0])

                        if _show_detector_data:
                            self.cmd_logger.print_cmd('--------interval {0}: speed {1}; count {2}'.format(interval,
                                                                                                          avg_data[
                                                                                                              det_name][
                                                                                                              0][-1],
                                                                                                          avg_data[
                                                                                                              det_name][
                                                                                                              1][-1]))
                            # print_cmd('----Detector {0} data:{1}'.format(det.getName(), avg_data[det.getName()])

        return avg_data

    def __get_averaged_detector_data(self, avg_result, plugin):
        """
        This function computes the average data from the simulation
        :param model: GK Model
        :param avg_result: avg_result replication
        :param plugin: the plugin from AIMSUN for computing the average
        :return: a dict, avg_data[detector_name] = [[speed mph],[count /5min]]
        """
        # compute the result
        calculate_status = plugin.calculateResult(avg_result)

        if calculate_status == GGetramModule.eOKCalculateResult:
            self.cmd_logger.print_cmd('Retrieving average data finished.')
        elif calculate_status == GGetramModule.eFailCalculateResult:
            # at 5th iteration failed.
            self.cmd_logger.print_cmd('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
            self.cmd_logger.print_cmd('$$$$$$$$$$$ ERROR: Retrieving average data failed.')
            self.cmd_logger.print_cmd('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
        elif calculate_status == GGetramModule.eSimRequiredCalculateResult:
            self.cmd_logger.print_cmd('$$$$$$$$$$$ ERROR: Retrieving average data failed. Simulation required.')

        # Not sure if this line is needed
        plugin.readResult(avg_result)

        # read the detector data out
        avg_data = self.__read_detector_data([avg_result])

        return avg_data

    def __evaluate_obj_val_RMSE(self, avg_data, valid_dict, det_weight):
        """
        This function evaluate the objective function using RMSE
        :param avg_data: the averaged simulation data speed (mph), count (every interval)
        :param valid_dict: the validation data
                    valid_dict[det_name] = [[speeds (mph)], [counts every interval]]
        :param det_weight: the weight assigned for each detector
        :return: RMSE_speed (mph), RMSE_all_speeds(mph)
                 the first returned value is a float which gives the average rmse of speeds for all detectors
                 the second returned value is a list of RMSE mph of detectors
        """
        rms_speed = 0
        # rms_count = 0
        rmse_speeds = []

        self.cmd_logger.print_cmd('\nEvaluating objective value...')

        # adds up the error for speed and count for each detector
        # Note the detector in the main entrance (EB3) is not used
        for key in valid_dict:
            valid_speed = np.array(valid_dict[key][0])
            # valid_count = np.array(valid_dict[key][1])
            avg_speed = np.array(avg_data[key][0])
            # avg_count = np.array(avg_data[key][1])

            # the following is added to deal with np.nan values
            # print 'before: {0}'.format(valid_speed-avg_speed)
            tmp_array = np.power(valid_speed - avg_speed, 2)
            # print 'after: {0}'.format(tmp_array)
            err = det_weight[key] * np.sqrt(np.nansum(tmp_array) / sum(~np.isnan(tmp_array)))
            rms_speed += err
            rmse_speeds.append(err)

            # tmp_array = np.power(valid_count - avg_count, 2)
            # rms_count += np.sqrt(np.nansum(tmp_array) / sum(~np.isnan(tmp_array)))

        self.cmd_logger.print_cmd('Evaluated objective: average RMSE speed: {0})'.format(rms_speed))

        return rms_speed, rmse_speeds


# ====================================================================================================
# ------------------- Class for cmd window logging and printing -----------------------
class CmdLogger:
    """
    This class logs the output in the cmd window to a text file
    """

    def __init__(self, path_str, start_time):
        """
        Start the logger
        :param path_str: the directory for saving the log file
        :param start_time: the start time which will be used for naming the file
        :return: a logger object
        """
        file_name = start_time.strftime("%Y%m%d_%H%M%S") + '_cmd_log'
        self.cmd_logger = open(path_str + file_name + '.txt', 'wb')

    def cmd_logger_header(self, description,
                          g_num_iter, g_num_rep, seed_list,
                          obj_func,
                          det_for_validation,
                          main_entrance_id, main_entrance_discount,
                          paras_list):
        self.print_cmd('Calibration experiment description:\n ---- {0}\n'.format(description))
        self.print_cmd('Calibration Configuration:')
        self.print_cmd(
            '---- Objective function is:           Minimize    {0}xRMS_speed + {1}xRMS_count'.format(obj_func[0],
                                                                                                     obj_func[1]))
        self.print_cmd('---- Detectors used for validation:   {0}'.format(det_for_validation))
        self.print_cmd('---- Number of iterations:            {0}'.format(g_num_iter))
        self.print_cmd('---- Main Entrance is:                {0}'.format(main_entrance_id))
        self.print_cmd('---- Main Entrance flow = VerMac EB3 x{0}'.format(main_entrance_discount))
        self.print_cmd('---- Number of replications {0}, with seeds: {1}\n'.format(g_num_rep, seed_list))

        self.print_cmd('\nParameters:')
        for para_tup in paras_list:
            self.print_cmd('-- {0} paras:'.format((para_tup[0])))
            for key in para_tup[1].keys():
                self.print_cmd('---- {0}: {1}'.format(key, para_tup[1][key]))

    def print_cmd(self, line_str):
        """
        This function writes the output into the cmd window as well as the text file
        :param line_str:
        :return:
        """
        self.cmd_logger.write(line_str + '\n')
        print line_str

    def stop(self):
        """
        Stop the logger and save the text file
        :return:
        """
        self.cmd_logger.close()

    def print_results(self, paras_list, obj_val_list):
        """
        This function prints out the result
        :param paras_list: a list of parameters [true_paras, default_paras]
        :param obj_val_list: [true_obj, default_obj], each obj is a tuple ('true', float_RMSE_mph)
        :param cmd_logger: the logger for cmd output
        :return:
        """
        self.print_cmd('\n\n===========================================================================')
        self.print_cmd('====================Calibration Finished===================================')
        # Forget about beautiful printout. Just log information
        self.print_cmd('Parameters: \n')
        for para in paras_list:
            self.print_cmd('---- {0}_paras:'.format(para[0]))
            for key in para[1].keys():
                self.print_cmd('-------- {0}:    {1}'.format(key, para[1][key]))

        self.print_cmd('\nObjective values: \n')
        for obj_val in obj_val_list:
            self.print_cmd('{0} objective value:\n---- RMS_speed: {1}'.format(obj_val[0], obj_val[1]))

    def print_opt_steps(self, opt_steps, baseline):
        """
        This function prints out the optimization steps
        :param opt_steps: the list updated by Opt_Logger
        :param baseline: the objective function value of the baseline
        :param cmd_logger: logger for cmd output
        :return:
        """
        self.print_cmd('Optimization steps: {0}'.format(opt_steps))
        self.print_cmd('        Optimal by {0} steps: {1} mph, default: {1} mph'.format(len(opt_steps),
                                                                                                   np.min(opt_steps),
                                                                                                   baseline))


# ------------------- Class for handling optimization logger -----------------------
class Opt:
    """
    This class:
        - Handles the communication with the optimization
        - logs the optimization steps into a text file
    """

    def __init__(self, path_str, start_time):
        """
        Start the logger
        :param path_str: string of the dir for saving the log file
        :param start_time: the start time of the optimization which is used as the log file name
        :return: a logger object
        """
        file_name = start_time.strftime("%Y%m%d_%H%M%S") + '_opt_log'
        self.opt_logger_file = open(path_str + file_name + '.csv', 'wb')

    def log_opt_step(self, opt_steps, iteration_counter, paras, obj_values):
        """
        This function logs the optimization step into the log file
            iteration_counter; car_para1; car_para2..., truck_para1, truck_para2; RMS_speed_avg; RMS_speeds
        :param opt_steps: a list for keeping a copy of optimal objective function values
        :param iteration_counter: int, iteration counter.
        :param paras: the parameters simulated in this iteration.
        :param obj_values: avg_speed_error (mph), [speed_error (mph) for each detector]
        :return:
        """
        string = str(iteration_counter) + ';'

        # first append key, then parameters
        for key in paras.keys():
            # add the key
            string += key + ':'
            string += ','.join([str(v) for v in paras[key]])
            string += ';'

        # Last two are reserved for the RMS result
        string += 'RMSE_avg:' + str(obj_values[0]) + ';'
        string += 'RMSEs:' + ','.join([str(v) for v in obj_values[1]])

        # only append the objective value
        opt_steps.append(obj_values[0])

        # write in file
        self.opt_logger_file.write(string + '\n')

    def stop(self):
        """
        This function stops the logger and saves the file.
        :return:
        """
        self.opt_logger_file.close()

    # ------------------- Functions for interfacing with optimization solver -----------------------
    def read_opt_paras(self, parafile, cmd_logger):
        """
        This function reads the new parameters.
        Note: MODIFY set_new_paras if the keys changed.
        :param parafile: the parameter file, each line:
            parameter_name,value1,value2,value3
        :param cmd_logger: logger for cmd output
        :return: dict. [parameter_name]:[value1, value2, value3]
        """
        paras = OrderedDict()

        wait_time = 0
        timeout = 0
        while not exists(parafile):
            time.sleep(0.1)
            wait_time += 1
            timeout += 1
            if wait_time >= 10:  # sleep 1 second
                cmd_logger.print_cmd('Waiting for paras...')
                wait_time = 0

            if timeout >= 200:  # 20 s
                # assume finished optimization
                return None

        if exists(parafile):
            paras = self.__paras_reader(parafile)

        # delete the file once read
        os.remove(parafile)

        # print_cmd('Have read paras:\n')
        # for key in paras.keys():
        #     print_cmd('---- {0}: {1}'.format(key, paras[key]))

        return paras

    def read_opt_solution(self, solutionfile, cmd_logger):
        """
        This function reads the final solution from optquest.
        The solver will stop at its max iteration.
        :param solutionfile: same format as the parafile
        :param cmd_logger: logger for cmd output
        :return: dict. same as read_optquest_paras
        """
        # Have not finished optimization
        if not exists(solutionfile):
            return None

        # a solution has been converged
        else:

            solution = self.__paras_reader(solutionfile)

            # delete the file once read
            # os.remove(solutionfile)

            cmd_logger.print_cmd('Have read solution:\n')
            for key in solution.keys():
                cmd_logger.print_cmd('---- {0}: {1}'.format(key, solution[key]))

            return solution

    def write_simval(self, simval, simvalfile, cmd_logger):
        """
        This function writes the simulation result to file
        :param simval: double, the optimal value
        :param simvalfile: the file string.
        :param cmd_logger: logger for cmd output
        :return:
        """
        f = open(simvalfile, 'w')
        f.write(str.format("{0:.16f}", simval))
        f.close()
        cmd_logger.print_cmd('Wrote objective value: {0}'.format(simval))

        return 0

    def __paras_reader(self, parafile):
        """
        A universal paras reader assuming the file exist.
        :param parafile: the file to read
        :return: a dict, key-value
        """
        paras = OrderedDict()

        f = open(parafile, 'r')

        for line in f:
            line = line.strip()
            items = line.split(',')

            # first item is the key
            paras[items[0]] = []

            for i in range(1, len(items)):
                paras[items[0]].append(float(items[i]))

        f.close()

        return paras

    # ------------------- Functions for loading and processing optimization solutions ---------------
    def load_previous_opt_solutions(self, opt_step_file):
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

    def try_get_obj_from_previous_solutions(self, solution_list, new_para):
        """
        Check if the new_para has been previously simulated.
        :param solution_list: the paras_list from load_previous_opt_solutions
        :param new_para: the new parameters
        :return:
        """

        # if not solution list
        if solution_list is None:
            return None

        for i in range(0, len(solution_list[0])):

            if new_para == solution_list[0][i]:
                return solution_list[1][i]

        # if none found
        return None

    def save_solution_data(self, solution, data, path_name, start_time, name):
        """
        This function saves the optimization result (parameters and detector data) to a file with standard format.
            First row of file is the parameters, then
            detector 1 speed (mph)
            detector 2 count
            ...
        :param solution: the parameters
        :param data: the detector data
        :param path_name: the path to save to
        :param start_time: the start time of the simulation
        :param name: the name of the set of parameters
        :return:
        """
        file_name = start_time.strftime("%Y%m%d_%H%M%S") + '_sol_'
        f = open(path_name + file_name + name + '.csv', 'wb')
        writer = csv.writer(f)
        # the solution: key1, value, key2, value
        if solution is not None:
            list = []
            for key in solution.keys():
                list.append(key)
                for item in solution[key]:
                    list.append(str(item))
            writer.writerow(list)
        else:
            # just to save the validation data. true parameters are unknown
            writer.writerow(['Data saved with unknown parameters: {0}'.format(name)])

        for key in data:
            tmp_line = []
            tmp_line.append(key)
            # write speed
            for item in data[key][0]:
                tmp_line.append(item)
            writer.writerow(tmp_line)

            tmp_line = []
            tmp_line.append(key)
            # write count
            for item in data[key][1]:
                tmp_line.append(item)
            writer.writerow(tmp_line)
            # writer.writerow([str(key), str(data[key][0]), str(data[key][1])])

        f.close()


# -------------------------- Other utility functions -------------------------------
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

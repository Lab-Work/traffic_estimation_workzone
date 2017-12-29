# IDOT-SmartWorkzone
Yanning Li, Juan Carlos Martinez Mori, August, 2016

## 1) Overview
This repository contains the source code developed during the project *"Improving the effectiveness of smart work zones"* (IDOT R27-155). This repository also supports a journal paper *"Estimating traffic conditions from smart work zone systems"* was submitted to *Journal of Intelligent Transportation Systems* during this project. The preprint of the journal can be downloaded in this [link](https://www.dropbox.com/s/0p4s5amhcjjou5h/LiMoriWork2016.pdf?dl=0). The final project report will be available upon the completion of the project by Oct 15, 2016.

In this project, two work zones in Illinois were modeled and calibrated in a micro calibration software AIMSUN. The source code for calibrating AIMSUN using its API and the calibration data can be found in the AutoCalibration folder. **The calibrated micro parameters can be found in the files:[I80_calibrated_AIMSUN_parameters.pdf](https://github.com/Lab-Work/IDOT-SmartWorkzone/blob/master/I80_calibrated_AIMSUN_parameters.pdf) and [I57_calibrated_AIMSUN_parameters.pdf](https://github.com/Lab-Work/IDOT-SmartWorkzone/blob/master/I80_calibrated_AIMSUN_parameters.pdf).**

In the calibrated work zones, a variety of sensor network configurations and traffic estimation algorithms were evaluated in terms of traffic estimation accuracy. The source code and data supporting the cross evaluation can be found in the [Cross_Evaluation folder](https://github.com/Lab-Work/IDOT-SmartWorkzone/tree/master/Cross_Evaluation). 

## 2) License

This software is licensed under the *University of Illinois/NCSA Open Source License*:

**Copyright (c) 2016 The Board of Trustees of the University of Illinois. All rights reserved**

**Developed by: Department of Civil and Environmental Engineering University of Illinois at Urbana-Champaign**

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal with the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimers. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimers in the documentation and/or other materials provided with the distribution. Neither the names of the Department of Civil and Environmental Engineering, the University of Illinois at Urbana-Champaign, nor the names of its contributors may be used to endorse or promote products derived from this Software without specific prior written permission.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS WITH THE SOFTWARE.

## 3) Folders
This repository consists of mainly two separate componets: autocalibration and cross-evaluation. It is organized as follows:

### [./I80_calibrated_AIMSUN_parameters.pdf](https://github.com/Lab-Work/IDOT-SmartWorkzone/blob/master/I80_calibrated_AIMSUN_parameters.pdf)
This file contains the AIMSUN parameters calibrated for the I-80 work zone.

### [./I57_calibrated_AIMSUN_parameters.pdf](https://github.com/Lab-Work/IDOT-SmartWorkzone/blob/master/I80_calibrated_AIMSUN_parameters.pdf)
This file contains the AIMSUN parameters calibrated for the I-80 work zone.

### [./AutoCalibration/](https://github.com/Lab-Work/IDOT-SmartWorkzone/tree/master/AutoCalibration)
This folder contains the source code and data for automated calibration of AIMSUN parameters using a nonlinear optimization program NOMAD.

### [./Cross_Evaluation](https://github.com/Lab-Work/IDOT-SmartWorkzone/tree/master/Cross_Evaluation)
This folder contains the source code for the generation of a variety of sensor network configurations, the implementation of three traffic estimation algorithms (spatial interpolation, spatio-temporal filtering, and ensemble Kalman filter), visualization and evaluation of the estimation results. 
- [./Cross_Evaluation/Ixx_topology_####.txt](https://github.com/Lab-Work/IDOT-SmartWorkzone/tree/master/Cross_Evaluation). 
These are the text files containing the network topology which is required for extracting the virtual sensor data and true states. 
- [./Cross_Evaluation/Ixx_configurations_input.txt](https://github.com/Lab-Work/IDOT-SmartWorkzone/tree/master/Cross_Evaluation). These are the configuration file used for the estimators. Each configuration file includes blocks separated by a blank line. Each block describes one sensor network configuration (with sensor parameters), and multiple algorithms (with parameters) to be applied to the sensor network configuration. The main program will generate estimation results for all configurations.
- [./Cross_Evaluation/src/](https://github.com/Lab-Work/IDOT-SmartWorkzone/tree/master/Cross_Evaluation/src). This folder contains all the source code for the cross evaluation of sensor networks and algorthms. See next section for instructions on how to run the code.
- [./Cross_Evaluation/Trajectory_data/](https://github.com/Lab-Work/IDOT-SmartWorkzone/tree/master/Cross_Evaluation/Trajectory_data). This folder contains the simulated trajectory data extracted from AIMSUN sqlite database (table MIVEHDETAILEDTRAJECTORY), which is used for generating realistic sensor measurements and true states. The complete trajectory data can be downloaded from this [link](https://uofi.box.com/s/2bkwejveuew2fospxcn4pqi0df7wrrfj). 
- [./Cross_Evaluation/Virtual_sensor_data/](https://github.com/Lab-Work/IDOT-SmartWorkzone/tree/master/Cross_Evaluation/Virtual_sensor_data). This folder contains the virtual sensor data generated for two work zones. 
  * [./Cross_Evaluation/Virtual_sensor_data/Ixx_rep####_to_generate.txt](https://github.com/Lab-Work/IDOT-SmartWorkzone/tree/master/Cross_Evaluation/Virtual_sensor_data). These files contains the parameters for sensors to be generated. Each row desribes the information for each sensor.
  * [./Cross_Evaluation/Virtual_sensor_data/Ixx_rep####_generated.txt](https://github.com/Lab-Work/IDOT-SmartWorkzone/tree/master/Cross_Evaluation/Virtual_sensor_data). These are text files that logs the generated virtual sensors to aviod re-generating data.
  * [./Cross_Evaluation/Virtual_sensor_data/I80/](https://github.com/Lab-Work/IDOT-SmartWorkzone/tree/master/Cross_Evaluation/Virtual_sensor_data/I80). This folder contains all the virtual sensor data generated for different replications in I-80 work zone. In each subfolder (rep####), each text (e.g. *RADARLoc700m.txt*) file saves the virtual sensor data for the sensor (e.g. sensor type *radar* located with absolute distance *700 m* to the entrance of the modeled section). Data format in each file is: timestamps (s), speed (kph), count.
  * [./Cross_Evaluation/Virtual_sensor_data/I57/](https://github.com/Lab-Work/IDOT-SmartWorkzone/tree/master/Cross_Evaluation/Virtual_sensor_data/I57). This folder contains all the virtual sensor data for I-57 work zone. 
- [./Cross_Evaluation/Estimation_results/](https://github.com/Lab-Work/IDOT-SmartWorkzone/tree/master/Cross_Evaluation/Estimation_results). This folder contains the estiamtion resutls. The complete estimation results can be downloaded from this [link](https://uofi.box.com/s/g6j15tbc2nlnftd0hgv87r2eru6xon6b). 
  * [./Cross_Evaluation/Estimation_results/Ixx_rep###_res5s##m_generated.txt](https://github.com/Lab-Work/IDOT-SmartWorkzone/tree/master/Cross_Evaluation/Estimation_results). These are the text files that log the genenrated estimation results. E.g., the file *I80\_rep41368\_res5s50m\_generated.txt* logs the generated estimation results for configurations and algorithms for replication *41368* in *I-80* work zone at a resolution *5 s* and *50 m*. It should be noted that true states are all generated at resolution *5s50m*; the spatial interpolation and spatio-temporal filtering are applied to grid *5s50m* while the ensemble Kalman filter is applied to *5s200m* to ensure numerical stability. For final cross comparison, the ensemble Kalman filter results should be copied to 5s50m\_generated.txt and 5s50m folder. Here are the general structure of the folder.
  * [./Cross_Evaluation/Estimation_results/I80/rep41368/5s50m/](https://github.com/Lab-Work/IDOT-SmartWorkzone/tree/master/Cross_Evaluation/Estimation_results/I80/rep41368/5s50m). This folder contains the estimation results for a specific sensor network configuration and algorithm. All unites in this folder are in **m/s**(speed), **veh/m**(density), **m**(queue), **s**(travel time). E.g., the file *configHOMO_11_RADAR_linearFILL_speed.txt* is the velocity field estimation in **m/s** with 11 radars deployed and using linear spatial interpolation algorithm. The file contains a matrix with num\_steps(row) x num\_space(col) with the left top as the origin. The estimation results for ensemble Kalman filter (at 5s 200 m resolution) should be copied to this folder for final cross evaluation. 
  * [./Cross_Evaluation/Estimation_results/I80/rep41368/truestate/](https://github.com/Lab-Work/IDOT-SmartWorkzone/tree/master/Cross_Evaluation/Estimation_results/I80/rep41368/truestate). This folder contains the true states. All unites in this folder are in **m/s**(speed), **veh/m**(density), **m**(queue), **s**(travel time).
      + The file *xx_density.txt* gives the true density field in veh/m.
      + The file *xx_speed.txt* gives the true speed field in m/s.
      + The file *xx_queue.txt* gives the queue length in meters.
      + The file *xx_true_traveltime.txt* gives the true travel time in seconds. The true travel time at time t is defined as the travel time for the vehicle that **enters** the road at time t.
      + The file *xx_measured_traveltime.txt* gives the measured travel time (s) by Bluetooth sensors with 100% penetration rate. In this file, the travel time at time t is the measured travel time of the vehicle **exits** at time t.
      + The file *xx_trueinst_traveltime.txt* gives the instantaneous travel time (s) computed from the true speed field. 


## 4) Run the code

### Run autocalibration of AIMSUN
The autocalibration has been configured to run automatically by running `./AutoCalibration/*.bat` files. For instance, to run autocalibration of I80 work zone, simply
- double click `./AutoCalibration/I80_thread1_start_NOMAD.bat`, which will start NOMAD optimization solver and show status in the terminal.
- double click `./AutoCalibration/I80_thread1_start_simulation.bat`, which will start AIMSUN simulation and show status in the terminal.

AIMSUN supports running two instances of simulation (or more depending on the type of th license). Hence you may use another thread to optimize parameters in parallel by double clicking the bat files `./AutoCalibration/I80_thread2_start_NOMAD.bat` and `./AutoCalibration/I80_thread1_start_simulation.bat`. 

Please refer to the README.md file in the `./AutoCalibration/` folder to learn how to change the optimization parameters. 

### Run cross evaluation
Cross evaluation mainly consists of the following steps:

1. `cd` to the folder ./Cross\_Evaluation/src/.
2. Generate virtual sensor data by `python I80_generate_sensors.py`. Comment/uncomment sections in the source code to generate the virtual sensor data desired. 
3. Generate a configuration file by `python I80_generate_configurations.py`. Comment/uncomment sections in the source code to generate the desired configuration input file [./Cross_Evaluation/Ixx_configurations_input.txt](https://github.com/Lab-Work/IDOT-SmartWorkzone/tree/master/Cross_Evaluation). 
4. Run estimators by `python I80_run_cross_evaluation.py`. Comment/uncomment sections as needed. The section `cross_eval.run_estimators()` runs all the estimators for the configurations specified in the configuration input file. If the estimation result for a combination of a sensor network and algorithm has been previoiusly generated and logged in [./Cross_Evaluation/Estimation_results/Ixx_rep###_res5s##m_generated.txt](https://github.com/Lab-Work/IDOT-SmartWorkzone/tree/master/Cross_Evaluation/Estimation_results), the it will not be re-generated. Remove entries in [./Cross_Evaluation/Estimation_results/Ixx_rep###_res5s##m_generated.txt](https://github.com/Lab-Work/IDOT-SmartWorkzone/tree/master/Cross_Evaluation/Estimation_results) to regenerate estimation results for certain sensor networks and algorithms.
5. Visualize and compare estimation resutls by `python I80_run_cross_evaluation.py`. Comment/uncomment sections to visualize and compare the estimation results. 


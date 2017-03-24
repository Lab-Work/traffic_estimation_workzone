This folder contains the source code and data for autocalibration of microscopic traffic models using a state-of-art simulation based optimization toolbox NOMAD ().


The current version of AIMSUN supports running two instances simutaneously. Therefore, two copies of source codes and files exists in the folder marked by thread1 and thread2.
To start the autocalibration for each thread:
1. Generate the demand data in folder (./data/Ixx/demand_data). The demond data are obtained from the sensor located at the entrance. 
	- Use the script generate_demand.py to facilitate the generation of data in the requried data format. 
	- Use the script convet_to_OD.py to convert to the matrix format which can be loaded to the aimsun ang file.
2. Start NOMAD:
	- Go to folder ./src/nomad/Ixx/ 
	- Configure the starting point in file thread1_I80_8dim_param_with_con.txt.
	- If modified the cpp file, compile the source code in VS command prompt: cl.exe *.cpp /EHsc
	- Make sure ./Ixx_thread1_start_NOMAD.bat points to the right file. Double click it to run the optimizor. 
3. Start AIMSUN:
	- Configure the calibration parameters in Ixx_configuration.txt file.
	- Make sure ./Ixx_thread1_start_simulation.bat points to the right file. Double clike it to run the simulation.
4. Check status of optimization.
	- NOMAD and AIMSUN communicates in the folder E:\\sim_opt_com\\ by files thread1_sim_paras.txt, thread1_sim_val.txt. Check this folder to make sure they are correctly communicating.
	- The log file and optimization steps will be saved in ./data/Ixx/Logs/thread#/, the files are named by the starting time of the calibration.
	- The folder in E:\\workzone_tmp\\ is used for NOMAD to cache intermediate data.
	- The command windows outputs the current step of the optimization and the current optimal parameters.



Organization:

/src:
	./aimsun_api/:
		This folder contains the AIMSUN API for running the simulation in AIMSUN. 

	./generate_demand:
		This folder contains the freeway entrance data and manually calibrated on/off ramp flow data. Run the code to convert the data into OD matrix which can be loaded to AIMSUN.

	./nomad:
		This folder contains the C++ source code and configuration file for NOMAD optimizer.

	./sensitivity_analysis:
		This folder contains the source code for analyzing the sensitivty of parameters in a synthetic work zone I00.

	./deprecated:
		This folder contains the deprecated code, including another optimizor OptQuest.
		
/data:
	./I80:
		./demand_data:
			This folder contains the source code, the inflow data (collected by sensors) at the entrace, and the manually calibrated on/off ramp flow data. Use the code to convert the flow files to OD matrix which can be loaded to AIMSUN.

		./validation_data:
			This folder contains the validataion data collected by sensors, which is used for computing the RMSE of the simulated traffic compared to the true traffic measurement. 

		./aimsun_files:
			This folder contains the original I80_EB_default.ang files, and the calirbated final optimal aimsun file.

		./Logs:
			This folder contains the logs generated during the calibration process.









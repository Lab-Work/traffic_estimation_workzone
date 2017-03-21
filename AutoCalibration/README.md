This folder contains the source code and data for autocalibration of microscopic traffic models using a state-of-art simulation based optimization toolbox NOMAD ().


Run the calibration:
0. Generate the demand data in folder.
1. Configure the calibration parameters in file. Specify the bounds and resolution.
2. Edit file run_autocalibration.bat to specify the correct workzone and threads
3. Double click run_autocalibration.bat to run the autocalibration.



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









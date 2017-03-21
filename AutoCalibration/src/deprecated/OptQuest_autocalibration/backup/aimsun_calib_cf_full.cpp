// Yanning Li
// July 19, 2015
// This is a major revision from previous version.
// This code reads a configuration file, and then decide the decision variables and number of iteration correspondingly.
// This version of the code hasn't been finished. 
// We can configure in a file and make the solver universal. However, it takes too much effort, since those parameters can 
// not be easily passed to the evaluateSolution function.



#include "stdafx.h"
#include "aimsun_calib_cf_full.h"
#include "OptQuestClient.h"
#include "DecisionVariable.h"
#include "Constraint.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

COptQuestClient oq;


using namespace std;

int _tmain(int argc, TCHAR* argv[], TCHAR* envp[])
{
	std::string sol_file ("C:\\sim_com\\sim_sol.txt");
	std::string sol_5_file ("C:\\sim_com\\sim_sol_5.txt");
	
	
	std::string config_file ("C:\\sim_com\\optquest_config.txt");

	// some constants
	const int max_buffer_size = 12;
	const int max_num_paras = 100;

	// some counters
	int num_paras = 0;
	int para_counter = 0;
	int item_counter = 0;

	// buffer saves the configuration file reading
	int paras_len[max_num_paras];	// saves the lenth of each line
	string paras_name[max_num_paras];
	double paras_min[max_num_paras];
	double paras_max[max_num_paras];
	double paras_step[max_num_paras];
	string paras_output_format[max_num_paras];




	HRESULT hr = S_OK;

    hr = CoInitialize(0);
    ATLASSERT(SUCCEEDED(hr));

	oq.Init();
	oq.SetLicenseID(994109469);
	
	// read configuration
	std::ifstream config (config_file.c_str());
	std::string line;

	// config.open(config_file.c_str());

	if (!config.is_open())
	{
		printf("configuration file missing...\n");
		return -1;
	}
	else
	{
		para_counter = 0;
		while (std::getline(config, line))
		{   
			para_counter += 1;
			std::istringstream ss(line);
			// parse the line, [name, min, max, step, (output format, replace -1 with value)]
			std::string item;
			item_counter = 0;
			std::string format_str;
			while (std::getline(ss, item, ','))
			{
				item_counter += 1;
				if (item_counter == 1)
					{paras_name[para_counter-1] = item;   cout << paras_name[para_counter-1] << '\n' ; }
				else if (item_counter == 2)
					{paras_min[para_counter-1] = std::stod(item); cout << paras_min[para_counter-1] << '\n' ;}
				else if (item_counter == 3)
					{paras_max[para_counter-1] = std::stod(item); cout << paras_max[para_counter-1] << '\n' ;}
				else if (item_counter == 4)
					{paras_step[para_counter-1] = std::stod(item); cout << paras_step[para_counter-1] << '\n' ;}
				else if (item_counter == 5)
				{
					format_str.clear();
					format_str = format_str + item;
				}
				else 
				{
					format_str = format_str + ',' + item;
				}
			}
			paras_output_format[para_counter-1] = format_str;
			cout << paras_output_format[para_counter-1] << '\n';
			paras_len[para_counter-1] = item_counter;
		}

	}

	config.close();
	printf("configuration file read\n");

	Sleep(10000);	


	try
	{
		// ProblemSetup holds the OptQuest identifiers
		ProblemSetup setup; 
		oq.setup = &setup;

		printf("to add variables\n");
		// a discrete variable for the car reaction time
		DecisionVariable X1;
		X1.m_name = "x1";
		// This is the car reaction time, must be an integral times of simulation step (0.8s)
		X1.m_ID = oq.AddDiscreteVariable(X1.m_name, 0.8, 2.4, 0.8);
		// X1.m_ID = oq.AddContinuousVariable(X1.m_name, 0.9, 1.4);
		setup.m_decVars.Add(X1);
		printf("Added variable x1\n");

		// a discrete variable for the max_accel_mean of the car
		DecisionVariable X2;
		X2.m_name = "x2";
		X2.m_ID = oq.AddDiscreteVariable(X2.m_name, 2.6, 3.4, 0.01);
		// X2.m_ID = oq.AddContinuousVariable(X2.m_name, 0, 0.3);
		setup.m_decVars.Add(X2);
		printf("Added variable x2\n");

		// a discrete variable for the onramp 3 flow. Was 500 cars/hr
		DecisionVariable X3;
		X3.m_name = "x3";
		X3.m_ID = oq.AddDiscreteVariable(X3.m_name, 0.0, 500.0, 25.0);
		setup.m_decVars.Add(X3);
		printf("Added variable x3\n");

		// a discrete variable for the offramp 3 ratio. Was 0
		DecisionVariable X4;
		X4.m_name = "x4";
		X4.m_ID = oq.AddDiscreteVariable(X4.m_name, 0, 20, 1);
		setup.m_decVars.Add(X4);
		printf("Added variable x4\n");

		// a discrete variable for the offramp 4 ratio. Was 3
		DecisionVariable X5;
		X5.m_name = "x5";
		X5.m_ID = oq.AddDiscreteVariable(X5.m_name, 0, 20, 1);
		setup.m_decVars.Add(X5);
		printf("Added variable x5\n");

		// a discrete variable for the truck reaction time
		DecisionVariable Y1;
		Y1.m_name = "y1";
		Y1.m_ID = oq.AddDiscreteVariable(Y1.m_name, 0.8, 2.4, 0.8);
		// Y1.m_ID = oq.AddContinuousVariable(Y1.m_name, 0.7, 1.1);
		setup.m_decVars.Add(Y1);
		printf("Added variable y1\n");

		// a discrete variable for the truck max accel mean
		DecisionVariable Y2;
		Y2.m_name = "y2";
		Y2.m_ID = oq.AddDiscreteVariable(Y2.m_name, 0.6, 1.8, 0.01);
		//Y2.m_ID = oq.AddContinuousVariable(Y2.m_name, 0, 0.3);
		setup.m_decVars.Add(Y2);
		printf("Added variable y2\n");

		// Make sure they are put in the correct order
		printf("Added variables: %S, %S, %S, %S, %S, %S, %S \n", setup.m_decVars.GetAt(0).m_name,
			setup.m_decVars.GetAt(1).m_name, setup.m_decVars.GetAt(2).m_name, setup.m_decVars.GetAt(3).m_name,
			setup.m_decVars.GetAt(4).m_name, setup.m_decVars.GetAt(5).m_name, setup.m_decVars.GetAt(6).m_name);

		//second argument is setMaximize. If true, then maximize; if false, then minimize
		setup.m_objID = oq.AddUserControlledObjective("myObj", false);
		
		oq.SetMaximumIterations(2);

		// Just modify SolutionEvaluator to communicate with the AIMSUN python script and then output a single objective function value

		// The problem has been setup.  Call Optimize() to start the optimization.
		// See COptQuestClient::EvaluateEvent(). This method is called when a solution
		// is ready to be evaluated.
		// See COptQuestClient::MonitorStatusEvent().  This method is called when a solution
		// has completed evaluation.  If you are giving the user feedback on the progress of the
		// optimization, do it from the MonitorStatusEvent() method.
		// If there are errors in the problem setup, Optimize() will throw an exception.
		oq.Optimize();
		char buffer[_CVTBUFSIZE];
		char iterBuf[15];

		// The optimization completed without throwing an exception.  Get the termination reason. 
		long reason = oq.GetTerminationReason();
		printf("Optimization completed");
		_itoa_s(reason, iterBuf, 15,10);
		printf("Termination reason = ");
		/* Returns the reason the optimization stopped.  
		* 0 = The optimization has not started
		* 1 = The optimization is still running	 
        * 3 = The optimization was solved using an Linear/Integer/Mixed Integer Program.
		* 4 = The optimization stopped due to the Auto Stop feature
		* 5 = All solutions have been generated.
		* 6 = The optimization stopped when the maximum number of iterations was reached.
		* 7 = The optimization stopped when the maximum time was reached.
		* 8 = The optimization was stopped by the user.
		* 10 = The optimization stopped due to an exception.
		* 11 = The Meta-Heuristic factories have completed.
		* 12 = There are no solutions that satisfy the constraints.
		* 13 = New (different) solutions cannot be generated.
		*/
		printf(iterBuf);
		printf("\n\n");

		// Output results - Print the top 5 solutions
		double objValue;
		double LHS;
		double RHS;

		// The following are getting the best a few solutions
		// Print out the best solutions (parameters) OptQuest found for AIMSUN.
		for (int nth = 1; nth <= 5; nth++)
		{
			int solutionID = oq.GetNthBestSolution(nth);
			objValue = oq.GetSolutionObjectiveValue(solutionID, setup.m_objID);

			_gcvt_s(buffer, _CVTBUFSIZE, objValue, 12);
			_itoa_s(solutionID, iterBuf,15,10);
			printf( "Objective value = ");
			printf(buffer);
			printf(" at iteration ");
			printf(iterBuf);
			printf(" with parameters:\n ");
			for (int ivar = 0; ivar < setup.m_decVars.GetCount(); ivar++)
			{
				double varValue = oq.GetSolutionVariableValue(solutionID, setup.m_decVars.GetAt(ivar).m_ID);
				CString name = setup.m_decVars.GetAt(ivar).m_name;
				printf("%S, %f;\n", name.GetString(), varValue);
			}
			printf ("\n");

			// if nth == 1, save into a file
			// Now write the solution parameters to the sim_sol.txt file
			if (nth == 1){
				std::ofstream solfile;
				solfile.open(sol_file.c_str());

				double x1_sol = oq.GetSolutionVariableValue(solutionID, setup.m_decVars.GetAt(0).m_ID);
				double x2_sol = oq.GetSolutionVariableValue(solutionID, setup.m_decVars.GetAt(1).m_ID);
				double x3_sol = oq.GetSolutionVariableValue(solutionID, setup.m_decVars.GetAt(2).m_ID);
				double x4_sol = oq.GetSolutionVariableValue(solutionID, setup.m_decVars.GetAt(3).m_ID);
				double x5_sol = oq.GetSolutionVariableValue(solutionID, setup.m_decVars.GetAt(4).m_ID);
				double y1_sol = oq.GetSolutionVariableValue(solutionID, setup.m_decVars.GetAt(5).m_ID);
				double y2_sol = oq.GetSolutionVariableValue(solutionID, setup.m_decVars.GetAt(6).m_ID);


				if (solfile.is_open()){
					solfile << (x1_sol) << "," << (x2_sol) << "," << (x3_sol) << "," << (x4_sol) << "," << (x5_sol) <<"\n" ;
					solfile << (y1_sol) << "," << (y2_sol) << "\n" ;}
				else 
					printf("unable to open C:\\sim_com\\sim_paras.txt file\n");
					// {printf("unable to open C:/paras.txt \n");}
				solfile.close();
			}


		}
		char c = cin.get();


		

	}
	catch (long err)
	{
		// TODO: map error to something meaningful for the user.
		char buffer[15];
		_itoa_s(err, buffer,15,10);
		printf("Exception number = ");
		printf(buffer);
		char c = cin.get();

	}
	CoUninitialize();

	return 0;
}


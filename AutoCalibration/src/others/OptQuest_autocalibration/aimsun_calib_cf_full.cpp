// Yanning Li
// July 19, 2015
// This is a major revision from previous version.
// This code reads a configuration file, and then decide the decision variables and number of iteration correspondingly.

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

	HRESULT hr = S_OK;

    hr = CoInitialize(0);
    ATLASSERT(SUCCEEDED(hr));

	oq.Init();
	oq.SetLicenseID(994109469);

	try
	{
		// ProblemSetup holds the OptQuest identifiers
		ProblemSetup setup; 
		oq.setup = &setup;

		printf("to add variables\n");

		// a discrete variable for the car speed acceptance
		DecisionVariable X0;
		X0.m_name = "x0";
		X0.m_ID = oq.AddDiscreteVariable(X0.m_name, 0.85, 1.3, 0.01);
		setup.m_decVars.Add(X0);
		printf("Added variable x0\n");

		// a discrete variable for the car reaction time
		DecisionVariable X1;
		X1.m_name = "x1";
		// This is the car reaction time, must be an integral times of simulation step (0.8s)
		X1.m_ID = oq.AddDiscreteVariable(X1.m_name, 0.8, 1.6, 0.8);
		setup.m_decVars.Add(X1);
		printf("Added variable x1\n");

		// a discrete variable for the max_accel_mean of the car
		DecisionVariable X2;
		X2.m_name = "x2";
		X2.m_ID = oq.AddDiscreteVariable(X2.m_name, 2.6, 3.4, 0.02);
		setup.m_decVars.Add(X2);
		printf("Added variable x2\n");

		// a discrete variable for car minDist mean
		DecisionVariable X3;
		X3.m_name = "x3";
		X3.m_ID = oq.AddDiscreteVariable(X3.m_name, 0.5, 6, 0.1);
		setup.m_decVars.Add(X3);
		printf("Added variable x3\n");

		// a discrete variable for sensitivity factor
		DecisionVariable X4;
		X4.m_name = "x4";
		X4.m_ID = oq.AddDiscreteVariable(X4.m_name, 0.8, 1.1, 0.01);
		setup.m_decVars.Add(X4);
		printf("Added variable x4\n");

		// a discrete variable for truck speed acceptance 
		DecisionVariable Y0;
		Y0.m_name = "y0";
		Y0.m_ID = oq.AddDiscreteVariable(Y0.m_name, 0.85, 1.1, 0.01);
		setup.m_decVars.Add(Y0);
		printf("Added variable y0\n");

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
		Y2.m_ID = oq.AddDiscreteVariable(Y2.m_name, 0.6, 1.8, 0.02);
		//Y2.m_ID = oq.AddContinuousVariable(Y2.m_name, 0, 0.3);
		setup.m_decVars.Add(Y2);
		printf("Added variable y2\n");

		// a discrete variable for the truck minDist mean
		DecisionVariable Y3;
		Y3.m_name = "y3";
		Y3.m_ID = oq.AddDiscreteVariable(Y3.m_name, 1, 10, 0.2);
		//Y2.m_ID = oq.AddContinuousVariable(Y2.m_name, 0, 0.3);
		setup.m_decVars.Add(Y3);
		printf("Added variable y3\n");

		// a discrete variable for the truck sensitivity factor
		DecisionVariable Y4;
		Y4.m_name = "y4";
		Y4.m_ID = oq.AddDiscreteVariable(Y4.m_name, 0.85, 1.1, 0.01);
		//Y2.m_ID = oq.AddContinuousVariable(Y2.m_name, 0, 0.3);
		setup.m_decVars.Add(Y4);
		printf("Added variable y4\n");

		// Make sure they are put in the correct order
		// printf("Added variables: (%S, %S, %S, %S, %S, %S, %S \n", setup.m_decVars.GetAt(0).m_name,
		// 	setup.m_decVars.GetAt(1).m_name, setup.m_decVars.GetAt(2).m_name, setup.m_decVars.GetAt(3).m_name,
		//	setup.m_decVars.GetAt(4).m_name, setup.m_decVars.GetAt(5).m_name, setup.m_decVars.GetAt(6).m_name);

		//second argument is setMaximize. If true, then maximize; if false, then minimize
		setup.m_objID = oq.AddUserControlledObjective("myObj", false);
		
		oq.SetMaximumIterations(800);

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

				double x0_sol = oq.GetSolutionVariableValue(solutionID, setup.m_decVars.GetAt(0).m_ID);
				double x1_sol = oq.GetSolutionVariableValue(solutionID, setup.m_decVars.GetAt(1).m_ID);
				double x2_sol = oq.GetSolutionVariableValue(solutionID, setup.m_decVars.GetAt(2).m_ID);
				double x3_sol = oq.GetSolutionVariableValue(solutionID, setup.m_decVars.GetAt(3).m_ID);
				double x4_sol = oq.GetSolutionVariableValue(solutionID, setup.m_decVars.GetAt(4).m_ID);
				
				double y0_sol = oq.GetSolutionVariableValue(solutionID, setup.m_decVars.GetAt(5).m_ID);
				double y1_sol = oq.GetSolutionVariableValue(solutionID, setup.m_decVars.GetAt(6).m_ID);
				double y2_sol = oq.GetSolutionVariableValue(solutionID, setup.m_decVars.GetAt(7).m_ID);
				double y3_sol = oq.GetSolutionVariableValue(solutionID, setup.m_decVars.GetAt(8).m_ID);
				double y4_sol = oq.GetSolutionVariableValue(solutionID, setup.m_decVars.GetAt(9).m_ID);


				if (solfile.is_open())
				{
					// follow the format of paras
					solfile << "car_speedAcceptance," << (x0_sol) << ",0.1,0.85,1.3\n" ;
					solfile << "car_reactionTime," << (x1_sol) << ",1.2,1.6,1\n" ; 
					solfile << "car_maxAccel," << (x2_sol) << ",0.2,2.6,3.4\n" ;
					solfile << "car_minDist," << (x3_sol) << ",0.3,0.5,6\n" ;
					solfile << "car_sensitivityFactor," << (x4_sol) << ",0," << (x4_sol) << "," << (x4_sol) << '\n';

					solfile << "truck_speedAcceptance," << (y0_sol) << ",0.1,0.85,1.1\n" ;
					solfile << "truck_reactionTime," << (y1_sol) << ",1.3,1.7,1\n" ; 
					solfile << "truck_maxAccel," << (y2_sol) << ",0.5,0.6,1.8\n" ;
					solfile << "truck_minDist," << (y3_sol) << ",0.5,1,10\n" ;
					solfile << "truck_sensitivityFactor," << (y4_sol) << ",0," << (y4_sol) << "," << (y4_sol) << '\n';
				}
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


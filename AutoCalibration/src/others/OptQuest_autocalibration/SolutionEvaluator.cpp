// This is the function where 

#include "StdAfx.h"
#include "SolutionEvaluator.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <Windows.h>
// #include <sys/stat.h>


SolutionEvaluator::SolutionEvaluator(void)
{
}

SolutionEvaluator::~SolutionEvaluator(void)
{
}

bool SolutionEvaluator::EvaluateOptQuestSolution(ProblemSetup* setup)
{
	double x0_sol, x1_sol, x2_sol, x3_sol, x4_sol, y0_sol, y1_sol, y2_sol, y3_sol, y4_sol;
	int wait_time;

	std::string para_file ("C:\\sim_com\\sim_paras.txt");
	std::string val_file ("C:\\sim_com\\sim_val.txt");

	try
	{
		// The decision variables in the ProblemSetup have the current solution values.
		// (See COptQuestClient.EvaluateEvent().

		// According to previous assumption, decision variables will be set as: 
		// (x1, x2, x3, x4, x5 y1, y2) for freeflow; 
		// They will be calibrated seperately, meaning, only 4 decision variables here
		
		x0_sol = setup->m_decVars.GetAt(0).m_currentValue;
		x1_sol = setup->m_decVars.GetAt(1).m_currentValue;
		x2_sol = setup->m_decVars.GetAt(2).m_currentValue;
		x3_sol = setup->m_decVars.GetAt(3).m_currentValue;
		x4_sol = setup->m_decVars.GetAt(4).m_currentValue;

		y0_sol = setup->m_decVars.GetAt(5).m_currentValue;
		y1_sol = setup->m_decVars.GetAt(6).m_currentValue;
		y2_sol = setup->m_decVars.GetAt(7).m_currentValue;
		y3_sol = setup->m_decVars.GetAt(8).m_currentValue;
		y4_sol = setup->m_decVars.GetAt(9).m_currentValue;
		

		// Now write the parameters to a file, then wait for the evaluation from AIMSUN
		std::ofstream parasfile;
		parasfile.open(para_file.c_str());
		// parasfile.open("C:/paras.txt");

		if (parasfile.is_open()){
			// follow the format of paras
			// follow the format of paras
			// follow the format of paras
			parasfile << "car_speedAcceptance," << (x0_sol) << ",0.1,0.85,1.3\n" ;
			parasfile << "car_reactionTime," << (x1_sol) << ",1.2,1.6,1\n" ; 
			parasfile << "car_maxAccel," << (x2_sol) << ",0.2,2.6,3.4\n" ;
			parasfile << "car_minDist," << (x3_sol) << ",0.3,0.5,6\n" ;
			parasfile << "car_sensitivityFactor," << (x4_sol) << ",0," << (x4_sol) << "," << (x4_sol) << '\n';

			parasfile << "truck_speedAcceptance," << (y0_sol) << ",0.1,0.85,1.1\n" ;
			parasfile << "truck_reactionTime," << (y1_sol) << ",1.3,1.7,1\n" ; 
			parasfile << "truck_maxAccel," << (y2_sol) << ",0.5,0.6,1.8\n" ;
			parasfile << "truck_minDist," << (y3_sol) << ",0.5,1,10\n" ;
			parasfile << "truck_sensitivityFactor," << (y4_sol) << ",0," << (y4_sol) << "," << (y4_sol) << '\n';
		}
		else 
			printf("unable to open C:\\sim_com\\sim_paras.txt file\n");
			// {printf("unable to open C:/paras.txt \n");}
		parasfile.close();
		printf("Wrote new parameters (%f, %f, %f, %f, %f), (%f, %f, %f, %f, %f)\n", x0_sol, x1_sol, x2_sol, x3_sol, x4_sol, y0_sol, y1_sol, y2_sol, y3_sol, y4_sol);

		// Now wait and read the evaluation result from another file which is to be created and set by Python
		std::ifstream simvalfile;
		simvalfile.open(val_file.c_str());
		// simvalfile.open("C:/simval.txt");

		wait_time = 0;
		printf("waiting for simulation result...\n");
		while (!simvalfile.is_open())
		{
			//wait and open again
			Sleep(1000);	// sleep 1 second
			wait_time += 1;
			if (wait_time >= 600)	// simulation takes a long time; update every 600 s
			{
				printf("waiting for simulation result...\n");
				wait_time = 0;
			}

			simvalfile.open(val_file.c_str());
			// simvalfile.open("C:/simval.txt");
		}
		float simval;
		simvalfile >> simval;
		simvalfile.close();
		
		printf("Got simulation result: %f \n", simval);
		// remove the file, this will be created again when python script gets the simulation result
		if( remove(val_file.c_str()) != 0)
		    perror("Error deleting file C:\\sim_com\\sim_val.txt file");
		//if( remove("C:/simval.txt") != 0)
		//	perror("Error deleting file C:/simval.txt");

		// Evaluate the objective function
		setup->m_objValue = simval;
		// setup->m_objValue = 3*pow(1-x1,2)*exp(-pow(x1,2) - pow(x2+1,2)) - 10*(x1/5 - pow(x1,3) - pow(x2,5))*exp(-pow(x1,2)-pow(x2,2)) - 1/3*exp(-pow(x1+1,2) - pow(x2,2));


		// print out the evaluation result
		// printf("(x1, x2): %f,%f; obj_val: %f\n", x1, x2,setup->m_objValue);

		return true;
	}
	catch (...)
	{
		return false;
	}
	
}

/*
bool fileExists(const std::string& filename)
{
    struct stat buf;
    if (stat(filename.c_str(), &buf) != -1)
    {
        return true;
    }
    return false;
}
*/
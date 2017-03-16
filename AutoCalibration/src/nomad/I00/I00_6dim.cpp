#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h> 
#include <Windows.h>
using namespace std;


int main ( int argc , char ** argv ) {

  double f = 1e20;
  int wait_time;

  // The decision variables for cars
  double x[6];

  // the parameter file and the obj value file
  string para_file ("C:\\sim_com\\sim_paras.txt");
  string val_file ("C:\\sim_com\\sim_val.txt");

  if ( argc >= 2 ) {

    ifstream in ( argv[1] );

    // ===========================================================
    // Fetch the parameters.
    // set the value for the variables
    for ( int i = 0 ; i < 6 ; i++ ) {
      in >> x[i];
    }
    
    // ===========================================================
    // The blackbox AIMSUN simulator
    // -----------------------------------------
    // Write parameters to files which will trigger AIMSUN to simulate
    ofstream parasfile;
    parasfile.open(para_file.c_str());
    if (parasfile.is_open()){
      // Follow the format of setting the parameters
      // mean, std, min, max
      parasfile << "car_speedAcceptance," << (x[0]/100) << ",0.1,0.85,1.1\n" ;
      parasfile << "truck_maxAccel," << (x[1]/100) << ",0.5,0.6,1.8\n" ;
      parasfile << "car_sensitivityFactor," << (x[2]/100) << ",0," << (x[2]/100) << "," << (x[2]/100) << '\n';
      parasfile << "truck_sensitivityFactor," << (x[3]/100) << ",0," << (x[3]/100) << "," << (x[3]/100) << '\n';
      parasfile << "car_reactionTime," << (x[4]*0.8) << ",1.2,1.6,1\n" ;
      parasfile << "truck_reactionTime," << (x[5]*0.8) << ",1.3,1.7,1\n" ;
       
    }
    else
      perror("unable to open C:\\sim_com\\sim_paras.txt file\n");

    parasfile.close();
    // printf("Wrote new parameters (%f, %f)\n", x[0]/100, x[1]/100);

    // -----------------------------------------
    // Now wait and read simulated value from AIMSUN
    ifstream simvalfile;
    simvalfile.open(val_file.c_str());

    wait_time = 0;
    // printf("waiting for simulation result...\n");
    while (!simvalfile.is_open())
    {
      //wait and open again
      Sleep(1000);  // sleep 1 second
      wait_time += 1;
      if (wait_time >= 600) // simulation takes a long time; update every 600 s
      {
        //printf("waiting for simulation result...\n");
        wait_time = 0;
      }

      simvalfile.open(val_file.c_str());
    }
    simvalfile >> f;
    simvalfile.close();

    if( remove(val_file.c_str()) != 0)
        perror("Error deleting file C:\\sim_com\\sim_val.txt file");

    //printf("Got simulation result: %f \n", f);
    // remove the file, this will be created again when python script gets the simulation result

    if ( in.fail() ){
      f = 1e20;}

    in.close();

  }

  cout << f << endl;

  return 0;
}



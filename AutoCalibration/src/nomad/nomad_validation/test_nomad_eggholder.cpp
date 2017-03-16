#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h> 
#include <Windows.h>
using namespace std;

double compute_obj_eggholder(double x, double y){

  double val;

  val = -(y+47)*sin( sqrt( abs( x/2 + (y+47) ) ) ) - x*sin( sqrt( abs( x-(y+47) ) ) );

  return val;
}

double compute_obj_peaks(double x, double y){

  double val;

  val = 3*pow(1-x, 2)*exp(-pow(x, 2) - pow(y+1, 2)) 
    - 10*(x/5 - pow(x, 3) - pow(y, 5))*exp(-pow(x,2)-pow(y, 2))
         - (1/3)*exp(-pow(x+1, 2) - pow(y, 2));

  return val;
}

int main ( int argc , char ** argv ) {

  double f = 1e20;
  double x[2];

  if ( argc >= 2 ) {

    ifstream in ( argv[1] );

    for ( int i = 0 ; i < 2 ; i++ ) {
      in >> x[i];
    }

    f = compute_obj_eggholder( double(x[0])/100, double(x[1])/100 );

    if ( in.fail() )
      f = 1e20;

    in.close();

    // remove the in file
    // if( remove(argv[1]) != 0 )
    //     perror("Can not remove file");

  }

  cout << f << endl;

  return 0;
}



#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
using namespace std;

int main ( int argc , char ** argv ) {

  double z = 1e20;

  if ( argc >= 2 ) {

    ifstream in ( argv[1] );

    double x , y;

    in >> x >> y;

    if ( !in.fail() )
      z = -x*sin(sqrt(fabs(x))) - y*sin(sqrt(fabs(y)));

    in.close();
  }

  cout.setf(ios::fixed);
  cout.precision ( 20 );

  cout << z << endl;
  
  return 0;
}


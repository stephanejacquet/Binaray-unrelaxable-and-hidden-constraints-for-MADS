/*------------------------------------------------*/
/*                   ROSENBROCK                   */
/*------------------------------------------------*/
#include <fstream>
#include <iostream>
#include <cmath>
using namespace std;

#define N   2


int main ( int argc , char ** argv ) {

  double x [N];

  int i;

  if ( argc == 2 ) {

    ifstream in ( argv[1] );

    for ( i = 0 ; i < N ; ++i )
      in >> x[i];

    if ( in.fail() ) {
      cout << 1e20 << endl;
      in.close();
      return 1;
    }

    in.close();
  }

  else  {
    cout << 1e20 << endl;
    return 1;
  }

  double f = 0.0;

  f += pow ( -13.0 + x[0] + ((5.0 -x[1])*x[1] -2.0)*x[1] , 2.0 );
  f += pow ( -29.0 + x[0] + ((x[1]+1.0)*x[1] -14.0)*x[1] , 2.0 );

  cout.precision(15);

  cout << f << endl;
  

  return 0;
}

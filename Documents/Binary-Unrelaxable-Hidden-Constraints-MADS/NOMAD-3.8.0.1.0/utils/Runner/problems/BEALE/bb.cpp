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

  f += pow ( 1.5 - x[0]*(1-pow(x[1],1.0)) , 2.0 );
  f += pow ( 2.25 - x[0]*(1-pow(x[1],2.0)) , 2.0 );
  f += pow ( 2.625 - x[0]*(1-pow(x[1],3.0)) , 2.0 );

    
  cout.precision(25);

  cout << f << endl;
  

  return 0;
}

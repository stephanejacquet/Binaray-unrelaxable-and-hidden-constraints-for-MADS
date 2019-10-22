#include <string> 
#include "nomad.hpp"

using namespace std;
using namespace NOMAD;

/*------------------------------------------------*/
/*               The problem                      */
/*------------------------------------------------*/
/*       n=4, m=3                                 */
/*       PressureVessel Mixed Case                       */
/*------------------------------------------------*/
class My_Evaluator : public Evaluator {

public:

  // ctor:
  My_Evaluator  ( const Parameters & p ) :
    Evaluator ( p ) {}

  // dtor:
  ~My_Evaluator ( void ) {}

  // evaluation of a point:
  bool eval_x ( Eval_Point          & x          ,
		const NOMAD::Double & h_max      ,
		bool                & count_eval   ) const 
  {

    double f, g1, g2, g3,x1,x2, x3, x4;
    const double g_Pi = 3.14159265358979323846;

    x1=x[0].value();
    x2=x[1].value();
    x3=x[2].value();
    x4=x[3].value();

// Cost function
    f = 0.6224 * 0.0625 * x1 * x3 * x4 + 1.7781 * 0.0625 * x2 * x3 * x3 + 3.1661 * 0.0625 * 0.0625 * x1 * x1 * x4 + 19.84 * 0.0625 * 0.0625 * x1 * x1 * x3;

    // Constraints
    g1 = - 0.0625 * x1 + 0.0193 * x3;   // <= 0
    g2 = - 0.0625 * x2 + 0.00954 * x3;  // <= 0
    g3 = - g_Pi * x3 * x3 * x4 - 4.0 / 3.0 * g_Pi * x3 * x3 * x3 + 1296000.0;  // <= 0

    count_eval = true; // count a black-box evaluation

    x.set_bb_output  ( 0 , f  ); // objective value
    x.set_bb_output  ( 1 , g1  ); 
    x.set_bb_output  ( 2 , g2  ); 
    x.set_bb_output  ( 3 , g3  ); 

    return true;       // the evaluation succeeded
  }
};

/*------------------------------------------*/
/*            NOMAD main function           */
/*------------------------------------------*/
int main ( int argc , char ** argv ) {

  // NOMAD initializations:
  begin ( argc , argv );

  // display:
  Display out ( std::cout );
  out.precision ( DISPLAY_PRECISION_STD );

  //  parameters creation:
  Parameters p ( out );

  int dim = 4;
  int cont = 3;
  int nbpoints = 1;
  int nbiter = 100;

  p.set_DIMENSION (dim);             // number of variables

  vector<bb_output_type> bbot (cont+1); // definition of
  bbot[0] = OBJ;                   // output types
  bbot[1] = PB;                   
  bbot[2] = PB;                   
  bbot[3] = PB;                   
  p.set_BB_OUTPUT_TYPE ( bbot );
  
  Point x0(dim); // starting points
  string temp;
  int temp2;
  ifstream fich("DOE.plan");
  for (int k=0; k<7; k++)
    {
      getline(fich,temp);
    }
  for (int k=0; k<nbpoints; k++)
    {
      fich >> temp2 >> temp;
      fich >> x0[0]; 
      fich >> x0[1];
      fich >> x0[2]; 
      fich >> x0[3];
      getline(fich,temp);
      p.set_X0 ( x0 );
    }
  fich.close();

  x0[0]=1;
  x0[1]=1;
  x0[2]=10.0;
  x0[3]=10.0;
  p.set_LOWER_BOUND ( x0 );
  x0[0]=99;
  x0[1]=99;
  x0[2]=200.0;
  x0[3]=200.0;
  p.set_UPPER_BOUND ( x0 );

  p.set_BB_INPUT_TYPE (0 , INTEGER);
  p.set_BB_INPUT_TYPE (1 , INTEGER);
    
    
    p.set_SEED(3183);

  // p.set_MAX_BB_EVAL ( nbpoints+nbiter );

  // p.set_DIRECTION_TYPE(NOMAD::ORTHO_2N);

//  p.set_OPPORTUNISTIC_EVAL(false);

  //p.set_SPECULATIVE_SEARCH ( false );

  //p.set_MODEL_SEARCH( NO_MODEL );

  // p.set_DISPLAY_DEGREE ( 3 );
  //p.set_DISPLAY_DEGREE ( "0300" );// display only the search step

  // p.set_DISPLAY_ALL_EVAL ( true );

  // p.set_DISPLAY_STATS ( "bbe sol obj" );

//  p.set_STATS_FILE ( "stat.txt","BBE SOL %=.10eOBJ %.10eBBO");

 // p.set_ADD_SEED_TO_FILE_NAMES ( false );

  //p.set_SOLUTION_FILE ("sol.txt");

 // p.set_HISTORY_FILE ("hist.txt");

  // parameters validation:
  p.check();

  // custom evaluator creation:
  My_Evaluator ev ( p );

  // algorithm creation:
  Mads mads ( p , &ev );

  // algorithm execution:
  mads.run();

  Slave::stop_slaves ( out );
  end();

  return EXIT_SUCCESS;
}
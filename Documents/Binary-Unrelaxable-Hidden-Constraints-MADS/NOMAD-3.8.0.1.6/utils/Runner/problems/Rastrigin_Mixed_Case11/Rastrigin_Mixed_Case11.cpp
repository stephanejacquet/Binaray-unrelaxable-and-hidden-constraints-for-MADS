#include <string>
#include "Rastrigin_Mixed_Case11.hpp"

using namespace std;
using namespace NOMAD;

/*------------------------------------------------*/
/*               The problem                      */
/*------------------------------------------------*/
/*       n=2, m=1                                 */
/*       Rastrigin Mixed Case 1                     */
/*------------------------------------------------*/
/*----------------------------------------------*/
/*                  constructor                 */
/*----------------------------------------------*/
Rastrigin_Mixed_Case11::Rastrigin_Mixed_Case11 ( void )
: Problem ( "Rastrigin_Mixed_Case11" , "Rastrigin_Mixed_Case11" , "xe.txt" , 2 , 1 ) {
    
    
    set_bbot ( 0 , NOMAD::OBJ );

    int dim = 2;
    
    // TAKE ONLY THE FIRST POINT FROM DOE
    int nbpoints = 1;
    
    Point x0(dim); // starting points
    string temp;
    int temp2;
    ifstream fich("problems/Rastrigin_Mixed_Case11/DOE.plan");
    for (int k=0; k<7; k++)
    {
        getline(fich,temp);
    }
    for (int k=0; k<nbpoints; k++)
    {
        fich >> temp2 >> temp;
        fich >> x0[0];
        if (x0[0].value()==-5)
            x0[0]=0;
        else if (x0[0].value()==-3)
            x0[0]=1;
        else if (x0[0].value()==-1)
            x0[0]=2;
        else if (x0[0].value()==0)
            x0[0]=3;
        else if (x0[0].value()==1)
            x0[0]=4;
        else if (x0[0].value()==3)
            x0[0]=5;
        else if (x0[0].value()==5)
            x0[0]=6;
        fich >> x0[1];
        getline(fich,temp);
        set_x0 ( x0 );
    }
    fich.close();
    
    Point xl(dim);
    xl[0]=0;
    xl[1]=-5.0;
    
    Point xu(dim);
    xu[0]=6;
    xu[1]=5.0;
    
    set_bounds ( xl , xu );
    
    set_bbit (0 , INTEGER );
    
    add_keyword ( "published"    );
    add_keyword ( "constrained" );
    add_keyword ( "anne_sophie_crelot");
}


bool Rastrigin_Mixed_Case11::eval_x ( Eval_Point          & x          ,
                                  bool                & count_eval   ) const
{
    
    double f,x1,x2;
    double pi = 3.1415926535897932385;
    int n=2;
    
    switch (static_cast<int>(x[0].value())) {
        case 0:
            x1=-5;
            break;
        case 1:
            x1=-3;
            break;
        case 2:
            x1=-1;
            break;
        case 3:
            x1=0;
            break;
        case 4:
            x1=1;
            break;
        case 5:
            x1=3;
            break;
        case 6:
            x1=5;
            break;
    }
    x2=x[1].value();
    
    f = 10.0 * n + x1*x1-10.0*cos(2.0*pi*x1) + x2*x2-10.0*cos(2.0*pi*x2);
    
    count_eval = true; // count a black-box evaluation
    
    x.set_bb_output  ( 0 , f  ); // objective value
    
    return true;       // the evaluation succeeded
}


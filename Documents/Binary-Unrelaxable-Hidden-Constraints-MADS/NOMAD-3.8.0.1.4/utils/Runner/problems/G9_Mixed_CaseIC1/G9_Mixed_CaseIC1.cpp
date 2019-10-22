#include <string>
#include "G9_Mixed_CaseIC1.hpp"

using namespace std;
using namespace NOMAD;

/*------------------------------------------------*/
/*               The problem                      */
/*------------------------------------------------*/
/*       n=7, m=5                                 */
/*       G9 Mixed Case                            */
/*------------------------------------------------*/
/*----------------------------------------------*/
/*                  constructor                 */
/*----------------------------------------------*/
G9_Mixed_CaseIC1::G9_Mixed_CaseIC1 ( void )
: Problem ( "G9_Mixed_CaseIC1" , "G9_Mixed_CaseIC1" , "xe.txt" , 7 , 5 ) {
    
    
    set_bbot ( 0 , NOMAD::OBJ );
    set_bbot ( 1 , NOMAD::PB );
    set_bbot ( 2 , NOMAD::PB );
    set_bbot ( 3 , NOMAD::PB );
    set_bbot ( 4 , NOMAD::PB );



    int dim = 7;
    
    // TAKE ONLY THE FIRST POINT FROM DOE
    int nbpoints = 1;
    
    Point x0(dim); // starting points
    string temp;
    int temp2;
    ifstream fich("problems/G9_Mixed_CaseIC1/DOE.plan");
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
        fich >> x0[4];
        fich >> x0[5];
        fich >> x0[6];
        getline(fich,temp);
        set_x0 ( x0 );
    }
    fich.close();
    
    Point xl(dim);
    xl[0]=-10.0;
    xl[1]=-10.0;
    xl[2]=-10.0;
    xl[3]=-10.0;
    xl[4]=-10.0;
    xl[5]=-10.0;
    xl[6]=-10.0;
    
    Point xu(dim);
    xu[0]=10.0;
    xu[1]=10.0;
    xu[2]=10.0;
    xu[3]=10.0;
    xu[4]=10.0;
    xu[5]=10.0;
    xu[6]=10.0;
    
    set_bounds ( xl , xu );
    
    set_bbit (0 , CATEGORICAL );
    set_bbit (1 , CATEGORICAL );
    set_bbit (2 , CATEGORICAL );
    
    
    
    add_keyword ( "published"    );
    add_keyword ( "constrained" );
    add_keyword ( "anne_sophie_crelot");
    add_keyword ( "categorical" );
}


bool G9_Mixed_CaseIC1::eval_x ( Eval_Point          & x          ,
                                  bool                & count_eval   ) const
{
    Double f;
    Double g1, g2, g3, g4;
    
    /* objective function */
    f = (x[0] - 10.0) * (x[0] - 10.0) + 5.0 * (x[1] - 12.0) * (x[1] - 12.0) + std::pow (x[2].value(), 4) + 3.0 * (x[3] - 11.0) * (x[3] - 11.0) + 10.0 * std::pow (x[4].value(), 6) + 7.0 * x[5] * x[5] + std::pow (x[6].value(), 4) - 4.0 * x[5] * x[6] - 10.0 * x[5] - 8.0 * x[6];
    /* constraints g<=0 */
    g1 = -127.0 + 2 * x[0] * x[0] + 3.0 * std::pow (x[1].value(), 4) + x[2] + 4.0 * x[3] * x[3] + 5.0 * x[4];
    g2 = -282.0 + 7.0 * x[0] + 3.0 * x[1] + 10.0 * x[2] * x[2] + x[3] - x[4];
    g3 = -196.0 + 23.0 * x[0] + x[1] * x[1] + 6.0 * x[5] * x[5] - 8.0 * x[6];
    g4 = 4.0 * x[0] * x[0] + x[1] * x[1] - 3.0 * x[0] * x[1] + 2.0 * x[2] * x[2] + 5.0 * x[5] - 11.0 * x[6];
    
    count_eval = true; // count a black-box evaluation
    
    x.set_bb_output  ( 0 , f  ); // objective value
    x.set_bb_output  ( 1 , g1  );
    x.set_bb_output  ( 2 , g2  );
    x.set_bb_output  ( 3 , g3  );
    x.set_bb_output  ( 4 , g4  );
    
    return true;       // the evaluation succeeded
}


/*--------------------------------------*/
/*  construct the extended poll points  */
/*      (categorical neighborhoods)     */
/*--------------------------------------*/
void EP_G9_Mixed_CaseIC1::construct_extended_points ( const Eval_Point & x )
{
    for (int k=0; k<3; k++)
    {
        vector<int> other_types;
        if (x[k].value()==-10)
        {
            other_types.push_back(-9);
            other_types.push_back(-8);
            other_types.push_back(-7);
            other_types.push_back(-6);
            other_types.push_back(-5);
            other_types.push_back(-4);
            other_types.push_back(-3);
            other_types.push_back(-2);
            other_types.push_back(-1);
            other_types.push_back(0);
            other_types.push_back(1);
            other_types.push_back(2);
            other_types.push_back(3);
            other_types.push_back(4);
            other_types.push_back(5);
            other_types.push_back(6);
            other_types.push_back(7);
            other_types.push_back(8);
            other_types.push_back(9);
            other_types.push_back(10);
        }
        else if (x[k].value()==-9)
        {
            other_types.push_back(-10);
            other_types.push_back(-8);
            other_types.push_back(-7);
            other_types.push_back(-6);
            other_types.push_back(-5);
            other_types.push_back(-4);
            other_types.push_back(-3);
            other_types.push_back(-2);
            other_types.push_back(-1);
            other_types.push_back(0);
            other_types.push_back(1);
            other_types.push_back(2);
            other_types.push_back(3);
            other_types.push_back(4);
            other_types.push_back(5);
            other_types.push_back(6);
            other_types.push_back(7);
            other_types.push_back(8);
            other_types.push_back(9);
            other_types.push_back(10);
        }
        else if (x[k].value()==-8)
        {
            other_types.push_back(-10);
            other_types.push_back(-9);
            other_types.push_back(-7);
            other_types.push_back(-6);
            other_types.push_back(-5);
            other_types.push_back(-4);
            other_types.push_back(-3);
            other_types.push_back(-2);
            other_types.push_back(-1);
            other_types.push_back(0);
            other_types.push_back(1);
            other_types.push_back(2);
            other_types.push_back(3);
            other_types.push_back(4);
            other_types.push_back(5);
            other_types.push_back(6);
            other_types.push_back(7);
            other_types.push_back(8);
            other_types.push_back(9);
            other_types.push_back(10);
        }
        else if (x[k].value()==-7)
        {
            other_types.push_back(-10);
            other_types.push_back(-9);
            other_types.push_back(-8);
            other_types.push_back(-6);
            other_types.push_back(-5);
            other_types.push_back(-4);
            other_types.push_back(-3);
            other_types.push_back(-2);
            other_types.push_back(-1);
            other_types.push_back(0);
            other_types.push_back(1);
            other_types.push_back(2);
            other_types.push_back(3);
            other_types.push_back(4);
            other_types.push_back(5);
            other_types.push_back(6);
            other_types.push_back(7);
            other_types.push_back(8);
            other_types.push_back(9);
            other_types.push_back(10);
        }
        else if (x[k].value()==-6)
        {
            other_types.push_back(-10);
            other_types.push_back(-9);
            other_types.push_back(-8);
            other_types.push_back(-7);
            other_types.push_back(-5);
            other_types.push_back(-4);
            other_types.push_back(-3);
            other_types.push_back(-2);
            other_types.push_back(-1);
            other_types.push_back(0);
            other_types.push_back(1);
            other_types.push_back(2);
            other_types.push_back(3);
            other_types.push_back(4);
            other_types.push_back(5);
            other_types.push_back(6);
            other_types.push_back(7);
            other_types.push_back(8);
            other_types.push_back(9);
            other_types.push_back(10);
        }
        else if (x[k].value()==-5)
        {
            other_types.push_back(-10);
            other_types.push_back(-9);
            other_types.push_back(-8);
            other_types.push_back(-7);
            other_types.push_back(-6);
            other_types.push_back(-4);
            other_types.push_back(-3);
            other_types.push_back(-2);
            other_types.push_back(-1);
            other_types.push_back(0);
            other_types.push_back(1);
            other_types.push_back(2);
            other_types.push_back(3);
            other_types.push_back(4);
            other_types.push_back(5);
            other_types.push_back(6);
            other_types.push_back(7);
            other_types.push_back(8);
            other_types.push_back(9);
            other_types.push_back(10);
        }
        else if (x[k].value()==-4)
        {
            other_types.push_back(-10);
            other_types.push_back(-9);
            other_types.push_back(-8);
            other_types.push_back(-7);
            other_types.push_back(-6);
            other_types.push_back(-5);
            other_types.push_back(-3);
            other_types.push_back(-2);
            other_types.push_back(-1);
            other_types.push_back(0);
            other_types.push_back(1);
            other_types.push_back(2);
            other_types.push_back(3);
            other_types.push_back(4);
            other_types.push_back(5);
            other_types.push_back(6);
            other_types.push_back(7);
            other_types.push_back(8);
            other_types.push_back(9);
            other_types.push_back(10);
        }
        else if (x[k].value()==-3)
        {
            other_types.push_back(-10);
            other_types.push_back(-9);
            other_types.push_back(-8);
            other_types.push_back(-7);
            other_types.push_back(-6);
            other_types.push_back(-5);
            other_types.push_back(-4);
            other_types.push_back(-2);
            other_types.push_back(-1);
            other_types.push_back(0);
            other_types.push_back(1);
            other_types.push_back(2);
            other_types.push_back(3);
            other_types.push_back(4);
            other_types.push_back(5);
            other_types.push_back(6);
            other_types.push_back(7);
            other_types.push_back(8);
            other_types.push_back(9);
            other_types.push_back(10);
        }
        else if (x[k].value()==-2)
        {
            other_types.push_back(-10);
            other_types.push_back(-9);
            other_types.push_back(-8);
            other_types.push_back(-7);
            other_types.push_back(-6);
            other_types.push_back(-5);
            other_types.push_back(-4);
            other_types.push_back(-3);
            other_types.push_back(-1);
            other_types.push_back(0);
            other_types.push_back(1);
            other_types.push_back(2);
            other_types.push_back(3);
            other_types.push_back(4);
            other_types.push_back(5);
            other_types.push_back(6);
            other_types.push_back(7);
            other_types.push_back(8);
            other_types.push_back(9);
            other_types.push_back(10);
        }
        else if (x[k].value()==-1)
        {
            other_types.push_back(-10);
            other_types.push_back(-9);
            other_types.push_back(-8);
            other_types.push_back(-7);
            other_types.push_back(-6);
            other_types.push_back(-5);
            other_types.push_back(-4);
            other_types.push_back(-3);
            other_types.push_back(-2);
            other_types.push_back(0);
            other_types.push_back(1);
            other_types.push_back(2);
            other_types.push_back(3);
            other_types.push_back(4);
            other_types.push_back(5);
            other_types.push_back(6);
            other_types.push_back(7);
            other_types.push_back(8);
            other_types.push_back(9);
            other_types.push_back(10);
        }
        else if (x[k].value()==0)
        {
            other_types.push_back(-10);
            other_types.push_back(-9);
            other_types.push_back(-8);
            other_types.push_back(-7);
            other_types.push_back(-6);
            other_types.push_back(-5);
            other_types.push_back(-4);
            other_types.push_back(-3);
            other_types.push_back(-2);
            other_types.push_back(-1);
            other_types.push_back(1);
            other_types.push_back(2);
            other_types.push_back(3);
            other_types.push_back(4);
            other_types.push_back(5);
            other_types.push_back(6);
            other_types.push_back(7);
            other_types.push_back(8);
            other_types.push_back(9);
            other_types.push_back(10);
        }
        else if (x[k].value()==1)
        {
            other_types.push_back(-10);
            other_types.push_back(-9);
            other_types.push_back(-8);
            other_types.push_back(-7);
            other_types.push_back(-6);
            other_types.push_back(-5);
            other_types.push_back(-4);
            other_types.push_back(-3);
            other_types.push_back(-2);
            other_types.push_back(-1);
            other_types.push_back(0);
            other_types.push_back(2);
            other_types.push_back(3);
            other_types.push_back(4);
            other_types.push_back(5);
            other_types.push_back(6);
            other_types.push_back(7);
            other_types.push_back(8);
            other_types.push_back(9);
            other_types.push_back(10);
        }
        else if (x[k].value()==2)
        {
            other_types.push_back(-10);
            other_types.push_back(-9);
            other_types.push_back(-8);
            other_types.push_back(-7);
            other_types.push_back(-6);
            other_types.push_back(-5);
            other_types.push_back(-4);
            other_types.push_back(-3);
            other_types.push_back(-2);
            other_types.push_back(-1);
            other_types.push_back(0);
            other_types.push_back(1);
            other_types.push_back(3);
            other_types.push_back(4);
            other_types.push_back(5);
            other_types.push_back(6);
            other_types.push_back(7);
            other_types.push_back(8);
            other_types.push_back(9);
            other_types.push_back(10);
        }
        else if (x[k].value()==3)
        {
            other_types.push_back(-10);
            other_types.push_back(-9);
            other_types.push_back(-8);
            other_types.push_back(-7);
            other_types.push_back(-6);
            other_types.push_back(-5);
            other_types.push_back(-4);
            other_types.push_back(-3);
            other_types.push_back(-2);
            other_types.push_back(-1);
            other_types.push_back(0);
            other_types.push_back(1);
            other_types.push_back(2);
            other_types.push_back(4);
            other_types.push_back(5);
            other_types.push_back(6);
            other_types.push_back(7);
            other_types.push_back(8);
            other_types.push_back(9);
            other_types.push_back(10);
        }
        else if (x[k].value()==4)
        {
            other_types.push_back(-10);
            other_types.push_back(-9);
            other_types.push_back(-8);
            other_types.push_back(-7);
            other_types.push_back(-6);
            other_types.push_back(-5);
            other_types.push_back(-4);
            other_types.push_back(-3);
            other_types.push_back(-2);
            other_types.push_back(-1);
            other_types.push_back(0);
            other_types.push_back(1);
            other_types.push_back(2);
            other_types.push_back(3);
            other_types.push_back(5);
            other_types.push_back(6);
            other_types.push_back(7);
            other_types.push_back(8);
            other_types.push_back(9);
            other_types.push_back(10);
        }
        else if (x[k].value()==5)
        {
            other_types.push_back(-10);
            other_types.push_back(-9);
            other_types.push_back(-8);
            other_types.push_back(-7);
            other_types.push_back(-6);
            other_types.push_back(-5);
            other_types.push_back(-4);
            other_types.push_back(-3);
            other_types.push_back(-2);
            other_types.push_back(-1);
            other_types.push_back(0);
            other_types.push_back(1);
            other_types.push_back(2);
            other_types.push_back(3);
            other_types.push_back(4);
            other_types.push_back(6);
            other_types.push_back(7);
            other_types.push_back(8);
            other_types.push_back(9);
            other_types.push_back(10);
        }
        else if (x[k].value()==6)
        {
            other_types.push_back(-10);
            other_types.push_back(-9);
            other_types.push_back(-8);
            other_types.push_back(-7);
            other_types.push_back(-6);
            other_types.push_back(-5);
            other_types.push_back(-4);
            other_types.push_back(-3);
            other_types.push_back(-2);
            other_types.push_back(-1);
            other_types.push_back(0);
            other_types.push_back(1);
            other_types.push_back(2);
            other_types.push_back(3);
            other_types.push_back(4);
            other_types.push_back(5);
            other_types.push_back(7);
            other_types.push_back(8);
            other_types.push_back(9);
            other_types.push_back(10);
        }
        else if (x[k].value()==7)
        {
            other_types.push_back(-10);
            other_types.push_back(-9);
            other_types.push_back(-8);
            other_types.push_back(-7);
            other_types.push_back(-6);
            other_types.push_back(-5);
            other_types.push_back(-4);
            other_types.push_back(-3);
            other_types.push_back(-2);
            other_types.push_back(-1);
            other_types.push_back(0);
            other_types.push_back(1);
            other_types.push_back(2);
            other_types.push_back(3);
            other_types.push_back(4);
            other_types.push_back(5);
            other_types.push_back(6);
            other_types.push_back(8);
            other_types.push_back(9);
            other_types.push_back(10);
        }
        else if (x[k].value()==8)
        {
            other_types.push_back(-10);
            other_types.push_back(-9);
            other_types.push_back(-8);
            other_types.push_back(-7);
            other_types.push_back(-6);
            other_types.push_back(-5);
            other_types.push_back(-4);
            other_types.push_back(-3);
            other_types.push_back(-2);
            other_types.push_back(-1);
            other_types.push_back(0);
            other_types.push_back(1);
            other_types.push_back(2);
            other_types.push_back(3);
            other_types.push_back(4);
            other_types.push_back(5);
            other_types.push_back(6);
            other_types.push_back(7);
            other_types.push_back(9);
            other_types.push_back(10);
        }
        else if (x[k].value()==9)
        {
            other_types.push_back(-10);
            other_types.push_back(-9);
            other_types.push_back(-8);
            other_types.push_back(-7);
            other_types.push_back(-6);
            other_types.push_back(-5);
            other_types.push_back(-4);
            other_types.push_back(-3);
            other_types.push_back(-2);
            other_types.push_back(-1);
            other_types.push_back(0);
            other_types.push_back(1);
            other_types.push_back(2);
            other_types.push_back(3);
            other_types.push_back(4);
            other_types.push_back(5);
            other_types.push_back(6);
            other_types.push_back(7);
            other_types.push_back(8);
            other_types.push_back(10);
        }
        else
        {
            other_types.push_back(-10);
            other_types.push_back(-9);
            other_types.push_back(-8);
            other_types.push_back(-7);
            other_types.push_back(-6);
            other_types.push_back(-5);
            other_types.push_back(-4);
            other_types.push_back(-3);
            other_types.push_back(-2);
            other_types.push_back(-1);
            other_types.push_back(0);
            other_types.push_back(1);
            other_types.push_back(2);
            other_types.push_back(3);
            other_types.push_back(4);
            other_types.push_back(5);
            other_types.push_back(6);
            other_types.push_back(7);
            other_types.push_back(8);
            other_types.push_back(9);
        }
        
        for ( size_t j = 0 ; j < other_types.size() ; j++ )
        {
            Point y = x ;
            y[k] = other_types[j];
            add_extended_poll_point ( y , *(x.get_signature())  );
        }
    }
}



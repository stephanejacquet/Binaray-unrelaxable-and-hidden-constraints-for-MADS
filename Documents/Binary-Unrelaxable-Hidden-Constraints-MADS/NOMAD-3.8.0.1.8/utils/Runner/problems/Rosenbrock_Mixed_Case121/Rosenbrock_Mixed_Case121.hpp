#ifndef __ROSENBROCKMIXEDCASE121__
#define __ROSENBROCKMIXEDCASE121__

#include "../../Problem.hpp"

class Rosenbrock_Mixed_Case121: public Problem {

public:

  Rosenbrock_Mixed_Case121 ( void );

  virtual ~Rosenbrock_Mixed_Case121 ( void ) {}

  virtual bool eval_x ( NOMAD::Eval_Point & x          ,
 			bool              & count_eval   ) const;
};

/*--------------------------------------------------*/
/*  user class to define categorical neighborhoods  */
/*--------------------------------------------------*/
class EP_Rosenbrock_Mixed_Case121 : public NOMAD::Extended_Poll
{
public:
    
    // constructor:
    EP_Rosenbrock_Mixed_Case121 ( NOMAD::Parameters & p):
    NOMAD::Extended_Poll ( p    ) {}
    
    // destructor:
    virtual ~EP_Rosenbrock_Mixed_Case121 ( void ) {}
    
    // construct the extended poll points:
    virtual void construct_extended_points ( const NOMAD::Eval_Point & );
    
};

#endif

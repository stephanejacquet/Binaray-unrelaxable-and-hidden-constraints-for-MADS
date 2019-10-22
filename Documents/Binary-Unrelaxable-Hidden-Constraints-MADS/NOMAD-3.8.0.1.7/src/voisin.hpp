/* This is a modified version from Nomad 3.8.0. All changes were made by Stéphane Jacquet. The original version of Nomad 3.8.0 can be downloaded on https://www.gerad.ca/nomad/ */
/* 
 * File:   voisin.hpp
 * Author: jacqstep
 *
 * Created on 21 février 2017, 14:26
 */

#ifndef VOISIN_HPP
#define	VOISIN_HPP
#include "Double.hpp"
class Voisin
{
  public:
    int			origin;
    
    NOMAD::Double 	distance; 	//right side of violated Inequality
		
    
    Voisin(){};
    
    //bool operator< (const Voisin&) const;
  
};



bool Voisin::operator<(const Voisin& c1) const{
  return (distance <= c1.distance);	
}

#endif	/* VOISIN_HPP */


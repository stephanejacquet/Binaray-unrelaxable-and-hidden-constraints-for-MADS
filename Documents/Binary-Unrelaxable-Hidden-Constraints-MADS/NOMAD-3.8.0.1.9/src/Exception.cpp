/* This is a modified version from Nomad 3.8.0. All changes were made by Stéphane Jacquet. The original version of Nomad 3.8.0 can be downloaded on https://www.gerad.ca/nomad/ */
/*-------------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct search - version 3.8.0      */
/*                                                                                     */
/*                                                                                     */
/*  NOMAD - version 3.8.0 has been created by                                          */
/*                 Charles Audet        - Ecole Polytechnique de Montreal              */
/*                 Sebastien Le Digabel - Ecole Polytechnique de Montreal              */
/*                 Christophe Tribes    - Ecole Polytechnique de Montreal              */
/*                                                                                     */
/*  The copyright of NOMAD - version 3.8.0 is owned by                                 */
/*                 Sebastien Le Digabel - Ecole Polytechnique de Montreal              */
/*                 Christophe Tribes    - Ecole Polytechnique de Montreal              */
/*                                                                                     */
/*  NOMAD v3 has been funded by AFOSR and Exxon Mobil.                                 */
/*                                                                                     */
/*  NOMAD v3 is a new version of NOMAD v1 and v2. NOMAD v1 and v2 were created and     */
/*  developed by Mark Abramson, Charles Audet, Gilles Couture and John E. Dennis Jr.,  */
/*  and were funded by AFOSR and Exxon Mobil.                                          */
/*                                                                                     */
/*                                                                                     */
/*  Contact information:                                                               */
/*    Ecole Polytechnique de Montreal - GERAD                                          */
/*    C.P. 6079, Succ. Centre-ville, Montreal (Quebec) H3C 3A7 Canada                  */
/*    e-mail: nomad@gerad.ca                                                           */
/*    phone : 1-514-340-6053 #6928                                                     */
/*    fax   : 1-514-340-5665                                                           */
/*                                                                                     */
/*  This program is free software: you can redistribute it and/or modify it under the  */
/*  terms of the GNU Lesser General Public License as published by the Free Software   */
/*  Foundation, either version 3 of the License, or (at your option) any later         */
/*  version.                                                                           */
/*                                                                                     */
/*  This program is distributed in the hope that it will be useful, but WITHOUT ANY    */
/*  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A    */
/*  PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.   */
/*                                                                                     */
/*  You should have received a copy of the GNU Lesser General Public License along     */
/*  with this program. If not, see <http://www.gnu.org/licenses/>.                     */
/*                                                                                     */
/*  You can find information on the NOMAD software at www.gerad.ca/nomad               */
/*-------------------------------------------------------------------------------------*/
/**
 \file   Exception.cpp
 \brief  custom class for exceptions (implementation)
 \author Sebastien Le Digabel
 \date   2010-03-29
 \see    Exception.hpp
 */
#include "Exception.hpp"

/*----------------------------------------------------------------*/
/*                     NOMAD::Exception::what()                   */
/*----------------------------------------------------------------*/
const char * NOMAD::Exception::what ( void ) const throw()
{
    std::ostringstream oss;
    if ( ! _file.empty() || _line > 0 )
        oss << "NOMAD::Exception thrown (" << _file << ", " << _line << ")";
    if ( !_what.empty() )
        oss << " " << _what;
    _what = oss.str();
    return _what.c_str();
}

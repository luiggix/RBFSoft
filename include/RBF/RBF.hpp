
/*------------------------------------------------------------------------
 *  Copyright (C) 2008  Luis M. de la Cruz
 *
 *  This file is part of TUNA::RBF
 *
 *  TUNA::RBF is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  TUNA::RBF is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ------------------------------------------------------------------------*/

#ifndef _RBF_HPP_
#define _RBF_HPP_

/*! 
 ***************************************************************************
 *  \namespace  RBF
 *  Implements several RBFs.
 *  This namespace contains implementations to evaluate several Radial Basis 
 *  Functions and its derivatives.
 *  \author  Luis M. de la Cruz [ Mon Nov 12 11:52:45 GMT 2007 ]
 ***************************************************************************
 */

namespace RBF {

/*!
 *  Function to evaluate the solution using an RBF in 3D.
 *  \param x coordinate where the solution will be evaluated.
 *  \param y coordinate where the solution will be evaluated.
 *  \param z coordinate where the solution will be evaluated.
 *  \param xv vector of x-coords where the RBF was calculated.
 *  \param yv vector of y-coords where the RBF was calculated.
 *  \param zv vector of z-coords where the RBF was calculated.
 *  \param lambda vector of coefficients for the RBF.
 *  \param NPol 0 -> not polynomial (default), 6 -> using a degree 2 pol.
 *                   (just for testing, could be eliminated...)
 */
template<typename T_RBF, typename Tprec>
inline
Tprec eval(Tprec x, Tprec y, Tprec z, 
	   const Vec& xv, const Vec& yv, const Vec& zv, 
	   const Vec& lambda, T_RBF f, int NPol = 0)
{
  //    int N = xv.length();    
    int N = xv.size();    
    Tprec interp = 0.0;    
    
    for(int i = 0; i < N; ++i) 
	interp += lambda(i) * f( x, y, z, xv(i), yv(i), zv(i) );
    
    if (NPol) { // Include polinomyal ???	
	interp += lambda(N + 0) * x * x +
	    lambda(N + 1) * y * y +
	    lambda(N + 2) * x * y +
	    lambda(N + 3) * x +
	    lambda(N + 4) * y +
	    lambda(N + 5);  
    }
    
    return interp;  
}


/*!
 *  Function to evaluate the solution using an RBF in 2D.
 *  \param x coordinate where the solution will be evaluated.
 *  \param y coordinate where the solution will be evaluated.
 *  \param xv vector of x-coords where the RBF was calculated.
 *  \param yv vector of y-coords where the RBF was calculated.
 *  \param lambda vector of coefficients for the RBF.
 *  \param NPol 0 -> not polynomial (default), 6 -> using a degree 2 pol.
 *                   (just for testing, could be eliminated...)
 */
template<typename T_RBF, typename Tprec>
inline
Tprec eval(Tprec x, Tprec y, 
	   const Vec& xv, const Vec& yv, const Vec& lambda, T_RBF f, 
	   int NPol = 0)
{
  //    int N = xv.length();    
    int N = xv.size();    
    Tprec interp = 0.0;    
    
    for(int i = 0; i < N; ++i) 
	interp += lambda(i) * f( x, y, xv(i), yv(i) );
    
    if (NPol) { // Include polinomyal ???	
	interp += lambda(N + 0) * x * x +
	    lambda(N + 1) * y * y +
	    lambda(N + 2) * x * y +
	    lambda(N + 3) * x +
	    lambda(N + 4) * y +
	    lambda(N + 5);  
    }
    
    return interp;  
}

/*!
 *  Function to evaluate the solution using an RBF in 1D.
 *  \param x coordinate where the solution will be evaluated.
 *  \param xv vector of x-coords where the RBF was calculated.
 *  \param lambda vector of coefficients for the RBF.
 */
template<typename T_RBF, typename Tprec>
inline
Tprec eval(Tprec x, const Vec& xv, const Vec& lambda, T_RBF f)
{
  //    int N = xv.length();     
    int N = xv.size();     
    Tprec interp = 0.0;    

    for(int i = 0; i < N; ++i)
	interp += lambda(i) * f( x, xv(i) );  
    return interp;  
}

} // namespace RBF

#endif

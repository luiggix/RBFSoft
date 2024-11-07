
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

#ifndef _TPS_HPP_
#define _TPS_HPP_

#include <cmath>

namespace RBF {

/*! 
 ***************************************************************************
 *  Generic functor to evaluate Thin Plate Spline.
 *  All the functionality of these classes is in the \c operator().
 *  This operator implements the complete TPS-RBF function. In this case:
 * \f[ \phi(r) = r^4\log(r) \f]
 *  where
 *  \f$ r = |x - x_j| \f$ in 1D, 
 *  \f$ r = \sqrt{(x - x_j)^2 + (y - y_j)^2} \f$ in 2D, and
 *  \f$ r = \sqrt{(x - x_j)^2 + (y - y_j)^2 + (z - z_j)^2} \f$ in 3D.
 *
 *  \author  Luis M. de la Cruz [ Mon Nov 12 14:02:37 GMT 2007 ]
 ***************************************************************************
 */
template<typename T, int Dim> class TPS { };

/*! 
 ***************************************************************************
 *  Thin Plate Spline specialization for 2D.
 *  \author  Luis M. de la Cruz [ Mon Nov 12 14:02:37 GMT 2007 ]
 ***************************************************************************
 */
template<typename T>
class TPS<T, 2> {
public:
    T operator() (T x, T y, T xj, T yj) {
	T r = sqrt( (x-xj) * (x-xj) + (y-yj) * (y-yj) );  
	if (r == 0) return 0;
	else return pow(r, 4) * log(r);	
    }
};
/*! 
 ***************************************************************************
 *  Generic functor to evaluate the 1st derivative of TPS with respect to x.
 *  All the functionality of these classes are in the \c operator().
 *  This operator implements the first derivative of TPS-RBF function with
 *  respect to x. In this case:
 *  \f[ \frac{\partial \phi(r)}{\partial x} = 
 *                                2 r^2(x - x_j)log(r^2) + r^2(x - x_j) \f]
 *  where
 *  \f$ r = |x - x_j| \f$ in 1D, 
 *  \f$ r = \sqrt{(x - x_j)^2 + (y - y_j)^2} \f$ in 2D, and
 *  \f$ r = \sqrt{(x - x_j)^2 + (y - y_j)^2 + (z - z_j)^2} \f$ in 3D.
 *
 *  \author  Luis M. de la Cruz [ Mon Nov 12 14:02:37 GMT 2007 ]
 ***************************************************************************
 */
template<typename T, int Dim> class TPS_1DX { };
/*! 
 ***************************************************************************
 *  1st derivative of TPS with respect to x for 2D.
 *  \author  Luis M. de la Cruz [ Mo ]
 ***************************************************************************
 */
template<typename T>
class TPS_1DX<T, 2> {
public:
    T operator() (T x, T y, T xj, T yj) {
	T r2, x_xj, y_yj;	
	x_xj = x - xj;	    
	y_yj = y - yj;
	r2 = x_xj * x_xj + y_yj * y_yj;	
	if( r2 == 0)
	    return 0.0;
	else {    
	    return  2 * r2 * x_xj * log(r2) + r2 * x_xj;
	}	
    }
};
/*! 
 ***************************************************************************
 *  Generic functor to evaluate the 1st derivative of TPS with respect to y.
 *  All the functionality of these classes are in the \c operator().
 *  This operator implements the first derivative of TPS-RBF function with
 *  respect to y. In this case:
 *  \f[ \frac{\partial \phi(r)}{\partial y} =
 *                                2 r^2(y - y_j)log(r^2) + r^2(y - y_j) \f]
 *  where
 *  \f$ r = \sqrt{(x - x_j)^2 + (y - y_j)^2} \f$ in 2D, and
 *  \f$ r = \sqrt{(x - x_j)^2 + (y - y_j)^2 + (z - z_j)^2} \f$ in 3D.
 *
 *  \author  Luis M. de la Cruz [ Mon Nov 12 14:02:37 GMT 2007 ]
 ***************************************************************************
 */
template<typename T, int Dim> class TPS_1DY { };
/*! 
 ***************************************************************************
 *  1st derivative of TPS with respect to y for 2D.
 *  \author  Luis M. de la Cruz [ Fri Apr 25 13:28:47 BST 2008 ]
 ***************************************************************************
 */
template<typename T>
class TPS_1DY<T, 2> {
public:
    T operator() (T x, T y, T xj, T yj) {
	T r2, x_xj, y_yj;	
	x_xj = x - xj;	    
	y_yj = y - yj;
	r2 = x_xj * x_xj + y_yj * y_yj;	
	if( r2 == 0)
	    return 0.0;
	else {    
	    return  2 * r2 * y_yj * log(r2) + r2 * y_yj;
	}	
    }
};

/*! 
 ***************************************************************************
 *  Generic functor to evaluate the 2nd derivative of TPS with respect to x.
 *  All the functionality of these classes are in the \c operator().
 *  This operator implements the second derivative of TPS-RBF function with
 *  respect to x. In this case:
 *  \f[ \frac{\partial^2 \phi(r)}{\partial x^2} = 
 *       4(x-x_i)^2\log(r^2) + 2r^2\log(r^2) + 6(x - x_j)^2 + r^2 \f]
 *  where
 *  \f$ r = |x - x_j| \f$ in 1D, 
 *  \f$ r = \sqrt{(x - x_j)^2 + (y - y_j)^2} \f$ in 2D, and
 *  \f$ r = \sqrt{(x - x_j)^2 + (y - y_j)^2 + (z - z_j)^2} \f$ in 3D.
 *
 *  \author  Luis M. de la Cruz [ Mon Nov 12 14:02:37 GMT 2007 ]
 ***************************************************************************
 */
template<typename T, int Dim> class TPS_2DX { };
/*! 
 ***************************************************************************
 *  2nd derivative of TPS with respect to x for 2D.
 *  \author  Luis M. de la Cruz [ Mon Nov 12 14:02:37 GMT 2007 ]
 ***************************************************************************
 */
template<typename T>
class TPS_2DX<T, 2> {
public:
    T operator() (T x, T y, T xj, T yj) {
	T r2, x_xj_2 ;
	
	x_xj_2 = (x - xj) * (x - xj);	    
	r2 = (x - xj) * (x - xj) + (y - yj) * (y - yj);	
	if( r2 == 0)
	    return 0.0;
	else {    
	    return  4 * x_xj_2 * log(r2) + 2 * r2 * log(r2) + 6 * x_xj_2 + r2;
	}	
    }  
};    
/*! 
 ***************************************************************************
 *  Generic functor to evaluate the 2nd derivative of TPS with respect to y.
 *  All the functionality of these classes are in the \c operator().
 *  This operator implements the second derivative of TPS-RBF function with
 *  respect to y. In this case:
 *  \f[ \frac{\partial^2 \phi(r)}{\partial y^2} = 
 *       4(y-y_i)^2\log(r^2) + 2r^2\log(r^2) + 6(y - y_j)^2 + r^2 \f]
 *  where
 *  \f$ r = \sqrt{(x - x_j)^2 + (y - y_j)^2} \f$ in 2D, and
 *  \f$ r = \sqrt{(x - x_j)^2 + (y - y_j)^2 + (z - z_j)^2} \f$ in 3D.
 *
 *  \author  Luis M. de la Cruz [ Mon Nov 12 14:02:37 GMT 2007 ]
 ***************************************************************************
 */
template<typename T, int Dim> class TPS_2DY { };
/*! 
 ***************************************************************************
 *  2nd derivative of TPS with respect to y for 2D.
 *  \author  Luis M. de la Cruz [ Mon Nov 12 14:02:37 GMT 2007 ]
 ***************************************************************************
 */
template<typename T>
class TPS_2DY<T, 2> {
public:
    T operator() (T x, T y, T xj, T yj) {
	T r2, y_yj_2 ;

	y_yj_2 = (y - yj) * (y - yj);
	r2 = (x - xj) * (x - xj) + (y - yj) * (y - yj);
	if( r2 == 0)
	    return 0.0;
	else {    
	    return  4 * y_yj_2 * log(r2) + 2 * r2 * log(r2) + 6 * y_yj_2 + r2;
	}	
    }  
};    

} // RBF namespace
 

#endif

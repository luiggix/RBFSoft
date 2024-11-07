
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

#ifndef _MQ_HPP_
#define _MQ_HPP_

#include <cmath>

namespace RBF {

/*! 
 ***************************************************************************
 *  Generic functor to evaluate Multiquadrics.
 *  All the functionality of these classes is in the \c operator().
 *  This operator implements the complete MQ-RBF function. In this case:
 *  \f[ \phi(r) = \sqrt{r^2+c^2} \f]
 *  where \f$ c \f$ is the shape parameter, and the \f$r\f$ is defined as
 *  \f$r = |x - x_j|\f$ in 1D, 
 *  \f$r = \sqrt{(x - x_j)^2 + (y - y_j)^2} \f$ in 2D, and
 *  \f$r = \sqrt{(x - x_j)^2 + (y - y_j)^2 + (z - z_j)^2} \f$ in 3D.
 *
 *  \author  Luis M. de la Cruz [ Mon Nov 12 14:02:37 GMT 2007 ]
 ***************************************************************************
 */
template<typename T, int Dim> class MQ { };

/*! 
 ***************************************************************************
 *  Multiquadric specialization for 1D.
 *  \author  Luis M. de la Cruz [ Mon Nov 12 14:02:37 GMT 2007 ]
 ***************************************************************************
 */
template<typename T>  
class MQ<T, 1> {
    T c2; /*!< Squared shape parameter */
 public:
    MQ(T c) : c2(c * c) { };

    inline 
    T operator() (T x, T xj) {
	return sqrt( (x-xj) * (x-xj) + c2 );
    }
};

/*! 
 ***************************************************************************
 *  Multiquadric specialization for 2D.
 *  \author  Luis M. de la Cruz [ Mon Nov 12 14:02:37 GMT 2007 ]
 ***************************************************************************
 */
template<typename T>    
class MQ<T, 2> { 
    T c2; /*!< Squared shape parameter */
public:
    MQ(T c) : c2(c * c) { };

    inline 
    T operator() (T x, T y, T xj, T yj) {
	return sqrt( (x-xj) * (x-xj) + (y-yj) * (y-yj) + c2 );
    }
};

/*! 
 ***************************************************************************
 *  Multiquadric specialization for 3D.
 *  \author  Luis M. de la Cruz [ Mon Nov 12 14:02:37 GMT 2007 ]
 ***************************************************************************
 */
template<typename T>    
class MQ<T, 3> { 
    T c2; /*!< Squared shape parameter */
public:
    MQ(T c) : c2(c * c) { };
    
    inline T operator() (T x, T y, T z, T xj, T yj, T zj) {
	return sqrt((x-xj) * (x-xj) + (y-yj) * (y-yj) + (z-zj) * (z-zj) + c2);
    }
};

/*! 
 ***************************************************************************
 *  Generic functor to evaluate the 1st derivative of MQ with respect to x.
 *  All the functionality of these classes are in the \c operator().
 *  This operator implements the first derivative of MQ-RBF function with
 *  respect to x. In this case:
 *  \f[ \frac{\partial \phi(r)}{\partial x} =
 *      \frac{(x-x_j)}{\sqrt{r^2+c^2}} \f]
 *  where
 *  \f$ c \f$ is the shape parameter,
 *  \f$ r = |x - x_j| \f$ in 1D, 
 *  \f$ r = \sqrt{(x - x_j)^2 + (y - y_j)^2} \f$ in 2D, and
 *  \f$ r = \sqrt{(x - x_j)^2 + (y - y_j)^2 + (z - z_j)^2} \f$ in 3D.
 *
 *  \author  Luis M. de la Cruz [ Mon Nov 12 14:02:37 GMT 2007 ]
 ***************************************************************************
 */
template<typename T, int Dim> class MQ_1DX { };

/*! 
 ***************************************************************************
 *  1st derivative of MQ with respect to x for 1D.
 *  \author  Luis M. de la Cruz [ Mon Nov 12 14:02:37 GMT 2007 ]
 ***************************************************************************
 */
template<typename T>
class MQ_1DX<T, 1> {
    T c2; /*!< Squared shape parameter */
public:
    MQ_1DX(T c) : c2(c * c) { };
    
    inline T operator() (T x, T xj) {
	return (x-xj) / sqrt( (x-xj) * (x-xj) + c2 );
    }
};

/*! 
 ***************************************************************************
 *  1st derivative of MQ with respect to x for 2D.
 *  \author  Luis M. de la Cruz [ Mon Nov 12 14:02:37 GMT 2007 ]
 ***************************************************************************
 */
template<typename T>
class MQ_1DX<T, 2> {
    T c2; /*!< Squared shape parameter */
public:
    MQ_1DX(T c) : c2(c * c) { };
    
    inline T operator() (T x, T y, T xj, T yj) {
	return (x-xj) / sqrt( (x-xj) * (x-xj) + (y-yj) * (y-yj) + c2 );
    }
};

/*! 
 ***************************************************************************
 *  1st derivative of MQ with respect to x for 3D.
 *  \author  Luis M. de la Cruz [ Mon Nov 12 14:02:37 GMT 2007 ]
 ***************************************************************************
 */
template<typename T>
class MQ_1DX<T, 3> {
    T c2; /*!< Squared shape parameter */
public:
    MQ_1DX(T c) : c2(c * c) { };
    
    inline T operator() (T x, T y, T z, T xj, T yj, T zj) {
	return (x-xj) / 
	    sqrt( (x-xj) * (x-xj) + (y-yj) * (y-yj) + (z-zj) * (z-zj) + c2 );
    }
};

/*! 
 ***************************************************************************
 *  Generic functor to evaluate the 2nd derivative of MQ with respect to x.
 *  All the functionality of these classes are in the \c operator().
 *  This operator implements the second derivative of MQ-RBF function with
 *  respect to x. In this case:
 *  \f[ \frac{\partial^2 \phi(r)}{\partial x^2} =
 *      \frac{r^2 + c^2 - (x-x_j)^2 }{(r^2+c^2)^{3/2}} \f]
 *  where
 *  \f$ c \f$ is the shape parameter, 
 *  \f$ r = |x - x_j| \f$ in 1D, 
 *  \f$ r = \sqrt{(x - x_j)^2 + (y - y_j)^2} \f$ in 2D, and
 *  \f$ r = \sqrt{(x - x_j)^2 + (y - y_j)^2 + (z - z_j)^2} \f$ in 3D.
 *
 *  \author  Luis M. de la Cruz [ Mon Nov 12 14:02:37 GMT 2007 ]
 ***************************************************************************
 */
template<typename T, int Dim> class MQ_2DX { };

/*! 
 ***************************************************************************
 *  2nd derivative of MQ with respect to x for 1D.
 *  \author  Luis M. de la Cruz [ Mon Nov 12 14:02:37 GMT 2007 ]
 ***************************************************************************
 */
template<typename T>  
class MQ_2DX<T, 1> {
    T c2; /*!< Squared shape parameter */
public:
    MQ_2DX(T c) : c2(c * c) { };

    inline T operator() (T x, T xj) {
	T r2 = (x-xj) * (x-xj);  
	return  c2 / ( sqrt(r2 + c2) * (r2 + c2) ); 
    }  
};

/*! 
 ***************************************************************************
 *  2nd derivative of MQ with respect to x for 2D.
 *  \author  Luis M. de la Cruz [ Mon Nov 12 14:02:37 GMT 2007 ]
 ***************************************************************************
 */
template<typename T>
class MQ_2DX<T, 2> {
    T c2; /*!< Squared shape parameter */
public:
    MQ_2DX(T c) : c2(c * c) { };

    inline T operator() (T x, T y, T xj, T yj) {
	T r2 = (x-xj) * (x-xj) + (y-yj) * (y-yj);  
	return  ((y - yj) * (y - yj) + c2) / 
	    (sqrt(r2 + c2) * (r2 + c2));
    }  
};

/*! 
 ***************************************************************************
 *  2nd derivative of MQ with respect to x for 3D.
 *  \author  Luis M. de la Cruz [ Mon Nov 12 14:02:37 GMT 2007 ]
 ***************************************************************************
 */
template<typename T>
class MQ_2DX<T, 3> {
    T c2; /*!< Squared shape parameter */
public:
    MQ_2DX(T c) : c2(c * c) { };

    inline T operator() (T x, T y, T z, T xj, T yj, T zj) {
	T r2 = (x-xj) * (x-xj) + (y-yj) * (y-yj) + (z-zj) * (z-zj);  
	return  ((y - yj) * (y - yj) + (z-zj) * (z-zj) + c2) / 
	    (sqrt(r2 + c2) * (r2 + c2)); 
    }  
};

/*! 
 ***************************************************************************
 *  Generic functor to evaluate the 1st derivative of MQ with respect to y.
 *  All the functionality of these classes are in the \c operator().
 *  This operator implements the first derivative of MQ-RBF function with
 *  respect to y. In this case:
 *  \f[ \frac{\partial \phi(r)}{\partial x} =
 *      \frac{(y-y_j)}{\sqrt{r^2+c^2}} \f]
 *  where
 *  \f$ c \f$ is the shape parameter,
 *  \f$ r = \sqrt{(x - x_j)^2 + (y - y_j)^2} \f$ in 2D, and
 *  \f$ r = \sqrt{(x - x_j)^2 + (y - y_j)^2 + (z - z_j)^2} \f$ in 3D.
 *  \par NOTE: The 1D version is not needed in this case.
 *
 *  \author  Luis M. de la Cruz [ Mon Nov 12 14:02:37 GMT 2007 ]
 ***************************************************************************
 */
template<typename T, int Dim> class MQ_1DY { };

 /*! 
 ***************************************************************************
 *  1st derivative of MQ with respect to y for 2D.
 *  \author  Luis M. de la Cruz [ Mon Nov 12 14:02:37 GMT 2007 ]
 ***************************************************************************
 */ 
template<typename T>
class MQ_1DY<T, 2> {
    T c2;
public:
    MQ_1DY(T c) : c2(c * c) { };
    
    inline T operator() (T x, T y, T xj, T yj) {
	T r2 = (x-xj) * (x-xj) + (y-yj) * (y-yj);
	return (y - yj) / sqrt(r2 + c2);
    }
};

/*! 
 ***************************************************************************
 *  1st derivative of MQ with respect to y for 3D.
 *  \author  Luis M. de la Cruz [ Mon Nov 12 14:02:37 GMT 2007 ]
 ***************************************************************************
 */
template<typename T>
class MQ_1DY<T, 3> {
    T c2; /*!< Squared shape parameter */
public:
    MQ_1DY(T c) : c2(c * c) { };
    
    inline T operator() (T x, T y, T z, T xj, T yj, T zj) {
	T r2 = (x-xj) * (x-xj) + (y-yj) * (y-yj) + (z-zj) * (z-zj);
	return (y - yj) / sqrt(r2 + c2);
    }
};
/*! 
 ***************************************************************************
 *  Generic functor to evaluate the 2nd derivative of MQ with respect to y.
 *  All the functionality of these classes are in the \c operator().
 *  This operator implements the second derivative of MQ-RBF function with
 *  respect to y. In this case:
 *  \f[ \frac{\partial^2 \phi(r)}{\partial x^2} =
 *      \frac{r^2 + c^2 - (y-y_j)^2}{(r^2+c^2)^{3/2}} \f]
 *  where
 *  \f$ c \f$ is the shape parameter,
 *  \f$ r = \sqrt{(x - x_j)^2 + (y - y_j)^2} \f$ in 2D, and
 *  \f$ r = \sqrt{(x - x_j)^2 + (y - y_j)^2 + (z - z_j)^2} \f$ in 3D.
 *  \par NOTE: The 1D version is not needed in this case.
 *
 *  \author  Luis M. de la Cruz [ Mon Nov 12 14:02:37 GMT 2007 ]
 ***************************************************************************
 */
template<typename T, int Dim> class MQ_2DY { };

/*! 
 ***************************************************************************
 *  2nd derivative of MQ with respect to y for 2D.
 *  \author  Luis M. de la Cruz [ Mon Nov 12 14:02:37 GMT 2007 ]
 ***************************************************************************
 */
template<typename T>
class MQ_2DY<T, 2> {
    T c2;
public:
    MQ_2DY(T c) : c2(c * c) { };

    inline 
    T operator() (T x, T y, T xj, T yj) {
	T r2 = (x-xj) * (x-xj) + (y-yj) * (y-yj);
	return  ((x - xj) * (x - xj) + c2) / (sqrt(r2 + c2) * (r2 + c2)) ;
    }  
};


/*! 
 ***************************************************************************
 *  2nd derivative of MQ with respect to y for 3D.
 *  \author  Luis M. de la Cruz [ Mon Nov 12 14:02:37 GMT 2007 ]
 ***************************************************************************
 */
template<typename T>
class MQ_2DY<T, 3> {
    T c2; /*!< Squared shape parameter */
public:
    MQ_2DY(T c) : c2(c * c) { };

    inline T operator() (T x, T y, T z, T xj, T yj, T zj) {
	T r2 = (x-xj) * (x-xj) + (y-yj) * (y-yj) + (z-zj) * (z-zj);  
	return  ((x - xj) * (x - xj) + (z-zj) * (z-zj) + c2) /
	    sqrt(r2 + c2) * (r2 + c2);
    }  
};

/*! 
 ***************************************************************************
 *  Generic functor to evaluate the 1st derivative of MQ with respect to z.
 *  All the functionality of these classes are in the \c operator().
 *  This operator implements the first derivative of MQ-RBF function with
 *  respect to z. In this case:
 *  \f[ \frac{\partial \phi(r)}{\partial x} =
 *      \frac{(z-z_j)}{\sqrt{r^2+c^2}} \f]
 *  where
 *  \f$ c \f$ is the shape parameter,
 *  \f$ r = \sqrt{(x - x_j)^2 + (y - y_j)^2 + (z - z_j)^2} \f$ in 3D.
 *  \par NOTE: The 1D and 2D versions are not needed in this case.
 *
 *  \author  Luis M. de la Cruz [ Mon Nov 12 14:02:37 GMT 2007 ]
 ***************************************************************************
 */
template<typename T, int Dim> class MQ_1DZ { };
/*! 
 ***************************************************************************
 *  1st derivative of MQ with respect to z for 3D.
 *  \author  Luis M. de la Cruz [ Mon Nov 12 14:02:37 GMT 2007 ]
 ***************************************************************************
 */
template<typename T>
class MQ_1DZ<T, 3> {
    T c2; /*!< Squared shape parameter */
public:
    MQ_1DZ(T c) : c2(c * c) { };
    
    inline T operator() (T x, T y, T z, T xj, T yj, T zj) {
	T r2 = (x-xj) * (x-xj) + (y-yj) * (y-yj) + (z-zj) * (z-zj);
	return (z - zj) / sqrt(r2 + c2);
    }
};
/*! 
 ***************************************************************************
 *  Generic functor to evaluate the 2nd derivative of MQ with respect to z.
 *  All the functionality of these classes are in the \c operator().
 *  This operator implements the second derivative of MQ-RBF function with
 *  respect to z. In this case:
 *  \f[ \frac{\partial^2 \phi(r)}{\partial x^2} =
 *      \frac{r^2 + c^2 - (z-z_j)^2}{(r^2+c^2)^{3/2}} \f]
 *  where
 *  \f$ c \f$ is the shape parameter,
 *  \f$ r = \sqrt{(x - x_j)^2 + (y - y_j)^2 + (z - z_j)^2} \f$ in 3D.
 *  \par NOTE: The 1D and 2D versions are not needed in this case.
 *
 *  \author  Luis M. de la Cruz [ Mon Nov 12 14:02:37 GMT 2007 ]
 ***************************************************************************
 */

template<typename T, int Dim> class MQ_2DZ { };

/*! 
 ***************************************************************************
 *  2nd derivative of MQ with respect to z for 3D.
 *  \author  Luis M. de la Cruz [ Mon Nov 12 14:02:37 GMT 2007 ]
 ***************************************************************************
 */
template<typename T>
class MQ_2DZ<T, 3> {
    T c2; /*!< Squared shape parameter */
public:
    MQ_2DZ(T c) : c2(c * c) { };

    inline T operator() (T x, T y, T z, T xj, T yj, T zj) {
	T r2 = (x-xj) * (x-xj) + (y-yj) * (y-yj) + (z-zj) * (z-zj);  
	return  ((x - xj) * (x - xj) + (y-yj) * (y-yj) + c2) /
	    ( sqrt(r2 + c2) * (r2 + c2) );
    }  
};

}  // RBF Namespace

#endif // _MQ_HPP_

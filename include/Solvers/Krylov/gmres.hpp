
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

#ifndef __INCLUDE_GMRES_H__
#define __INCLUDE_GMRES_H__

/*! 
 ***************************************************************************
 *  \file gmres.hpp
 *  \brief   GMRES IMPLEMENTATIONS
 *  \author  Luis M. de la Cruz
 *  \date Mon Nov 12 11:52:45 GMT 2007
 ***************************************************************************
 */

#include <cmath>

using namespace std;

/*!
 *  GMRES, preconditioners, based on FLENS.
 */
namespace flens {

/*! 
 ***************************************************************************
 *  Type definition for GMRES.
 *  \author  Luis M. de la Cruz [ Wed Nov 28 10:06:43 GMT 2007 ]
 ***************************************************************************
 */
template <typename A>
struct gmres_
{
    typedef A AuxVector;
    typedef typename A::ElementType T;
};
/*! 
 ***************************************************************************
 *  Dense Vector specialization.
 *  \author  Luis M. de la Cruz [ Wed Nov 28 10:06:43 GMT 2007 ]
 ***************************************************************************
 */
template <typename I>
struct gmres_<DenseVector<I> >
{
    typedef typename DenseVector<I>::NoView AuxVector;
    typedef typename DenseVector<I>::ElementType T;
};

/*! 
 ***************************************************************************
 *  GMRES NON PRECONDITIONED.
 *  \author  Luis M. de la Cruz [ Mon Nov 12 11:52:45 GMT 2007 ]
 *  \param A coefficients matrix.
 *  \param x solution vector.
 *  \param b right-hand vector.
 *  \param maxiter maximum number or GMRES iterations.
 *  \param tol tolerance.
 *  \return number of GMRES iterations to achieve the prescribed tolerance.
 ***************************************************************************
 */
template <class Mat, class Vec>
int gmres ( Mat const & A, 
	    Vec & x, 
	    const Vec & b, 
	    unsigned int const maxiter, 
	    double const tol) 
{
    typedef GeMatrix<FullStorage<double, ColMajor> >  GEMatrix;
    typedef typename gmres_<Vec>::AuxVector vec_t;
    typedef typename gmres_<Vec>::T prec_t;

    unsigned int N = b.length();

    GEMatrix  V(N          , maxiter + 1);  // Orthonormal basis
    GEMatrix  H(maxiter + 1, maxiter    );  // Hessemberg Matrix

    vec_t v_kp1(N);        // v_k+1
    vec_t r(N);            // residual
    vec_t g(maxiter + 1);  // rhs of R * y = g 
    vec_t c(maxiter + 1), s(maxiter + 1);  // Given's rotations coefficients

    prec_t tol_2 = 10e-6;
    prec_t beta, nu, rho = tol + 1;
    prec_t Htemp, h_ii;
    
    r = b - A * x;       // Initial residual
    beta = ::sqrt(r * r);  // || r ||
    V(_,1) = r / beta;   // First vector of the orthonormal base: v_1
    g(1) = beta;         // First element of rhs of: || g - R_k * y_k ||

// ----- GMRES MAIN LOOP : BEGIN -----
    unsigned int k = 0;
    while (k < maxiter && rho > tol){
	++k;
	v_kp1 = A * V(_,k);

// ----- Arnoldi's Orthogonalization with Modified Gram-Schmidt
        for (unsigned int i = 1; i <= k; ++i) {
            H(i, k) = v_kp1 * V(_,i);
            v_kp1 -=  H(i, k) * V(_, i);
        }
	H(k+1, k) = ::sqrt( v_kp1 * v_kp1 ); // || v_kp1 ||

// ----- Check orthogonality : ||A * v_k|| + d * ||v_k+1|| = ||A * v_k||
	if ( H(k+1, k) < tol_2) {
	    for (unsigned int i = 1; i <= k; ++i) {
		Htemp = v_kp1 * V(_, i) ;
		v_kp1 -= Htemp * V(_, i);
	    }
	    H(k+1, k) = ::sqrt( v_kp1 * v_kp1 ); // || v_kp1 ||
	}

	V(_, k+1) = v_kp1 / H(k+1, k);     // v_k+1 of the orthonormal base

// ----- if k > 1 apply Gives rotations to the kth column of H
        for (unsigned int i = 1; i < k; ++i) {
            h_ii  = c(i) * H(i,k) + s(i) * H(i+1,k);
            H(i+1, k)  = s(i) * H(i,k) - c(i) * H(i+1,k);
            H(i  , k) = h_ii;
        }

// ----- Given's rotation coefficients
        nu  = ::sqrt( H(k, k) * H(k, k) + H(k+1, k) * H(k+1, k) );	
        c(k) = H(k  , k) / nu;
	s(k) = H(k+1, k) / nu;	
        H(k, k) = nu; 

// ----- Apply Given's rotation to rhs
        g(k+1) = s(k) * g(k);
        g(k  ) = c(k) * g(k);
        rho = g(k+1) / beta; // g(k+1) == || beta * e_1 - H_m * y || 
                
    } 
// ----- GMRES MAIN LOOP : END -----

    c = 0; // Re-use c to hold y_k solution

// ----- Solve the triangular R * y_k = g system (R = H)
    if (k > 1 ) {
	for (int i = k; i >= 1; --i) {
	    c(i) = g(i) / H(i,i);
	    for (int j = k; j > i; --j)
		c(i) -= H(i,j) * c(j) / H(i,i);
	}
// ----- Update the solution x_k = x_0 + V_k * y_k (recall y = c)
	x += V * c;
    }
    r = b - A * x;
    double residual = ::sqrt ( r * r);
//    std::cout << "\n r = b - Ax = " <<  residual << std::endl;  

    return k;
}

/*! 
 ***************************************************************************
 *  GMRES PRECONDITIONED.
 *  \author  Luis M. de la Cruz [ Wed Nov 28 10:10:09 GMT 2007 ]
 *  \param A coefficients matrix.
 *  \param x solution vector.
 *  \param b right-hand vector.
 *  \param M preconditioner.
 *  \param maxiter maximum number or GMRES iterations.
 *  \param tol tolerance.
 *  \return number of GMRES iterations to achieve the prescribed tolerance.
 ***************************************************************************
 */
template <class Mat, class Vec, class Precond>
int gmres ( Mat const & A, 
	    Vec & x, 
	    const Vec & b, 
	    const Precond & M,
	    unsigned int const maxiter, 
	    double const tol) 
{
    typedef GeMatrix<FullStorage<double, ColMajor> >  GEMatrix;
    typedef typename gmres_<Vec>::AuxVector vec_t;
    typedef typename gmres_<Vec>::T prec_t;

    unsigned int N = b.length();

    GEMatrix  V(N          , maxiter + 1);  // Orthonormal basis
    GEMatrix  H(maxiter + 1, maxiter    );  // Hessemberg Matrix

    vec_t v_kp1(N);        // v_k+1
    vec_t r(N);            // residual
    vec_t g(maxiter + 1);  // rhs of R * y = g 
    vec_t c(maxiter + 1), s(maxiter + 1);  // Given's rotations coefficients

    prec_t tol_2 = 10e-6;
    prec_t beta, nu, rho = tol + 1;
    prec_t Htemp, h_ii;
    
    r = M * ( b - A * x );       // Initial residual (preconditioned)

    beta = ::sqrt(r * r);  // || r ||
    V(_,1) = r / beta;   // First vector of the orthonormal base: v_1
    g(1) = beta;         // First element of rhs of: || g - R_k * y_k ||

// ----- GMRES MAIN LOOP : BEGIN -----
    unsigned int k = 0;
    while (k < maxiter && rho > tol){
	++k;
	v_kp1 = M * ( A * V(_,k) ); // Solve the preconditioned problem

// ----- Arnoldi's Orthogonalization with Modified Gram-Schmidt
        for (unsigned int i = 1; i <= k; ++i) {
            H(i, k) = v_kp1 * V(_,i);
            v_kp1 -=  H(i, k) * V(_, i);
        }
	H(k+1, k) = ::sqrt( v_kp1 * v_kp1 ); // || v_kp1 ||

// ----- Check orthogonality : ||A * v_k|| + d * ||v_k+1|| = ||A * v_k||
	if ( H(k+1, k) < tol_2) {
	    for (unsigned int i = 1; i <= k; ++i) {
		Htemp = v_kp1 * V(_, i) ;
		v_kp1 -= Htemp * V(_, i);
	    }
	    H(k+1, k) = ::sqrt( v_kp1 * v_kp1 ); // || v_kp1 ||
	}

	V(_, k+1) = v_kp1 / H(k+1, k);     // v_k+1 of the orthonormal base

// ----- if k > 1 apply Gives rotations to the kth column of H
        for (unsigned int i = 1; i < k; ++i) {
            h_ii  = c(i) * H(i,k) + s(i) * H(i+1,k);
            H(i+1, k)  = s(i) * H(i,k) - c(i) * H(i+1,k);
            H(i  , k) = h_ii;
        }

// ----- Given's rotation coefficients
        nu  = ::sqrt( H(k, k) * H(k, k) + H(k+1, k) * H(k+1, k) );	
        c(k) = H(k  , k) / nu;
	s(k) = H(k+1, k) / nu;	
        H(k, k) = nu; 

// ----- Apply Given's rotation to rhs
        g(k+1) = s(k) * g(k);
        g(k  ) = c(k) * g(k);
        rho = g(k+1) / beta; // g(k+1) == || beta * e_1 - H_m * y || 
                
    } 
// ----- GMRES MAIN LOOP : END -----

    c = 0; // Re-use c to hold y_k solution

// ----- Solve the triangular R * y_k = g system (R = H)
    if (k > 1 ) {
	for (int i = k; i >= 1; --i) {
	    c(i) = g(i) / H(i,i);
	    for (int j = k; j > i; --j)
		c(i) -= H(i,j) * c(j) / H(i,i);
	}
// ----- Update the solution x_k = x_0 + V_k * y_k (recall y = c)
	x += V * c;
    }
    r = b - A * x;
    double residual = ::sqrt ( r * r);
//    std::cout << "\n r = b - Ax = " <<  residual << std::endl;  

    return k;
}



} // Namespace flens

#endif

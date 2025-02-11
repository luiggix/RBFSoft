
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

#ifndef _GAUSS_H_
#define _GAUSS_H_

/*! 
 ***************************************************************************
 *  \namespace  Solver
 *  Linear system solvers.
 *  This namespace contains implementations for solving general linear 
 *  systems.
 *  \author  Luis M. de la Cruz [ Tue Sep 11 17:06:45 BST 2007 ]
 ***************************************************************************
 */

#include <algorithm>
#include "Traits.hpp"

namespace Solver {

/*! 
 *  Gauss elimination with partial pivoting algorithm for dense linear systems.
 *  \author Daniel Cervantes & Luis M. de la Cruz [ Tue Sep 11 17:06:45 BST 2007 ]
 *  \param A matrix of coefficients.
 *  \param b right term vector.
 *  \return the x vector solution of the linear system A * x = b.
 *  \todo Check the matrix access, remember that I'm using ColMajor
 *        for all the matrices (see FLENS manual). 
 */
template<typename T>
Vec Gauss(Mat& A, Vec &b) {

// FLENS specific : using numRows() to get the number of rows of the matrix
//    int N = A.numRows(); // Number of rows

    int N = A.rows(); // Number of rows

    Mat a(N, N+1); 
    Vec x(N);

    T  pivot, temp, q;  
    int ipivot,ip1;

    for(int i = 0; i < N; ++i)
	for(int j = 0; j < N; ++j) 
	    a(i, j) = A(i, j);

    for(int i = 0; i < N; ++i)
	a(i, N) = b(i);

    //    std::cout << "\n a = \n" << a << std::endl;


    for(int i = 0; i < N-1; ++i){
	pivot = 0.0;
	for(int j = i; j < N; ++j){      
	    temp = fabs(a(j,i));      
	    if(pivot < temp) {
		pivot = temp;
		ipivot = j;
	    }            
	}

	if(fabs(pivot) < 1e-10)
	{
	    std::cout << "\n Solver::Gauss : Error : Singular Matrix \n";
	    exit(1);
	}
	
	if(ipivot != i) 
	    for(int k = i; k < N; ++k)	
		std::swap( a(i,k), a(ipivot,k) ); // from STL algorithm
	
	ip1 = i + 1;
	for(int k = ip1; k < N; ++k){
	    q = -a(k, i) / a(i, i);
	    a(k, i) = 0.0;
	    for(int j = ip1; j < N; j++)
		a(k, j) += q * a(i, j);
	}                        
    }


    x(N-1) = a(N-1, N) / a(N-1, N-1);
    for(int k = N-2; k >= 0; --k) {
        q = 0.0;
        for(int j = k+1; j < N; ++j)
            q += a(k, j) * x(j);
        x(k) = (a(k, N+1) - q) / (a(k,k));
    }

    return x;  
} 

} // Solver namespace

#endif //_GAUSS_H_

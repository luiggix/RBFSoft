
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

#ifndef __INCLUDE_ACBFPREC_H__
#define __INCLUDE_ACBFPREC_H__

#include "Solvers/Gauss.hpp"

using namespace std;

namespace flens {

// ACBFPrec declaration.
template <typename MA, typename TD> class ACBFPrec;

template <typename MA, typename TD>
struct TypeInfo<ACBFPrec<MA, TD> >
{
    typedef ACBFPrec<MA, TD>  Impl;
    typedef double    ElementType;
};

/*! 
 ***************************************************************************
 *  ACBF Preconditioner.
 *  \author  Luis M. de la Cruz [ Tue Dec 11 12:54:41 GMT 2007 ]
 ***************************************************************************
 */ 
template <typename MA, typename TD>
class ACBFPrec : public SymmetricMatrix<ACBFPrec<MA, TD> >
{

public:
    ACBFPrec(const MA &_A) : A(_A) { }

/*!
 *  This function constructs the preconditioner as defined in Kansa et al.
 *  [Adv. Comp. Math., 2005, 23, (1-2), 31-54].
 */
    bool construct(TD domain, int max_knots, int sp = 0, 
		   bool print_screen = false) 
	{
	    int N = A.numRows();	
	    
	    timer time;
	    double time_neigh = 0, time_Bi = 0, time_sol = 0, time_kdtree = 0;
	    int irow;

	    Mat B(max_knots+sp, N), 
		BT(N, max_knots+sp), 
		Brhs(N, 1),
		C(max_knots+sp, max_knots+sp);
	    
	    W.resize(max_knots+sp, N);
	    P.resize(max_knots+sp, N);
	    
	    Vec w(N), eaux(N);

	    time.tic();
	    domain.calcKDTree();
	    time_kdtree = time.toc();

//	    std::ofstream precfile ("precond.dat");	    

	    for(int j = W.firstCol(); j <= W.lastCol() ; ++j) {

//		std::cout << "\n Knot " << j;

		time.tic();

		domain.findNeighbors(j, P, max_knots, sp);
		time_neigh += time.toc();

//		precfile << "\n" << j << "\t";

		time.tic();
		for(int im = 1; im <= max_knots+sp; ++im) {
//		    precfile << P(im, j) << "  ";
// 		    precfile << j << "\t" << N - P(im, j) << "\n";
		    irow = P(im, j);
// For Gauss 
/* */
		    B(im, _ ) = A( irow , _ );
		    BT(_, im ) = B(im, _);
/* */
// For QR-SVD
/* *
		    BT(_, im ) = A( irow , _ );
/* */
		}
		time_Bi += time.toc();
/*!
 * \todo Optimize the Matrix-Matrix multiplication
 */
/* Solve the least square problem -----  BEGIN L-NE ------- */
		time.tic();
		C = B * BT;
		eaux =  B( _ , j);  // = B * e;
		W(_, j) = Solver::Gauss<prec_t> ( C, eaux );
		time_sol += time.toc();
/* Solve the least square problem -----  END  L-NE ------- */


/*!
 * \todo Optimize this process (L-QR-SVD), basically I need to 
 *       incorporate the algorithm described in Kansa & Ling, 
 *       Adv. Comp. Math. 2005, 23: 31-54, pag 39.
 */
/* Solve the least square problem -----  BEGIN L-QR or SVD -- *
		time.tic();
		Brhs = 0.0;
		Brhs(j, 1) =  1.0;
		ls(NoTrans, BT, Brhs); // QR
//		lss(BT, Brhs);         // SVD
		W(_, j) = Brhs(_(1,max_knots+sp), 1);
		time_sol += time.toc();
/* Solve the least square problem -----  END  L-QR or SVD -- */

	    }
//	    precfile.close();

	    if( print_screen ) {
		std::cout << " \n | KDTree constr = " << time_kdtree;
		std::cout << " \n | Find neigbors = " << time_neigh;
		std::cout << " \n | B_i matrices = " << time_Bi;
		std::cout << " \n | least-square = " << time_sol;
	    }
	}

    const MA  &A; ///< Coefficient matrix.
    iMat P; ///< Index matrix.
    MA W;  ///< Preconditioner coefficients.
};
  
/*! 
 ***************************************************************************
 *  Multiplication definition for the ACBF preconditioner.
 *  Here the multiplaction \c alpha * \c J * \c b + \c beta * \c x is defined.
 *  In this case \c alpha = 1 and \c beta = 0, so we do \c x = \c J * \c b .
 *  \param alpha scalar.
 *  \param J ACBF preconditioner.
 *  \param b vector.
 *  \param beta scalar.
 *  \param x vector where the result is stored.
 *  \author  Luis M. de la Cruz [ Tue Dec 11 12:54:41 GMT 2007 ]
 ***************************************************************************
 */ 
template<class Mat, class TD>
void
mv(double alpha, 
   const ACBFPrec<Mat, TD> &ACBF, 
   const DenseVector<Array<double> > &b,
   double beta, 
   DenseVector<Array<double> > &x)
{
    assert(alpha==1.);
    assert(beta==0.);
        
    if ( x.firstIndex() != b.firstIndex() )
	x.resize( b.length(), b.firstIndex() );
    
    int bi = x.firstIndex();
    int ei = x.lastIndex();
    int bj = ACBF.P.firstRow(); 
    int ej = ACBF.P.lastRow();

    int col;
    for (int i = bi; i <= ei; ++i) {
	x(i) = 0.0;
	for(int j = bj; j <= ej; ++j) {
	    col = ACBF.P(j,i);
	    x(i) += b(col) * ACBF.W(j,i);
	}
    }
}



} // Namespace FLENS

#endif

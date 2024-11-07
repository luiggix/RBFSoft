
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

#ifndef __INCLUDE_DIAGONALPREC_H__
#define __INCLUDE_DIAGONALPREC_H__

using namespace std;

namespace flens {

// DiagonalPrec declaration.
template <typename MA> class DiagonalPrec;

template <typename MA>
struct TypeInfo<DiagonalPrec<MA> >
{
    typedef DiagonalPrec<MA>  Impl;
    typedef double    ElementType;
};

/*! 
 ***************************************************************************
 *  Diagonal Preconditioner.
 *  \author  Luis M. de la Cruz [ Wed Nov 28 10:06:43 GMT 2007 ]
 ***************************************************************************
 */  
template <typename MA>
class DiagonalPrec : public SymmetricMatrix<DiagonalPrec<MA> >
{
public:
    DiagonalPrec(const MA &_A, const double eps = 1.0e-14) 
	: A(_A), EPS(eps) { }
    const MA  &A;  ///< Coefficient matrix.
    const double EPS; ///< almost zero.
};

/*! 
 ***************************************************************************
 *  Multiplication definition for the Diagonal preconditioner.
 *  Here the multiplaction \c alpha * \c J * \c b + \c beta * \c x is defined.
 *  In this case \c alpha = 1 and \c beta = 0, so we do \c x = \c J * \c b, but
 *  because this operation is intended for the Diagonal preconditioner, we
 *  finally do \c x = \c b / \c diag(A).
 *  \param alpha scalar.
 *  \param J Diagonal preconditioner.
 *  \param b vector.
 *  \param beta scalar.
 *  \param x vector where the result is stored.
 *  \author  Luis M. de la Cruz [ Wed Nov 28 10:06:43 GMT 2007 ]
 ***************************************************************************
 */   
template<class Mat>
void
mv(double alpha, 
   const DiagonalPrec<Mat> &J, 
   const DenseVector<Array<double> > &b,
   double beta, 
   DenseVector<Array<double> > &x)
{
    assert(alpha==1.);
    assert(beta==0.);
        
    if ( x.firstIndex() != b.firstIndex() )
	x.resize( b.length(), b.firstIndex() );
    
    for (int i = x.firstIndex(); i <= x.lastIndex(); ++i) {
	if (J.A(i,i) > J.EPS)
	    x(i) = b(i) / J.A(i,i);
	else
	    x(i) = 0.0;
    }
}

} // Namespace FLENS

#endif

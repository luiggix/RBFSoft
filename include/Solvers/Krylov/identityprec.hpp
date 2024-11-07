
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

#ifndef __INCLUDE_IDENTITYPREC_H__
#define __INCLUDE_IDENTITYPREC_H__

using namespace std;

namespace flens {

// IdentityPrec declaration.
template <typename MA> class IdentityPrec;

template <typename MA>
struct TypeInfo<IdentityPrec<MA> >
{
    typedef IdentityPrec<MA>  Impl;
    typedef double    ElementType;
};

/*! 
 ***************************************************************************
 *  Identity Preconditioner.
 *  \author  Luis M. de la Cruz [ Wed Nov 28 10:06:43 GMT 2007 ]
 ***************************************************************************
 */  
template <typename MA>
class IdentityPrec : public SymmetricMatrix<IdentityPrec<MA> >
{
public:
    IdentityPrec(const MA &_A)  { }
};
/*! 
 ***************************************************************************
 *  Multiplication definition for the Identity preconditioner.
 *  Here the multiplaction \c alpha * \c J * \c b + \c beta * \c x is defined.
 *  In this case \c alpha = 1 and \c beta = 0, so we do \c x = \c J * \c b, but
 *  because this operation is intended for the Identity preconditioner, we
 *  finally do \c x = \c b.
 *  \param alpha scalar.
 *  \param J Identity preconditioner.
 *  \param b vector.
 *  \param beta scalar.
 *  \param x vector where the result is stored.
 *  \author  Luis M. de la Cruz [ Wed Nov 28 10:06:43 GMT 2007 ]
 ***************************************************************************
 */  
template<class Mat>
void
mv(double alpha, 
   const IdentityPrec<Mat> &J, 
   const DenseVector<Array<double> > &b,
   double beta, 
   DenseVector<Array<double> > &x)
{
    assert(alpha==1.);
    assert(beta==0.);
    x = b;
}

} // Namespace FLENS

#endif

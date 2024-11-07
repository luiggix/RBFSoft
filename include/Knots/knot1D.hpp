
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

#ifndef _KNOTS1D_HPP_
#define _KNOTS1D_HPP_

#include <fstream>
#include <vector>
#include "Traits.hpp"
#include "kdtree++/kdtree.hpp"

/*!
 ***************************************************************************
 *  The leaves of the KDTree in 1D.
 *  Objects of this structure are used by kdtree++ as the leaves of the 
 *  complete KDTree. It requires some particular behavior, explained in 
 *  http://libkdtree.alioth.debian.org/
 *  \author Luis M. de la Cruz [ Mon Dec 17 10:30:21 GMT 2007 ]
 ***************************************************************************
 */
struct knot1D 
{    
    prec_t coord[1]; ///< (x,y) knots coordinates in 1D.
    int index; ///< index inside the global Matrix, used for the ACBFprecond.
/*!
 *  Overloaded [ ] operator required by kdtree++ to access the individual
 *  dimensional components.
 *  \param N the dimensional component.
 */
    inline prec_t operator[](size_t const N) const { return coord[N]; }
};
/*!
 *  Comparison, needed by kdtree++
 */
inline bool operator==(knot1D const& A, knot1D const& B) 
{
    return A.coord[0] == B.coord[0];
}
/*!
 *  Simple output operator overloading (<<)
 */ 
std::ostream& operator<<(std::ostream& out, knot1D const& T)
{
    return out << T.coord[0] << "\n";
}
/*!
 *  Function needed by kdtree++.
 */
prec_t tac( knot1D t, int k ) { return t[k]; }
/*!
 *  KDTree type definition in 1D.
 */
typedef 
KDTree::KDTree<1, 
	       knot1D, 
	       std::pointer_to_binary_function<knot1D, 
					       int, 
					       prec_t> >  tree1D_t;
//
//  Declaration of a tree and the array of nearest neigbohrs, 
//  these variables will be used in several functions of Knots class,
//  specially in Knots::kdtree() and  Knots::findNeighborsAll().
//
tree1D_t tree(std::ptr_fun(tac)); // Tree to store all knots in 1D. 
std::vector<knot1D> neighbors_array; // Array of neighbors.
std::vector<knot1D>::const_iterator neighbors_i; // Iterator.

#endif // _KNOTS1D_HPP_

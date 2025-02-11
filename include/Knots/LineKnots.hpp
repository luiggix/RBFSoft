
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

#ifndef _LINEKNOTS_HPP_
#define _LINEKNOTS_HPP_

#include <fstream>
#include "Traits.hpp"
#include "Knots.hpp"
#include "knot1D.hpp"

/*!
 ***************************************************************************
 *  Knots in a line.
 *  This class is aimed to generate knots along the x-axis.
 *  The function constructKnots() must be reimplemented here to generate the 
 *  knots.
 *  \author  Luis M. de la Cruz [ Mon Nov 12 14:53:34 GMT 2007 ]
 ***************************************************************************
 */
template < typename Tprec>
class LineKnots : public Knots< LineKnots<Tprec>, knot1D > 
{ 
// Precision type (Tprec is for template, prec_t is cool in emacs ;) )
    typedef Tprec prec_t;

// These variables are from the Knots base class and must be defined
// in each derived class.
    using Knots<LineKnots<prec_t>, knot1D >::N;   // Total points
    using Knots<LineKnots<prec_t>, knot1D >::NI;  // Interior points
    using Knots<LineKnots<prec_t>, knot1D >::NB;  // Boundary points
    using Knots<LineKnots<prec_t>, knot1D >::Dim; // Dimension
    using Knots<LineKnots<prec_t>, knot1D >::xyz; // Coordinates
    using Knots<LineKnots<prec_t>, knot1D >::max_upper_range; // used in
                                                              // KDTREE

private:
    prec_t lx;  ///< Horizontal length of the line
    int Nx;     ///< Number of points on x-axis

public:

/*!
 *  Default constructor.
 */
    LineKnots() : Knots< LineKnots<Tprec>, knot1D > () { }

/*!
 *  Constructor to generate the knots.
 *  \param l_x Length of the line.
 *  \param n_x Number of points in the line.
 */
    LineKnots(prec_t l_x, int n_x) : Knots< LineKnots<Tprec>, knot1D > (),  
	lx(l_x), Nx (n_x) 
	{
	    Dim = 1;            // A line is defined in 1D. 
	    NB = 2;
	    NI = Nx - 2;
	    N = Nx;
	    xyz = new knot1D[N];
	    
// Set the max uppper range for the binary search of the RANGE in kdtree
	    max_upper_range = lx * 0.5;
	} 
    
    ~LineKnots() { }
/*!
 *  Construct all the knots of the domain.
 */
    bool constructKnots();
/*!
 *  Add special points for the ACBF aproximation.
 */
    inline bool addSpecialPoints(iMat &, int, int);
/*!
 *  Print info to the standard output (screen).
 */ 
    bool print();
};

template<typename prec_t> 
bool LineKnots<prec_t>::print() 
{
    Knots<LineKnots, knot1D >::print();
    std::cout << "\n +-----[ Line information ]-----+ " 
	      << "\n | RBF : Length in x axis = " << lx
	      << "\n +-----+";
    return 0;
}

/*! 
 *  Construct regularly spaced knots on a line of length lx;
 */
template<typename prec_t>
bool LineKnots<prec_t>::constructKnots() 
{
    prec_t dx = lx / (Nx - 1);

    int index = 0;
// ----- Interior points.
    for (int i = 1; i < Nx; ++i, ++index) 
	xyz[index].coord[0] = i * dx; 

// ----- Boundary points
    xyz[Nx-2].coord[0] = 0;
    xyz[Nx-1].coord[0] = lx;

    return 0;
}

/*!
 *  Add special points located in the boundary of the line.
 *  \param P  indices of neighbors knots.
 *  \param j  actual knot.
 *  \param nsp number of special points.
 */
//NOTE: Remember that the indices in the matrix P begin in 1, and the
//      indices in xyz[] begin in 0. 
//NOTE: Now the indices in the matrix P starts at 0 for EIGEN.
template<typename prec_t>
inline bool LineKnots<prec_t>::addSpecialPoints(iMat &P, int j, int nsp)
{
    int sp[2], side_x = Nx - 2;
// Corners
    sp[0] = N-2;  // 1
    sp[1] = N-1;  // 2

    /*** THIS WORKS FOR FLENS
// Add 1 to sp[] to be in agreement with P() matrix indices.
    for(int i = 0; i < 2; ++i) sp[i] += 1;
    int bi = P.firstRow();
    int ei = P.lastRow();
    ***/

    /*** THIS WORKS NOW FOR EIGEN ***/
    int bi = 0;
    int ei = P.rows();

    if ( nsp ) {
	
	for(int i = bi; i <= ei; ++i) {
	    if ( P(i,j) == sp[0] ) sp[0] = -1; 
	    if ( P(i,j) == sp[1] ) sp[1] = -1;
	}
	
	ei -= 1;

	if (sp[0] != -1) P(ei,j) = sp[0];
	else P(ei,j) = N - 2;
	
	if (sp[1] != -1) P(ei+1,j) = sp[1];
	else P(ei+1,j) = 2;
    }

    return 0;
}

#endif // _LINEKNOTS_HPP_

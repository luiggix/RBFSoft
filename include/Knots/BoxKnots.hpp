
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

#ifndef _BOXKNOTS_HPP_
#define _BOXKNOTS_HPP_

#include <fstream>
#include "Traits.hpp"
#include "Knots.hpp"
#include "knot3D.hpp"

/*!
 ***************************************************************************
 *  Octoedro or Cuboide.
 *  This class is aimed to generate knots inside of 3D rectangular domains. 
 *  The function constructKnots() must be reimplemented here to generate the 
 *  knots on the rectangular domain.
 *  \author  Luis M. de la Cruz [ Mon Nov 12 14:53:34 GMT 2007 ]
 ***************************************************************************
 */
template < typename Tprec>
class BoxKnots : public Knots< BoxKnots<Tprec>, knot3D > 
{ 
// Precision type (Tprec is for template, prec_t is cool in emacs ;) )
    typedef Tprec prec_t;

// These variables are from the Knots base class and must be defined
// in each derived class.
    using Knots<BoxKnots<prec_t>, knot3D >::N;   // Total points
    using Knots<BoxKnots<prec_t>, knot3D >::NI;  // Interior points
    using Knots<BoxKnots<prec_t>, knot3D >::NB;  // Boundary points
    using Knots<BoxKnots<prec_t>, knot3D >::Dim; // Dimension
    using Knots<BoxKnots<prec_t>, knot3D >::xyz; // Coordinates
    using Knots<BoxKnots<prec_t>, knot3D >::max_upper_range; // used for
                                                             // KDTREE

private:
    prec_t hx;  ///< Horizontal length of the box
    prec_t hy;  ///< Vertical length of the box
    prec_t hz;  ///< Depth length of the box
    int Nx;     ///< Number of points on x-axis
    int Ny;     ///< Number of points on y-axis
    int Nz;     ///< Number of points on z-axis
    prec_t dx, dy, dz;
    prec_t ep;  ///< Randomness of the points distribution (default ep = 0)
    int cel;    ///< Used to keep uniform knots in one cell near the boundary (default cel = 0)

public:

/*!
 *  Default constructor.
 */
    BoxKnots() : Knots< BoxKnots<Tprec>, knot3D > () { }
/*!
 *  Constructor to generate knots in the centers of the "cells" of
 *  a gridded mesh. It allows some randomness of the points.
 *  \param l_x Length of the box in x-axis direction.
 *  \param l_y Length of the box in y-axis direction.
 *  \param l_z Length of the box in z-axis direction.
 *  \param n_x Number of points in x-axis direction.
 *  \param n_y Number of points in y-axis direction.
 *  \param n_z Number of points in z-axis direction.
 */
    BoxKnots(prec_t l_x, int n_x, prec_t l_y, int n_y, prec_t l_z, int n_z)
      : Knots< BoxKnots<Tprec>, knot3D > (),  
	hx(l_x), Nx (n_x), hy(l_y), Ny(n_y), hz(l_z), Nz(n_z) {

        Dim = 3;            // A box is defined in 3D.
//	ep = (1 - e) * 0.5; // e = 1 total randomness, 0 no randomness
 
	NB = 2 * (Nx * Ny + Ny * Nz + Nx * Nz) - 4 * (Nx + Ny + Nz) + 8;
	NI = (Nx - 2) * (Ny - 2) * (Nz - 2);
	N = NI + NB;
	xyz = new knot3D[N];
	
// Set the max uppper range for the binary search of the RANGE in kdtree
	max_upper_range = (hy > hx) ? hx : hy;
	max_upper_range = (max_upper_range > hz) ? hz : max_upper_range;
    } 
    ~BoxKnots() { }

/*!
 *  Return the size of the mesh in x-axis.
 */
    prec_t getDx() { return dx; }
/*!
 *  Return the size of the mesh in y-axis.
 */
    prec_t getDy() { return dy; }
/*!
 *  Return the size of the mesh in z-axis.
 */
    prec_t getDz() { return dz; }
/*!
 *  Return the length of the box in x-axis.
 */
    prec_t getLx() { return hx; }
/*!
 *  Return the length of the box in y-axis.
 */
    prec_t getLy() { return hy; }
/*!
 *  Return the length of the box in z-axis.
 */
    prec_t getLz() { return hz; }
/*!
 *  Return the number of points in x-axis.
 */
    prec_t getNx() { return Nx; }
/*!
 *  Return the number of points in y-axis.
 */
    prec_t getNy() { return Ny; }
/*!
 *  Return the number of points in z-axis.
 */
    prec_t getNz() { return Nz; }
/*!
 *  Construct all the knots of the domain.
 */
    bool constructKnots();
/*!
 *  Add special points for the ACBF aproximation (is not defined yet).
 */
    inline bool addSpecialPoints(iMat &, int, int);
/*!
 *  Print info to the standard output (screen).
 */ 
    bool print();
};

template<typename prec_t> 
bool BoxKnots<prec_t>::print() 
{
    std::string flag = "true";
    if (cel == 0) flag = "false";
    Knots<BoxKnots, knot3D >::print();
    std::cout << "\n +-----[ Box information ]-----+ " 
	      << "\n | RBF : Length in x axis = " << hx
	      << "\n | RBF : Length in y axis = " << hy
	      << "\n | RBF : Length in z axis = " << hz
	      << "\n | RBF : Boundary points in x-axis = " << Nx
	      << "\n | RBF : Boundary points in y-axis = " << Ny
	      << "\n | RBF : Boundary points in z-axis = " << Nz
//	      << "\n | RBF : Randomness (epsilon [0 to 1]) = " << 1 - 2 * ep
//	      << "\n | RBF : Uniform cell adjacent to boundary = " << flag
	      << "\n +-----+";
    return 0;
}

/*!
 *  Detailed documentation is missing, it is coming soon.
 *
 */
template<typename prec_t>
bool BoxKnots<prec_t>::constructKnots() 
{
//    if (ep < 0.5) srand( time(0) );  
    
    dx = hx / (Nx - 1);
    dy = hy / (Ny - 1);
    dz = hz / (Nz - 1);

// ----- Interior points.
    int index = 0; // The xyz is a C-array, so begins in 0.
    for (int k = 1; k < Nz - 1; ++k)   
	for (int j = 1; j < Ny - 1; ++j)
	    for (int i = 1; i < Nx - 1; ++i, ++index) {   
		xyz[index].coord[0] = i * dx;
		xyz[index].coord[1] = j * dy;
		xyz[index].coord[2] = k * dz;
	    }

// ----- Boundary points

// Faces 1 and 2: x = constant
   for (int i = 0; i <= 1; ++i)
	for (int k = 1; k < Nz-1; ++k) 
	    for (int j = 1; j < Ny-1; ++j, ++index) {	    
		xyz[index].coord[0] = hx * i;
		xyz[index].coord[1] = j * dy;
		xyz[index].coord[2] = k * dz;
	    }

// Faces 3 and 4: y = constant
   for (int j = 0; j <= 1; ++j)
	for (int k = 1; k < Nz-1; ++k) 
	    for (int i = 1; i < Nx-1; ++i, ++index) {	    
		xyz[index].coord[0] = i * dx;
		xyz[index].coord[1] = hy * j;
		xyz[index].coord[2] = k * dz;
	    }

// Faces 5 and 6: z = constant
    for (int k = 0; k <= 1; ++k)
	for (int j = 1; j < Ny-1; ++j) 
	    for (int i = 1; i < Nx-1; ++i, ++index) {
		xyz[index].coord[0] = i * dx;
		xyz[index].coord[1] = j * dy;
		xyz[index].coord[2] = hz * k;
	    }

// Edges 1, 2, 3 and 4: x = constant, y = constant
    for(int i = 0; i <= 1; ++i) 
	for(int j = 0; j <= 1; ++j) 		   
	    for(int k = 1; k < Nz - 1; ++k, ++index) {
		xyz[index].coord[0] = hx * i;
		xyz[index].coord[1] = hy * j;
		xyz[index].coord[2] = k * dz;	
	    }


// Edges 5, 6, 7 and 8: x = constant, z = constant
    for(int i = 0; i <= 1; ++i) 
	for(int k = 0; k <= 1; ++k) 
	    for(int j = 1; j < Ny-1; ++j, ++index) {
		xyz[index].coord[0] = hx * i;
		xyz[index].coord[1] = j * dy;
		xyz[index].coord[2] = hz * k;	
	    }

// Edges 9, 10, 11 and 12: y = constant, z = constant
    for(int j = 0; j <= 1; ++j) 
	for(int k = 0; k <= 1; ++k) 	    
	    for(int i = 1; i < Nx-1; ++i, ++index) {
		xyz[index].coord[0] = i * dx;
		xyz[index].coord[1] = hy * j;
		xyz[index].coord[2] = hz * k;	
	    }

// Corners of the box
    for(int i = 0; i <= 1; ++i) 
	for(int j = 0; j <= 1; ++j) 
	    for(int k = 0; k <= 1; ++k, ++index) { 
		xyz[index].coord[0] = i * hx;
		xyz[index].coord[1] = j * hy;
		xyz[index].coord[2] = k * hz;
	    }

    return 0;
}



#endif // _BOXKNOTS_HPP_

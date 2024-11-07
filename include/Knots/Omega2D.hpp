
/*------------------------------------------------------------------------
 *  Copyright (C) 2015  Luis M. de la Cruz
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

#ifndef _OMEGA2D_HPP_
#define _OMEGA2D_HPP_

#include <fstream>
#include "Traits.hpp"
#include "Knots.hpp"
#include "knot2D.hpp"

/*!
 ***************************************************************************
 *  Generic two dimensional domain.
 *  This class is aimed to generate knots inside of two dimensional domains. 
 *  The function constructKnots() must be reimplemented here to generate the 
 *  knots on the domain.
 *  \author  Luis M. de la Cruz [ Sun May  3 19:05:31 CDT 2015 ]
 ***************************************************************************
 */
template < typename Tprec>
class Omega2D : public Knots< Omega2D<Tprec>, knot2D > 
{ 
// Precision type (Tprec is for template, prec_t is cool in emacs ;) )
    typedef Tprec prec_t;

// These variables are from the Knots base class and must be defined
// in each derived class.
    using Knots<Omega2D<prec_t>, knot2D >::N;   // Total points
    using Knots<Omega2D<prec_t>, knot2D >::NI;  // Interior points
    using Knots<Omega2D<prec_t>, knot2D >::NB;  // Boundary points
    using Knots<Omega2D<prec_t>, knot2D >::Dim; // Dimension
    using Knots<Omega2D<prec_t>, knot2D >::xyz; // Coordinates
    using Knots<Omega2D<prec_t>, knot2D >::max_upper_range; // used for
                                                            // KDTREE

private:
    prec_t lx;   ///< Horizontal length of the rectangle
    prec_t ly;   ///< Vertical length of the rectangle
    int Nx;      ///< Number of points on x-axis
    int Ny;      ///< Number of points on y-axis
    random_t knots_dist;   
    prec_t epru; ///< Used to control randomness in regular distr.
    int NG;      ///< Number of ghost points around de domain.
    bool corners;
    
public:

/*!
 *  Default constructor.
 */
    Omega2D() : Knots< Omega2D<Tprec>, knot2D > () { }
/*!
 *  Constructor to generate knots in the centers of the "cells" of
 *  a gridded mesh, totally random knots or controlled random knots. 
 *  \param l_x Length of the rectangle in x-axis direction.
 *  \param l_y Length of the rectangle in y-axis direction.
 *  \param n_x Number of points in x-axis direction.
 *  \param n_y Number of points in y-axis direction.
 *  \param knd distribution of points (0 unif, 1 rand, 2 rand unif).
 *  \param layer 0 no layers around the domain, 1 just one layer.
 */
    Omega2D(prec_t l_x, int n_x, prec_t l_y, int n_y,
		   random_t knd = U, int layer = 0)
	: Knots< Omega2D<Tprec>, knot2D > (),  
	  lx(l_x), Nx (n_x), ly(l_y), Ny(n_y), knots_dist(knd) {
	
        Dim = 2;            // A rectangle is defined in 2D. 
	NB = 2 * (Nx + Ny) - 4;

	if (knots_dist == RU)
	    NI = (Nx - 1) * (Ny - 1); // Random wrapped in cells
	else
	    NI = (Nx - 2) * (Ny - 2); // Uniform or Random
	
	NG = 0;
	if (layer > 0) { // One layer outside the boundary
	    NG = 2 * (Nx + Ny);
	    NI += NG;
	}

	N = NI + NB;
	xyz = new knot2D[N];
	
// Set the max uppper range for the binary search of the RANGE in kdtree
	max_upper_range = (ly > lx) ? lx : ly;

// Default randomnesss value:
	epru = 0.0;

	corners = true;
    } 

    ~Omega2D() { }
/*!
 *  Construct all the knots of the domain.
 */
    bool constructKnots();

/*!
 *  Get the number of Ghost knots around the domain.
 */
    int getGhostKnots() { return NG; }

/*!
 *  Add special points for the ACBF aproximation.
 */
    inline bool addSpecialPoints(iMat &, int, int);
/*!
 *  Random parameter
 */ 
    void setRandomness(prec_t e) { epru = (1 - e) * 0.5 ; }
/*!
 *  Generate knots in a "regular grid".
 */
    bool gridKnots();
/*!
 *  Print info to the standard output (screen).
 */ 
    bool print();

    void setCorners(bool c) { 
	corners = c;
	if (!corners) { NB -= 4; N -= 4; }
	
    }
};

template<typename prec_t> 
bool Omega2D<prec_t>::print() 
{
    std::string dist_type = "UNIFORM";
    if (knots_dist == R) dist_type = "RANDOM";
    if (knots_dist == RU) dist_type = "RANDOM-CELLS";

    Knots<Omega2D, knot2D >::print();
    std::cout << "\n +-----[ Rectangle information ]-----+ " 
	      << "\n | RBF : Ghost points = " << NG
	      << "\n | RBF : Length in x axis = " << lx
	      << "\n | RBF : Length in y axis = " << ly
	      << "\n | RBF : Boundary points in x-axis = " << Nx
	      << "\n | RBF : Boundary points in y-axis = " << Ny
	      << "\n | RBF : Knots distribution = " << dist_type
	      << "\n | RBF : Randomness (epsilon [0 to 1]) = " << 1 - 2 * epru
	      << "\n +-----+";
    return 0;
}

/*!
 *  This function generate nodes in a uniform grid (knd = 0), randomly
 *  distributed (knd = 1) and random controlled ( knd = 2, 0 < epru < 0.5). 
 *  The random controlled is generated as follows: Taking the values of \c Nx 
 *  and \c Ny an "invisible grid" is defined inside the rectangle, this grid 
 *  will contains (\c Nx - 1) \c X (\c Ny - 1) "cells", for each cell a knot 
 *  is generated and this knot is allowed to be in a random position inside 
 *  its cell, the randomness depends on the \c ep parameter. The order of the 
 *  knots in the \c xyz array is as shown in the next figure:
 * \image html  Omega2D.png "Knots order" width=5cm 
 * \image latex Omega2D.eps "Knots order" width=5cm  
 */
template<typename prec_t>
bool Omega2D<prec_t>::constructKnots() 
{
    prec_t dx = lx / (Nx - 1);
    prec_t dy = ly / (Ny - 1);
    prec_t x0, x1, y0, y1;
    prec_t ep = dx > dy ? dy * 0.1 : dx * 0.1;
    
// ----- Interior points.
    int ei = Nx - 2, ej = Ny - 2;
    
    int index = 0; // The xyz is a C-array, so begins in 0.
    
    switch (knots_dist) 
    {
	case 0: // Uniform 
	    for (int j = 0; j < ej; ++j)
		for (int i = 0; i < ei; ++i, ++index) {
		    xyz[index].coord[0] = (i+1) * dx; 
		    xyz[index].coord[1] = (j+1) * dy;
		}
	    break;		
	case 1: // Random	
	    srand( time(0) );
	    for (int j = 0; j < ej ; ++j)
		for (int i = 0; i < ei; ++i, ++index) { 
		    xyz[index].coord[0] = ep +
			(1-2*ep) * lx * static_cast<prec_t>(rand())/RAND_MAX; 
		    xyz[index].coord[1] = ep +
			(1-2*ep) * ly * static_cast<prec_t>(rand())/RAND_MAX; 
		}
	    break;
	case 2:	// Random controlled 
	    srand( time(0) );
	    ++ei; ++ej;
	    for (int j = 0; j < ej; ++j)
		for (int i = 0; i < ei; ++i, ++index) {
		    x0 = (i + epru) * dx; x1 = (i + 1 - epru) * dx;// Mapping
		    y0 = (j + epru) * dy; y1 = (j + 1 - epru) * dy;// to a cell
		    xyz[index].coord[0] = x0 + 
			( static_cast<prec_t>(rand())/RAND_MAX ) * (x1 - x0);
		    xyz[index].coord[1] = y0 +
			( static_cast<prec_t>(rand())/RAND_MAX ) * (y1 - y0);
		}
	    break;
	default:
	    break;
    }
    
// ----- Ghost points    
    if ( NG > 0 ) {
// Generation of extra ghost-knots around the domain.
// The corners are avoided.
	
// Edge 1: x = 0 - dx
	for (int j = 0; j < Ny; ++j, ++index) {
	    xyz[index].coord[0] = -dx;
	    xyz[index].coord[1] = j * dy;
	}
// Edge 2: x = lx + dx
	for (int j = 0; j < Ny; ++j, ++index) {
	    xyz[index].coord[0] = lx + dx;
	    xyz[index].coord[1] = j * dy;
	}
	
// Edge 3: y = 0 - dy
	for (int i = 0; i < Nx; ++i, ++index) {
	    xyz[index].coord[0] = i * dx;
	    xyz[index].coord[1] = -dy;
	}
	
// Edge 4: y = ly + dy
	for (int i = 0; i < Nx; ++i, ++index) {
	    xyz[index].coord[0] = i * dx;
	    xyz[index].coord[1] = ly + dy;
	}
    }
    
// ----- Boundary points

// Edges 1 and 2: x = constant
    for (int i = 0; i <= 1; ++i)
	for (int j = 1; j < Ny-1; ++j, ++index) {
	    xyz[index].coord[0] = lx * i;
	    xyz[index].coord[1] = j * dy;
	}
    
// Edges 3 and 4: y = constant
    for (int j = 0; j <= 1; ++j)
	for (int i = 1; i < Nx-1; ++i, ++index) {
	    xyz[index].coord[0] = i * dx;
	    xyz[index].coord[1] = ly * j;
	}

// Corners    
    if (corners) {
	for (int i = 0; i <= 1; ++i)
	    for (int j = 0; j <= 1; ++j, ++index) {
		xyz[index].coord[0] = lx * i;
		xyz[index].coord[1] = ly * j;
	    }
    }
    return 0;
}

/*!
 *  Add special points located in the boundary of the rectangle.
 *  \param P  indices of neighbors knots.
 *  \param j  actual knot.
 *  \param nsp number of special points.
 */
//NOTE: Remember that the indices in the matrix P begin in 1, and the
//      indices in xyz[] begin in 0. 
template<typename prec_t>
inline bool Omega2D<prec_t>::addSpecialPoints(iMat &P, int j, int nsp)
{
    int sp[9], side_x = Nx - 2, side_y = Ny - 2;
// Corners
    sp[0] = N - 4;      // 1
    sp[1] = sp[0] + 1;  // 2
    sp[2] = sp[1] + 1;  // 3
    sp[3] = sp[2] + 1;  // 4

// Middle Edges
    sp[4] = NI + side_y / 2;  // Edge 1 
    sp[5] = sp[4] + side_y;   // Edge 2
    sp[6] = NI + 2 * side_y + side_x / 2;  // Edge 3 
    sp[7] = sp[6] + side_x;                // Edge 4

/* *
    std::cout << std::endl;
    for(int i = 0; i < 8; ++i) 
	std::cout << sp[i] << "\t"
		  << xyz[ sp[i] ].coord[0] << "\t"
		  << xyz[ sp[i] ].coord[1] << "\n";
    std::cout << std::endl;
/* */

/* *
    std::ofstream sfile ("special.dat");
    for(int i = 0; i < 8; ++i)
	sfile << xyz[ sp[i] ].coord[0] << "\t"
	      << xyz[ sp[i] ].coord[1] << "\n";
    sfile.close();
/* */

// Add 1 to sp[] to be in agreement with P() matrix indices.

/**** FOR EIGEN THIS IS NOT NEEDED
    for(int i = 0; i < 8; ++i) sp[i] += 1;
****/
/* *
    std::cout << std::endl;
    for(int i = 0; i < 8; ++i) 
	std::cout << sp[i] << "\t"
		  << xyz[ sp[i] ].coord[0] << "\t"
		  << xyz[ sp[i] ].coord[1] << "\n";
    std::cout << std::endl;
/* */


    int bi = 0; // P.firstRow();
    int ei = P.rows(); // P.lastRow();
    
    if ( nsp ) {

	for(int i = bi; i <= ei; ++i) {
	    if ( P(i,j) == sp[0] ) sp[0] = -1; 
	    if ( P(i,j) == sp[1] ) sp[1] = -1;
	    if ( P(i,j) == sp[2] ) sp[2] = -1;
	    if ( P(i,j) == sp[3] ) sp[3] = -1;
	    
	    if (nsp > 4) {
		if ( P(i,j) == sp[4] ) sp[4] = -1; 
		if ( P(i,j) == sp[5] ) sp[5] = -1;
		if ( P(i,j) == sp[6] ) sp[6] = -1;
		if ( P(i,j) == sp[7] ) sp[7] = -1;
	    }
	}
	
	if ( nsp == 4) ei -= 3;
	else ei -= 7;

	if (sp[0] != -1) P(ei,j) = sp[0];
	else P(ei,j) = sp[5] + side_y / 4;
	
	if (sp[1] != -1) P(ei+1,j) = sp[1];
	else P(ei+1,j) = sp[6] + side_x / 4;
	
	if (sp[2] != -1) P(ei+2,j) = sp[2];
	else P(ei+2,j) = sp[7] - side_x / 4;
	
	if (sp[3] != -1) P(ei+3,j) = sp[3];
	else P(ei+3,j) = sp[4] - side_y / 4;
	
	if (nsp > 4) {       
	    if (sp[4] != -1) P(ei+4,j) = sp[4];
	    else P(ei+4,j) = sp[5] - side_y / 4;
	    
	    if (sp[5] != -1) P(ei+5,j) = sp[5];
	    else P(ei+5,j) = sp[4] + side_y / 4;
	    
	    if (sp[6] != -1) P(ei+6,j) = sp[6];
	    else P(ei+6,j) = sp[7] + side_x / 4;
	    
	    if (sp[7] != -1) P(ei+7,j) = sp[7];
	    else P(ei+7,j) = sp[6] - side_x / 4;
	}
    }

/* *
    ei = P.lastRow();
    if ( nsp == 4) ei -= 3;
    else ei -= 7;

    std::ofstream sfile ("special.dat");
    for(int i = 0; i <= P.lastRow() - ei; ++i)
	sfile << xyz[ P(ei + i, j) - 1 ].coord[0] << "\t"
	      << xyz[ P(ei + i, j) - 1 ].coord[1] << "\n";
    sfile.close();
/* */

    return 0;
}

/*!
 *  This function is aimed to generate the knots in a Cartesian 
 *  uniform grid.
 */
template<typename prec_t>
bool Omega2D<prec_t>::gridKnots() 
{
    prec_t dx = lx / (Nx - 1);
    prec_t dy = ly / (Ny - 1);
    int index = 0;
    for (int j = 0; j < Ny; ++j)
	for (int i = 0; i < Nx; ++i, ++index) {   
	    xyz[index].coord[0] = i * dx; 
	    xyz[index].coord[1] = j * dy;
	}
    
    return 0;
}	




#endif // _OMEGA2D_HPP_

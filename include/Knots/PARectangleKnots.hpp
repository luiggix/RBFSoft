
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

#ifndef _PARECTANGLEKNOTS_HPP_
#define _PARECTANGLEKNOTS_HPP_

#include <fstream>
#include "Traits.hpp"
#include "DD/CartComm.hpp"
#include "Knots.hpp"
#include "knot2D.hpp"

/*!
 ***************************************************************************
 *  PARectangular domain (NO detailed documentation).
 *  This class is aimed to generate knots inside of rectangular domains. 
 *  The function constructKnots() must be reimplemented here to generate the 
 *  knots on the rectangular domain.
 *  \author  Luis M. de la Cruz [ Mon Nov 12 14:53:34 GMT 2007 ]
 ***************************************************************************
 */
template < typename Tprec>
class PARectangleKnots : public Knots< PARectangleKnots<Tprec>, knot2D > 
{ 
// Precision type (Tprec is for template, prec_t is cool in emacs ;) )
    typedef Tprec prec_t;

// These variables are from the Knots base class and must be defined
// in each derived class.
    using Knots<PARectangleKnots<prec_t>, knot2D >::N;   // Total points
    using Knots<PARectangleKnots<prec_t>, knot2D >::NI;  // Interior points
    using Knots<PARectangleKnots<prec_t>, knot2D >::NB;  // Boundary points
    using Knots<PARectangleKnots<prec_t>, knot2D >::Dim; // Dimension
    using Knots<PARectangleKnots<prec_t>, knot2D >::xyz; // Coordinates
    using Knots<PARectangleKnots<prec_t>, knot2D >::max_upper_range; // used for
                                                                   // KDTREE

private:
    prec_t lx;  ///< Horizontal length of the rectangle
    prec_t ly;  ///< Vertical length of the rectangle
    int Nx;     ///< Number of points on x-axis
    int Ny;     ///< Number of points on y-axis
    random_t knots_dist;   
    prec_t epru;  ///< Randomness of the points distribution (default ep = 0)
    int NG;
    int rank; // For parallel implementation, id of the procesor
    prec_t shift[2];
    prec_t dx, dy;

public:

/*!
 *  Default constructor.
 */
    PARectangleKnots() : Knots< PARectangleKnots<Tprec>, knot2D > () { }
/*!
 *  Constructor to generate knots in the centers of the "cells" of
 *  a gridded mesh. It allows some randomness of the points.
 *  \param l_x Length of the rectangle in x-axis direction.
 *  \param l_y Length of the rectangle in y-axis direction.
 *  \param n_x Number of points in x-axis direction.
 *  \param n_y Number of points in y-axis direction.
 *  \param e Degree of randomness of the points distribution (default ep = 0).
 */
    PARectangleKnots(prec_t l_x, int n_x, prec_t l_y, int n_y, 
		     CartComm& cart, random_t knd = U) 
	: Knots< PARectangleKnots<Tprec>, knot2D > (), knots_dist(knd) {

// ----- Get the correct lenght, num. of points and shift's 
// ----- every subdomain, based on its position in the CartComm
// The numeration of CartComm is as in matrices of Linear Algebra
// therefore the I goes in y direction, while J goes in x direction
// (I,J) \in [0, NP_I - 1] X [0, NP_J - 1]
//
	rank = cart.comm.Get_rank();
	int NP_I = cart.getNumProc_I();
	int NP_J = cart.getNumProc_J();
	int I = cart.get_I();
	int J = cart.get_J();
	if ( (n_x - 2) % NP_J ) {
	    std::cerr << "\n\n DD::CartComm : Error : num of nodes minus 2 (" 
		      << n_x - 2 
		      << ") is not divisible by num of procesors (" 
		      << NP_J
		      << ") \n\n";
	    std::exit(0);
	}
	if ( (n_y - 2) % NP_I ) {
	    std::cerr << "\n\n DD::CartComm : Error : num of nodes minus 2 (" 
		      << n_y - 2 
		      << ") is not divisible by num of procesors (" 
		      << NP_I
		      << ") \n\n";
	    std::exit(0);
	}

	dx = l_x / (n_x - 1);
	dy = l_y / (n_y - 1);
	Nx = (n_x - 2) / NP_J;
	Ny = (n_y - 2) / NP_I;
	lx = (l_x - (NP_J + 1) * dx) / NP_J;
	ly = (l_y - (NP_I + 1) * dy) / NP_I;
	shift[0] = (lx + dx) * J + dx; 
	shift[1] = (ly + dy) * (NP_I - I - 1) + dy;

// ----- Now define some parameters to construct the rectangle

        Dim = 2;            // A rectangle is defined in 2D. 
	NB = 2 * (Nx + Ny) - 4;

	if (knots_dist == RU)
	    NI = (Nx - 1) * (Ny - 1); // Random wrapped in cells
	else
	    NI = (Nx - 2) * (Ny - 2); // Uniform or Random

	NG = NB + 4 + 4; // I'm using the 4 corners 

	N = NI + NB + NG; 
	xyz = new knot2D[N];
	
// Set the max uppper range for the binary search of the RANGE in kdtree
	max_upper_range = (ly > lx) ? lx : ly;

// Default randomnesss value:
	epru = 0.1;
    } 

    ~PARectangleKnots() { }
/*!
 *  Construct all the knots of the domain.
 */
    bool constructKnots();

    bool constructKnotsDriven();

/*!
 *  Add special points for the ACBF aproximation.
 */
    inline bool addSpecialPoints(iMat &, int, int);
/*!
 *  Random parameter
 */ 
    inline void setRandomness(prec_t e) { epru = (1 - e) * 0.5 ; }
    inline prec_t getShapeParameter() { return 1.0 / sqrt(N); }
    inline int getGhostKnots() { return NG; }
    inline int getNx() { return Nx; }
    inline int getNy() { return Ny; }
    inline prec_t getLx() { return lx; }
    inline prec_t getLy() { return ly; }
    inline prec_t getDx() { return dx; }
    inline prec_t getDy() { return dy; }    
    inline prec_t getShift(axis_t a) { return shift[a]; }
/*!
 *  Generate knots in a "Cartesian grid".
 */
    bool gridKnots();
/*!
 *  Print info to the standard output (screen).
 */ 
    bool print();
};

template<typename prec_t> 
bool PARectangleKnots<prec_t>::print() 
{
    std::string dist_type = "UNIFORM";
    if (knots_dist == R) dist_type = "RANDOM";
    if (knots_dist == RU) dist_type = "RANDOM-CELLS";

    Knots<PARectangleKnots, knot2D >::print();
    std::cout << "\n +-----[ PARectangle information ]-----+ " 
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
 *  The way the knots are generated is tricky: Taking the values of \c Nx and
 *  \c Ny an "invisible grid" is defined inside the rectangle, this grid will
 *  contains (\c Nx - 1) \c X (\c Ny - 1) "cells", for each cell a knot is 
 *  generated and this knot is allowed to be in a random position inside its 
 *  cell, the randomness depends on the \c ep parameter. The order of the 
 *  knots are as shown in the next figure:
 * \image html  PARectangleKnots.png "Knots order" width=5cm 
 * \image latex PARectangleKnots.eps "Knots order" width=5cm  
 */
template<typename prec_t>
bool PARectangleKnots<prec_t>::constructKnots() 
{
/*
    prec_t dx = lx / (Nx - 1);
    prec_t dy = ly / (Ny - 1);
*/
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
		    xyz[index].coord[0] = (i+1) * dx + shift[0]; 
		    xyz[index].coord[1] = (j+1) * dy + shift[1];
		}
	    break;		
	case 1: // Random	
	    srand( time(0) + rank);
	    for (int j = 0; j < ej ; ++j)
		for (int i = 0; i < ei; ++i, ++index) { 
		    xyz[index].coord[0] = ep + shift[0] +
			(1-2*ep) * lx * static_cast<prec_t>(rand())/RAND_MAX; 
		    xyz[index].coord[1] = ep + shift[1] +
			(1-2*ep) * ly * static_cast<prec_t>(rand())/RAND_MAX; 
		}
	    break;
	case 2:	  
	    srand( time(0) + rank);
	    ++ei; ++ej;
	    for (int j = 0; j < ej; ++j)
		for (int i = 0; i < ei; ++i, ++index) {
		    x0 = (i + epru) * dx; x1 = (i + 1 - epru) * dx;// Mapping
		    y0 = (j + epru) * dy; y1 = (j + 1 - epru) * dy;// to a cell
		    xyz[index].coord[0] = x0 + shift[0] +
			( static_cast<prec_t>(rand())/RAND_MAX ) * (x1 - x0);
		    xyz[index].coord[1] = y0 + shift[1] +
			( static_cast<prec_t>(rand())/RAND_MAX ) * (y1 - y0);
		}
	    break;
	default:
	    break;
    }

// ----- "Internal" Ghost Boundaries to exchange information

// Edges 1: y = shift[1]
    
    for (int i = 0; i < Nx - 1; ++i, ++index) {
	xyz[index].coord[0] = i * dx + shift[0];
	xyz[index].coord[1] = shift[1];
    }

// Edges 2: x = lx + shift[0]
    for (int j = 0; j < Ny - 1; ++j, ++index) {
	xyz[index].coord[0] = lx + shift[0];
	xyz[index].coord[1] = j * dy + shift[1];
    }

// Edges 3: y = ly + shift[1]
    
    for (int i = Nx - 1; i >= 1; --i, ++index) {
	xyz[index].coord[0] = i * dx + shift[0];
	xyz[index].coord[1] = ly + shift[1];
    }

// Edges 4: x = shift[0]
    for (int j = Ny - 1; j >= 1; --j, ++index) {
	xyz[index].coord[0] = shift[0];
	xyz[index].coord[1] = j * dy + shift[1];
    }

// ----- Ghost Boundaries that acts as boundaries

// Edges 1: x = shift[0] - dx
    for (int j = 0; j < Ny; ++j, ++index) {
	xyz[index].coord[0] = shift[0] - dx;
	xyz[index].coord[1] = shift[1] + j * dy;
    }

// Edges 2: y = shift[1] + ly + dy
    for (int i = 0; i < Nx; ++i, ++index) {
	xyz[index].coord[0] = shift[0] + i * dx;
	xyz[index].coord[1] = shift[1] + ly + dy;
    }

// Edges 3: x = shift[0] + lx + dx
    for (int j = Ny - 1; j >= 0; --j, ++index) {
	xyz[index].coord[0] = shift[0] + lx + dx;
	xyz[index].coord[1] = shift[1] + j * dy;
    }

// Edges 4: y = shift[1] - dy
    for (int i = Nx - 1; i >= 0; --i, ++index) {
	xyz[index].coord[0] = shift[0] + i * dx;
	xyz[index].coord[1] = shift[1] - dy;
    }


// Corners:
// 1
    xyz[N-4].coord[0] = shift[0] - dx;
    xyz[N-4].coord[1] = shift[1] - dy;
// 2
    xyz[N-3].coord[0] = shift[0] - dx;
    xyz[N-3].coord[1] = shift[1] + ly + dy;

// 3
    xyz[N-2].coord[0] = shift[0] + lx + dx;
    xyz[N-2].coord[1] = shift[1] + ly + dy;
// 4
    xyz[N-1].coord[0] = shift[0] + lx + dx;
    xyz[N-1].coord[1] = shift[1] - dy;

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
inline bool PARectangleKnots<prec_t>::addSpecialPoints(iMat &P, int j, int nsp)
{
    int sp[9], side_x = Nx / 2, side_y = Ny / 2;
// Corners
    sp[0] = N - 4;      // 1
    sp[1] = sp[0] + 1;  // 2
    sp[2] = sp[1] + 1;  // 3
    sp[3] = sp[2] + 1;  // 4

    int sx_sy = side_x + side_y;
// Middle Edges
    sp[4] = NI + NB + side_y;  // Edge 1 
    sp[5] = sp[4] + sx_sy;     // Edge 2
    sp[6] = sp[5] + sx_sy;     // Edge 3 
    sp[7] = sp[6] + sx_sy;     // Edge 4

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
    for(int i = 0; i < 8; ++i) sp[i] += 1;

/* *
    std::cout << std::endl;
    for(int i = 0; i < 8; ++i) 
	std::cout << sp[i] << "\t"
		  << xyz[ sp[i] ].coord[0] << "\t"
		  << xyz[ sp[i] ].coord[1] << "\n";
    std::cout << std::endl;
/* */


    int bi = P.firstRow();
    int ei = P.lastRow();
    
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
	else P(ei,j) = NI + Nx + Ny;
	
	if (sp[1] != -1) P(ei+1,j) = sp[1];
	else P(ei+1,j) = NI + Nx;
	
	if (sp[2] != -1) P(ei+2,j) = sp[2];
	else P(ei+2,j) = 1;
	
	if (sp[3] != -1) P(ei+3,j) = sp[3];
	else P(ei+3,j) = NI + 2 * Nx + Ny;
	
	if (nsp > 4) {       
	    if (sp[4] != -1) P(ei+4,j) = sp[4];
	    else P(ei+4,j) = NI + Nx + side_y;
	    
	    if (sp[5] != -1) P(ei+5,j) = sp[5];
	    else P(ei+5,j) = NI + side_x;
	    
	    if (sp[6] != -1) P(ei+6,j) = sp[6];
	    else P(ei+6,j) = NI + 2 * Nx + Ny + side_x;
	    
	    if (sp[7] != -1) P(ei+7,j) = sp[7];
	    else P(ei+7,j) = NI + Nx + Ny + side_x;
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
bool PARectangleKnots<prec_t>::gridKnots() 
{
/*
    prec_t dx = lx / (Nx - 1);
    prec_t dy = ly / (Ny - 1);
*/
    int index = 0;
    for (int j = 0; j < Ny; ++j)
	for (int i = 0; i < Nx; ++i, ++index) {   
	    xyz[index].coord[0] = i * dx + shift[0]; 
	    xyz[index].coord[1] = j * dy + shift[1];
	}
    
    return 0;
}	


/* */

template<typename prec_t>
bool PARectangleKnots<prec_t>::constructKnotsDriven() 
{

    NG = (2 * Nx + 2 * Ny - 8);
    NB -= 4;
    N = NI + NB + NG;

    delete [] xyz;
    xyz = new knot2D[N];

/*
    prec_t dx = lx / (Nx - 1);
    prec_t dy = ly / (Ny - 1);
*/
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
	case 2:	  
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
// Edges 1: x = 0 - dx
    for (int j = 1; j < Ny-1; ++j, ++index) {
	xyz[index].coord[0] = -dx;
	xyz[index].coord[1] = j * dy;
    }
// Edges 2: x = lx + dx
    for (int j = 1; j < Ny-1; ++j, ++index) {
	xyz[index].coord[0] = lx + dx;
	xyz[index].coord[1] = j * dy;
    }

// Edges 3: y = 0 - dy
    for (int i = 1; i < Nx-1; ++i, ++index) {
	xyz[index].coord[0] = i * dx;
	xyz[index].coord[1] = -dy;
    }

// Edges 4: y = ly + dy
    for (int i = 1; i < Nx-1; ++i, ++index) {
	xyz[index].coord[0] = i * dx;
	xyz[index].coord[1] = ly + dy;
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
    
    return 0;
}
/* */








#endif // _PARECTANGLEKNOTS_HPP_

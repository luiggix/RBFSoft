
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

#ifndef _RECTBLOCKSKNOTS_HPP_
#define _RECTBLOCKSKNOTS_HPP_

#include <fstream>
#include "Traits.hpp"
#include "Knots.hpp"
#include "knot2D.hpp"

/*!
 ***************************************************************************
 *  Rectangular domain with blocks near boundaries.
 *  This class is aimed to generate knots inside of rectangular domains with
 *  regular distribution near the boundaries. 
 *  The function constructKnots() must be reimplemented here to generate the 
 *  knots on the rectangular domain.
 *  \author  Luis M. de la Cruz [ Tue Jan 22 15:24:58 GMT 2008 ]
 ***************************************************************************
 */
template < typename Tprec>
class RectBlocksKnots : public Knots< RectBlocksKnots<Tprec>, knot2D > 
{ 
// Precision type (Tprec is for template, prec_t is cool in emacs ;) )
    typedef Tprec prec_t;

// These variables are from the Knots base class and must be defined
// in each derived class.
    using Knots<RectBlocksKnots<prec_t>, knot2D >::N;   // Total points
    using Knots<RectBlocksKnots<prec_t>, knot2D >::NI;  // Interior points
    using Knots<RectBlocksKnots<prec_t>, knot2D >::NB;  // Boundary points
    using Knots<RectBlocksKnots<prec_t>, knot2D >::Dim; // Dimension
    using Knots<RectBlocksKnots<prec_t>, knot2D >::xyz; // Coordinates
    using Knots<RectBlocksKnots<prec_t>, knot2D >::max_upper_range; // used for
                                                                   // KDTREE

private:
    prec_t lx;  ///< Horizontal length of the rectangle
    prec_t ly;  ///< Vertical length of the rectangle
    int Nx;     ///< Number of points on x-axis
    int Ny;     ///< Number of points on y-axis
// R = right, L = left, N = north, S = south
    int NBR, NBL, NBN, NBS; //< Number of nodes on the "boundary layer"
    int NIx, NIy;           //< Nodes for the inner block
    random_t knots_dist;
    prec_t epru;  ///< Randomness of the points distribution (default ep = 0)
    int NG;       ///< Number of ghost points around de domain.

public:

/*!
 *  Default constructor.
 */
    RectBlocksKnots() : Knots< RectBlocksKnots<Tprec>, knot2D > () { }
/*!
 *  Constructor to generate knots regularly spaced near the boundaries
 *  and random our uniform knots in the internal block.
 *  \param l_x Length of the rectangle in x-axis direction.
 *  \param n_x Number of points in x-axis direction.
 *  \param l_y Length of the rectangle in y-axis direction.
 *  \param n_y Number of points in y-axis direction.
 *  \param knd distribution of points (0 unif, 1 rand, 2 rand unif).
 *  \param layer 0 no layers around the domain, 1 just one layer.
 *  \param nbl x-length of the block 1
 *  \param nbr x-length of the block 2
 *  \param nbs y-length of the block 3
 *  \param nbn y-length of the block 4
 *  \param nix Number of points in x-axis for the internal block
 *  \param niy Number of points in y-axis for the internal block 
 *  
 */
    RectBlocksKnots(prec_t l_x, int n_x, prec_t l_y, int n_y, 
		    random_t knd = U, int layer = 0,
		    int nbl = 0, int nbr = 0, int nbs = 0, int nbn = 0, 
		    int nix = -1, int niy = -1)
	: Knots< RectBlocksKnots<Tprec>, knot2D > (),  
	  lx(l_x), Nx (n_x), ly(l_y), Ny(n_y), knots_dist(knd),
	  NBL(nbl), NBR(nbr), NBS(nbs), NBN(nbn), 
	  NIx(nix), NIy(niy) {
	
	Dim = 2;            // A rectangle is defined in 2D.

	NB = 2 * (Nx + Ny) - 4;

	if (nix == -1 && niy == -1) { NIx = Nx - 4; NIy = Ny - 4; }

// Vertical
	NI = NIx * NIy + // Inner block
	    (Ny - 2) * (NBL + NBR) +  // 1 and 2 blocks
	    (Nx - 2 - NBL - NBR) * (NBS + NBN); // 3 and 4 blocks

/* 
// Horizontal
	NI = NIx * NIy + // Inner block
	    (Nx - 2) * (NBS + NBN) +  // 1 and 2 blocks
	    (Ny - 2 - NBS - NBN) * (NBL + NBR); // 3 and 4 blocks
*/

/*
	if ( !NBL ) NB += ( -Ny + NBS + NBN + NIy );
	if ( !NBR ) NB += ( -Ny + NBS + NBN + NIy );
	if ( !NBS ) NB += ( -Nx + NBL + NBR + NIx );
	if ( !NBN ) NB += ( -Nx + NBL + NBR + NIx );
*/

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
	epru = 0.1;
    } 


    ~RectBlocksKnots() { }
/*!
 *  Construct all the knots of the domain.
 */
    bool constructKnots();
/*!
 *  Add special points for the ACBF aproximation.
 */
    inline bool addSpecialPoints(iMat &, int, int);
/*!
 *  Random parameter
 */ 
    void setRandomness(prec_t e) { 	
	epru = (1 - e) * 0.5 ;
    }
/*!
 *  Print info to the standard output (screen).
 */ 
    bool print();
};

template<typename prec_t> 
bool RectBlocksKnots<prec_t>::print() 
{
    std::string dist_type = "UNIFORM";
    if (knots_dist == R) dist_type = "RANDOM";
    if (knots_dist == RU) dist_type = "RANDOM-CELLS"; 

    Knots<RectBlocksKnots, knot2D >::print();
    std::cout << "\n +-----[ RectBlocks information ]-----+ " 
	      << "\n | RBF : Ghost points = " << NG
	      << "\n | RBF : Length in x axis = " << lx
	      << "\n | RBF : Length in y axis = " << ly
	      << "\n | RBF : Boundary points in x-axis = " << Nx
	      << "\n | RBF : Boundary points in y-axis = " << Ny
	      << "\n | RBF : Knots distribution = " << dist_type
	      << "\n | RBF : Randomness (epsilon [0 to 1]) = " << 1 - 2 * epru
	      << "\n | RBF : NBL = " << NBL << "\t NBR = " << NBR
	      << "\n | RBF : NBS = " << NBS << "\t NBN = " << NBN
	      << "\n | RBF : NIx = " << NIx << "\t NIy = " << NIy
	      << "\n | RBF : Boundary block 1 = " << NBS * (Nx-2)
	      << "\n | RBF : Boundary block 2 = " << NBN * (Nx-2)
	      << "\n | RBF : Boundary block 3 = " << NBL * (Ny-2-NBN-NBS)
	      << "\n | RBF : Boundary block 4 = " << NBR * (Ny-2-NBN-NBS)
	      << "\n | RBF : Total block knots= " << 
NBS * (Nx-2) + NBN * (Nx-2) + NBL * (Ny-2-NBN-NBS) + NBR * (Ny-2-NBN-NBS)
	      << "\n | RBF : Inner block      = " << NIx * NIy
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
 *  its cell, the randomness depends on the \c ep parameter. Also, it is 
 *  possible to generate blocks of uniform-gridded points near the boundaries
 *  setting the parameters \c nbl , \c nbr , \c nbs and \c nbn . Look at the
 *  next figure: 
 * \image html  RectBlocksKnots.png "Knots order" width=5cm 
 * \image latex RectBlocksKnots.eps "Knots order" width=5cm  
 */
template<typename prec_t>
bool RectBlocksKnots<prec_t>::constructKnots() 
{
    prec_t dx = lx / (Nx - 1);
    prec_t dy = ly / (Ny - 1);
    prec_t shift_l = dx * NBL, shift_r = dx * NBR;
    prec_t shift_s = dy * NBS, shift_n = dy * NBN;

    prec_t lx_I = lx - shift_l - shift_r;
    prec_t ly_I = ly - shift_s - shift_n;
    prec_t dx_I = lx_I / NIx;
    prec_t dy_I = ly_I / NIy;
    prec_t x0, x1, y0, y1;    
    prec_t ep = dx > dy ? dy * 0.25 : dx * 0.25;

// ----- Interior points.
    int index = 0; // The xyz is a C-array, so begins in 0.

    switch (knots_dist) 
    {
	case 0: // Uniform
	    for (int j = 1; j <= NIy; ++j)
		for (int i = 1; i <= NIx; ++i, ++index) {
		    xyz[index].coord[0] = shift_l + dx_I * (i - 0.5);
		    xyz[index].coord[1] = shift_s + dy_I * (j - 0.5);
		}
	    break;
	case 1: // Random
	    srand( time(0) );	    
	    
	    for (int j = 1; j <= NIy; ++j)
		for (int i = 1; i <= NIx; ++i, ++index) {
		    xyz[index].coord[0] = shift_l + ep +
			(1-2*ep) * lx_I * static_cast<prec_t>(rand())/RAND_MAX;
		    xyz[index].coord[1] = shift_s + ep +
			(1-2*ep) * ly_I * static_cast<prec_t>(rand())/RAND_MAX;
		}	    
	    break;
	case 2: // Random controlled
	    srand( time(0) );
	    for (int j = 1; j <= NIy; ++j)
		for (int i = 1; i <= NIx; ++i, ++index) {
		    x0 = (i - 1 + epru) * dx_I; x1 = (i - epru) * dx_I;
		    y0 = (j - 1 + epru) * dy_I; y1 = (j - epru) * dy_I;
		    xyz[index].coord[0] = shift_l + x0 + 
			( static_cast<prec_t>(rand())/RAND_MAX ) * (x1 - x0);
		    xyz[index].coord[1] = shift_s + y0 +
			( static_cast<prec_t>(rand())/RAND_MAX ) * (y1 - y0);
		}	    
	    break;
	default:
	    break;
    }

// ----- Boundary blocks

/* */
// Block 1: vertical
    for (int i = 1; i <= NBL; ++i)
	for (int j = 1; j <= Ny - 2 ; ++j, ++index) {
	    xyz[index].coord[0] = i * dx;
	    xyz[index].coord[1] = j * dy;
	}

// Block 2: vertical
    prec_t shift = lx - shift_r;
    for (int i = 0; i < NBR; ++i)
	for (int j = 1; j <= Ny - 2; ++j, ++index) {
	    xyz[index].coord[0] = shift + i * dx;
	    xyz[index].coord[1] = j * dy;
	}

// Block 3: horizontal
    int ei = Nx - 2 - NBL - NBR;
    for (int j = 1; j <= NBS; ++j)
	for (int i = 1; i <= ei; ++i, ++index) {
	    xyz[index].coord[0] = shift_l + i * dx;
	    xyz[index].coord[1] = j * dy;
	}

// Block 4: horizontal
    shift = ly - shift_n;
    for (int j = 0; j < NBN; ++j)
	for (int i = 1; i <= ei; ++i, ++index) {
	    xyz[index].coord[0] = shift_l + i * dx;
	    xyz[index].coord[1] = shift + j * dy;
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

// Edge 1 and 2: 
    for(int i = 0; i <= 1 ; ++i)
	for(int j = 1; j <= Ny - 2; ++j, ++index) {
	    xyz[index].coord[0] = lx * i;
	    xyz[index].coord[1] = j * dy;
	}

// Edge 3 and 4:
    for(int j = 0; j <= 1 ; ++j)
	for(int i = 1; i <= Nx - 2; ++i, ++index) {
	    xyz[index].coord[0] = i * dx;
	    xyz[index].coord[1] = ly * j ;
	} 

// Corners:
// 1
    xyz[index].coord[0] = 0;
    xyz[index].coord[1] = 0;

// 2
    xyz[++index].coord[0] = 0;
    xyz[index].coord[1] = ly;
// 3
    xyz[++index].coord[0] = lx;
    xyz[index].coord[1] = 0;   
// 4
    xyz[++index].coord[0] = lx;
    xyz[index].coord[1] = ly;

/* * //---------------------------------------------------------------
// Block 1: horizontal
    for (int j = 1; j <= NBS; ++j)
	for (int i = 1; i <= Nx - 2; ++i, ++index) {
	    xyz[index].coord[0] = i * dx;
	    xyz[index].coord[1] = j * dy;
	}

// Block 2: horizontal
    prec_t shift = ly - shift_n; 
    for (int j = 0; j < NBN - 1; ++j)
	for (int i = 1; i <= Nx - 2; ++i, ++index) {
	    xyz[index].coord[0] = i * dx;
	    xyz[index].coord[1] = shift + j * dy;
	}

// Block 3: vertical
    int ej = Ny - 2 - NBN - NBS;
    for (int i = 1; i <= NBL; ++i)
	for (int j = 1; j <= ej; ++j, ++index) {
	    xyz[index].coord[0] = i * dx;
	    xyz[index].coord[1] = shift_s + j * dy;
	}

// Block 4: vertical
    shift = lx - shift_r;
    for (int i = 0; i < NBR; ++i)
	for (int j = 1; j <= ej; ++j, ++index) {
	    xyz[index].coord[0] = shift + i * dx;
	    xyz[index].coord[1] = shift_s + j * dy;
	}

// Edge 1 and 2:
    for(int j = 0; j <= 1 ; ++j)
	for(int i = 1; i <= Nx - 2; ++i, ++index) {
	    xyz[index].coord[0] = i * dx;
	    xyz[index].coord[1] = ly * j ;
	} 

// Edge 3 and 4: 
    for(int i = 0; i <= 1 ; ++i)
	for(int j = 1; j <= Ny - 2; ++j, ++index) {
	    xyz[index].coord[0] = lx * i;
	    xyz[index].coord[1] = j * dy;
	}


// Corners:
// 1
    xyz[index].coord[0] = 0;
    xyz[index].coord[1] = 0;
// 2
    xyz[++index].coord[0] = lx;
    xyz[index].coord[1] = 0;   
// 3
    xyz[++index].coord[0] = 0;
    xyz[index].coord[1] = ly;
// 4
    xyz[++index].coord[0] = lx;
    xyz[index].coord[1] = ly;
/* */
    return 0;
}

/*!
 *  Add special points located in the boundary of the rectangle and 
 *  in the geometric centroid.
 *  \param P  indices of neighbors knots.
 *  \param j  actual knot.
 *  \param nsp number of special points.
 */
//NOTE: Remember that the indices in the matrix P begin in 1, and the
//      indices in xyz[] begin in 0. 
template<typename prec_t>
inline bool RectBlocksKnots<prec_t>::addSpecialPoints(iMat &P, int j, int nsp)
{
    int sp[9], side_x = Nx - 2, side_y = Ny - 2;
// Corners
    sp[0] = N - 4; 
    sp[1] = sp[0] + 1;
    sp[2] = sp[1] + 1;
    sp[3] = sp[2] + 1;

// Middle Edges
    sp[4] = NI + side_y / 2;
    sp[5] = sp[4] + side_y; 
    sp[6] = NI + 2 * side_y + side_x / 2;
    sp[7] = sp[6] + side_x; 

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



#endif // _RECTBLOCKSKNOTS_HPP_

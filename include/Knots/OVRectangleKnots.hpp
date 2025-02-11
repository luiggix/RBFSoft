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


#ifndef _OVRECTANGLEKNOTS_HPP_
#define _OVRECTANGLEKNOTS_HPP_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "Traits.hpp"
#include "DD/CartComm.hpp"
#include "Knots.hpp"
#include "knot2D.hpp"

/*!
 ***************************************************************************
 *  OVRectangular domain.
 *  This class is aimed to generate knots inside overlapping rectangular 
 *  domains used by Domain Decomposition algorithms. The function 
 *  constructKnots() must be reimplemented here to generate the knots on the 
 *  rectangular domain.
 *  \author  Luis M. de la Cruz [ Mon Nov 12 14:53:34 GMT 2007 ]
 ***************************************************************************
 */
template < typename Tprec>
class OVRectangleKnots : public Knots< OVRectangleKnots<Tprec>, knot2D > 
{ 
// Precision type (Tprec is for template, prec_t is cool in emacs ;) )
    typedef Tprec prec_t;

// These variables are from the Knots base class and must be defined
// in each derived class.
    using Knots<OVRectangleKnots<prec_t>, knot2D >::N;   // Total points
    using Knots<OVRectangleKnots<prec_t>, knot2D >::NI;  // Interior points
    using Knots<OVRectangleKnots<prec_t>, knot2D >::NB;  // Boundary points
    using Knots<OVRectangleKnots<prec_t>, knot2D >::Dim; // Dimension
    using Knots<OVRectangleKnots<prec_t>, knot2D >::xyz; // Coordinates
    using Knots<OVRectangleKnots<prec_t>, knot2D >::max_upper_range; // used for
                                                                   // KDTREE

private:
    prec_t lx;  ///< Horizontal length of the rectangle
    prec_t ly;  ///< Vertical length of the rectangle
    prec_t dx;
    prec_t dy;
    int Nx;     ///< Number of points on x-axis
    int Ny;     ///< Number of points on y-axis    
    random_t knots_dist; ///< kind of point distribution
    prec_t epru;  ///< Randomness of the points distribution (default epru = 0)
    int rank; // For parallel implementation, id of the procesor
    int I, J, NP_I, NP_J;
    prec_t shift[2], overlap[2];
    Vec array_L, array_R, array_U, array_D;

public:

/*!
 *  Default constructor.
 */
    OVRectangleKnots() : Knots< OVRectangleKnots<Tprec>, knot2D > () { }
/*!
 *  Constructor to generate knots in the centers of the "cells" of
 *  a gridded mesh. It allows some randomness of the points.
 *  \param l_x Global length of the rectangle in x-axis direction.
 *  \param l_y Global length of the rectangle in y-axis direction.
 *  \param n_x Local number of points in x-axis 
 *  \param n_y Local number of points in y-axis direction.
 *  \param g_x overlapping in x-axis.
 *  \param g_y overlapping in y-axis.
 *  \param cart Communicator between processors (CartComm).
 *  \param knd kind of distribution for the points.
 */
    OVRectangleKnots(prec_t l_x, int n_x, prec_t l_y, int n_y, 
		     prec_t g_x, prec_t g_y,
		     CartComm& cart, random_t knd = U) 
	: Knots< OVRectangleKnots<Tprec>, knot2D > (), knots_dist(knd) {

// ----- Get the correct lenght, num. of points, shifts and overlapping of 
// ----- every subdomain, based on its position in the CartComm
// The numeration of CartComm is as in matrices of Linear Algebra
// therefore the I goes in -y direction, while J goes in x direction
// (I,J) \in [0, NP_I - 1] X [0, NP_J - 1]
//
	rank = cart.comm.Get_rank();
	NP_I = cart.getNumProc_I();
	NP_J = cart.getNumProc_J();
	I = cart.get_I();
	J = cart.get_J();

	prec_t llx = l_x / NP_J; // Size of the subdomain not
	prec_t lly = l_y / NP_I; // including the overlapping

	lx = l_x / NP_J + 2 * g_x; // Same size for all subdom, except for
	ly = l_y / NP_I + 2 * g_y; // the overlapping.
	overlap[0] = g_x;  
	overlap[1] = g_y;  
	Nx = n_x;        
	Ny = n_y;        

// Redimensioning the lenght depending on the location of the subdomain	
	if ( I == 0)         ly -= g_y;
	if ( I == NP_I - 1 ) ly -= g_y;
	if ( J == 0)         lx -= g_x;
	if ( J == NP_J - 1 ) lx -= g_x;

	dx = lx / (Nx - 1); 
	dy = ly / (Ny - 1); 
	shift[0] = llx * J; 
	shift[1] = lly * (NP_I - I - 1);

// Correction of the shift by the overlapping
	if ( I != NP_I - 1 ) { shift[1] -= g_y; }
	if ( J != 0 ) { shift[0] -= g_x; }

// ----- Now define some parameters to construct the rectangle

        Dim = 2;            // A rectangle is defined in 2D. 
	NB = 2 * (Nx + Ny) - 4;

	if (knots_dist == RU)
	    NI = (Nx - 1) * (Ny - 1); // Random wrapped in cells
	else
	    NI = (Nx - 2) * (Ny - 2); // Uniform or Random

	N = NI + NB; 
	xyz = new knot2D[N];
	
// Set the max uppper range for the binary search of the RANGE in kdtree
	max_upper_range = (ly > lx) ? lx : ly;

// Default randomnesss value:
	epru = 0.1;
    } 

    ~OVRectangleKnots() { }
/*!
 *  Construct all the knots of the domain.
 */
    bool constructKnots();
    bool constructKnotsDriven();

    inline bool constructOverlapping(CartComm&);

/*!
 *  Add special points for the ACBF aproximation.
 */
    inline bool addSpecialPoints(iMat &, int, int);
/*!
 *  Random parameter
 */ 
    inline void setRandomness(prec_t e) { epru = (1 - e) * 0.5 ; }
    inline prec_t getShapeParameter() { return 1.0 / sqrt(N); }
    inline int getNx() { return Nx; }
    inline int getNy() { return Ny; }
    inline prec_t getLx() { return lx; }
    inline prec_t getLy() { return ly; }
    inline prec_t getDx() { return dx; }
    inline prec_t getDy() { return dy; }    
    inline prec_t getShift(axis_t a) { return shift[a]; }
    inline prec_t getOverlap(axis_t a) {
	if ( I == NP_I - 1 && a == Y ) return 0;       
	if ( J == 0        && a == X ) return 0;       	
	return overlap[a];
    }

    inline Vec getArrayL() { return array_L; }
    inline Vec getArrayR() { return array_R; }
    inline Vec getArrayU() { return array_U; }
    inline Vec getArrayD() { return array_D; }

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
bool OVRectangleKnots<prec_t>::print() 
{
    std::string dist_type = "UNIFORM";
    if (knots_dist == R) dist_type = "RANDOM";
    if (knots_dist == RU) dist_type = "RANDOM-CELLS";

    Knots<OVRectangleKnots, knot2D >::print();
    std::cout << "\n +-----[ OVRectangle information ]-----+ " 
	      << "\n | RBF : Array of subdoms, rows x cols = " 
	      << NP_I << " x " << NP_J
	      << "\n | RBF : (I,J) = (" << I << "," << J << ")"
	      << "\t Rank = " << rank 
	      << "\n | RBF : Length in x axis = " << lx
	      << "\n | RBF : Length in y axis = " << ly
	      << "\n | RBF : Overlap in x-axis = " << overlap[0]
	      << "\n | RBF : Overlap in y-axis = " << overlap[1]
	      << "\n | RBF : Shift in x-axis = " << shift[0]
	      << "\n | RBF : Shift in y-axis = " << shift[1]
	      << "\n | RBF : Boundary points in x-axis = " << Nx
	      << "\n | RBF : Boundary points in y-axis = " << Ny
	      << "\n | RBF : Knots distribution = " << dist_type
	      << "\n | RBF : Randomness (epsilon [0 to 1]) = " << 1 - 2 * epru
	      << "\n +-----+";
    return 0;
}

/*!
 * Construction of the knots.
 */
template<typename prec_t>
bool OVRectangleKnots<prec_t>::constructKnots() 
{
    prec_t x0, x1, y0, y1;
    prec_t ep = dx > dy ? dy * 0.1 : dx * 0.1;

// ----- Interior points.
    int ei = Nx - 2, ej = Ny - 2;

    int index = 0; // The xyz is a C-array, so it begins in 0.

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

// ----- Ghost Boundaries or real boundaries, depends on the subdomain

// Edges 1: x = shift[0]
    for (int j = 1; j < Ny - 1; ++j, ++index) {
	xyz[index].coord[0] = shift[0];
	xyz[index].coord[1] = j * dy + shift[1];
    }

// Edges 2: x = lx + shift[0]
    for (int j = 1; j < Ny - 1; ++j, ++index) {
	xyz[index].coord[0] = lx + shift[0];
	xyz[index].coord[1] = j * dy + shift[1];
    }

// Edges 3: y = shift[1]    
    for (int i = 1; i < Nx - 1; ++i, ++index) {
	xyz[index].coord[0] = i * dx + shift[0];
	xyz[index].coord[1] = shift[1];
    }

// Edges 4: y = ly + shift[1]    
    for (int i = 1; i < Nx - 1; ++i, ++index) {
	xyz[index].coord[0] = i * dx + shift[0];
	xyz[index].coord[1] = ly + shift[1];
    }

// Corners:
// 1
    xyz[N-4].coord[0] = shift[0];
    xyz[N-4].coord[1] = shift[1];

// 2
    xyz[N-3].coord[0] = shift[0];
    xyz[N-3].coord[1] = shift[1] + ly;

// 3
    xyz[N-2].coord[0] = shift[0] + lx;
    xyz[N-2].coord[1] = shift[1];

// 4
    xyz[N-1].coord[0] = shift[0] + lx;
    xyz[N-1].coord[1] = shift[1] + ly;

    return 0;
}

template<typename prec_t>
inline bool OVRectangleKnots<prec_t>::constructOverlapping(CartComm& cart)
{
    MPI::Request req[8];
    MPI::Status  sta[8];

    const int *neighbor = cart.getNeighbors();

// ----- Exchange of the size of overlapping.
    prec_t gx_L, gx_R, gy_U, gy_D;
    req[0] = cart.comm.Isend(&overlap[0], 1, MPI::DOUBLE, neighbor[LEFT],  1);
    req[1] = cart.comm.Isend(&overlap[0], 1, MPI::DOUBLE, neighbor[RIGHT], 1);
    req[2] = cart.comm.Isend(&overlap[1], 1, MPI::DOUBLE, neighbor[UP],    1);
    req[3] = cart.comm.Isend(&overlap[1], 1, MPI::DOUBLE, neighbor[DOWN],  1);

    req[4] = cart.comm.Irecv(&gx_L, 1, MPI::DOUBLE, neighbor[LEFT],  1);
    req[5] = cart.comm.Irecv(&gx_R, 1, MPI::DOUBLE, neighbor[RIGHT], 1);
    req[6] = cart.comm.Irecv(&gy_U, 1, MPI::DOUBLE, neighbor[UP],    1);
    req[7] = cart.comm.Irecv(&gy_D, 1, MPI::DOUBLE, neighbor[DOWN],  1);

    MPI::Request::Waitall(8, req, sta);    

// ----- Exchange of the number of points for the overlapping.
    int Ny_L, Ny_R, Nx_U, Nx_D;
    req[0] = cart.comm.Isend(&Ny, 1, MPI::INT, neighbor[LEFT],  2);
    req[1] = cart.comm.Isend(&Ny, 1, MPI::INT, neighbor[RIGHT], 2);
    req[2] = cart.comm.Isend(&Nx, 1, MPI::INT, neighbor[UP],    2);
    req[3] = cart.comm.Isend(&Nx, 1, MPI::INT, neighbor[DOWN],  2);

    req[4] = cart.comm.Irecv(&Ny_L, 1, MPI::INT, neighbor[LEFT],  2);
    req[5] = cart.comm.Irecv(&Ny_R, 1, MPI::INT, neighbor[RIGHT], 2);
    req[6] = cart.comm.Irecv(&Nx_U, 1, MPI::INT, neighbor[UP],    2);
    req[7] = cart.comm.Irecv(&Nx_D, 1, MPI::INT, neighbor[DOWN],  2);

    MPI::Request::Waitall(8, req, sta);

// ----- Construction of overlapping zones for the neighbor subdomains.
//       array_L, array_R, array_U and array_D will contain the coordinates
//       of the overlapping lines. The last element of these arrays contains
//       the x- or y-coordinate, depending on the orientation of the line.
    prec_t dh, xya;
    if ( neighbor[LEFT] != -1)  {  // Edge 1
	array_L.resize(Ny_L + 1);
	xya = shift[0] + overlap[0] + gx_L;
	dh = ly / (Ny_L - 1);
	for(int i = 1; i <= Ny_L; ++i) 
	    array_L(i) = (i - 1) * dh + shift[1];
	array_L(Ny_L + 1) = xya; // the x-coordinate
    }

    if ( neighbor[RIGHT] != -1) {  // Edge 2
	array_R.resize(Ny_R + 1);
	xya = shift[0] + lx - overlap[0] - gx_R;
	dh = ly / (Ny_R - 1);
	for(int i = 1; i <= Ny_R; ++i)
	    array_R(i) = (i - 1) * dh + shift[1];
	array_R(Ny_R + 1) = xya; // the x-coordinate
    }

    if ( neighbor[UP] != -1)    {  // Edge 3
	array_U.resize(Nx_U + 1);
	xya = shift[1] + ly - overlap[1] - gy_U;
	dh = lx / (Nx_U - 1);
	for(int i = 1; i <= Nx_U; ++i)
	    array_U(i) = (i - 1) * dh + shift[0];
	array_U(Nx_U + 1) = xya; // the y-coordinate
    }

    if ( neighbor[DOWN] != -1)  {  // Edge 4
	array_D.resize(Nx_D + 1);
	xya = shift[1] + overlap[1] + gy_D;
	dh = lx / (Nx_D - 1);
	for(int i = 1; i <= Nx_D; ++i)
	    array_D(i) = (i - 1) * dh + shift[0];
	array_D(Nx_D + 1) = xya; // the y-coordinate
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
inline bool OVRectangleKnots<prec_t>::addSpecialPoints(iMat &P, int j, int nsp)
{
    int sp[9], side_x = Nx / 2, side_y = Ny / 2;
// Corners
    sp[0] = N - 4;      // 1
    sp[1] = sp[0] + 1;  // 2
    sp[2] = sp[1] + 1;  // 3
    sp[3] = sp[2] + 1;  // 4

    int sx_sy = side_x + side_y;
// Middle Edges
    sp[4] = NI + side_y;  // Edge 1 
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
bool OVRectangleKnots<prec_t>::gridKnots() 
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


#endif // _OVRECTANGLEKNOTS_HPP_

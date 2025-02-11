
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

#ifndef _SEMICIRCLEKNOTS_HPP_
#define _SEMICIRCLEKNOTS_HPP_

#include <fstream>
#include "Traits.hpp"
#include "Knots.hpp"
#include "knot2D.hpp"

#include <cstdio>
#include <cstdlib>


/*!
 ***************************************************************************
 *  SemiCircle domain.
 *  This class is aimed to generate knots inside of semicircle domains. 
 *  The function constructKnots() must be reimplemented here to generate the 
 *  knots in the semicircle domain.
 *  \author  Luis M. de la Cruz [ Fri Jan 25 13:48:24 GMT 2008 ]
 ***************************************************************************
 */
template < typename Tprec>
class SemiCircleKnots : public Knots< SemiCircleKnots<Tprec>, knot2D > 
{ 
// Precision type (Tprec is for template, prec_t is cool in emacs ;) )
    typedef Tprec prec_t;

// These variables are from the Knots base class and must be defined
// in each derived class.
    using Knots<SemiCircleKnots<prec_t>, knot2D >::N;   // Total points
    using Knots<SemiCircleKnots<prec_t>, knot2D >::NI;  // Interior points
    using Knots<SemiCircleKnots<prec_t>, knot2D >::NB;  // Boundary points
    using Knots<SemiCircleKnots<prec_t>, knot2D >::Dim; // Dimension
    using Knots<SemiCircleKnots<prec_t>, knot2D >::xyz; // Coordinates
    using Knots<SemiCircleKnots<prec_t>, knot2D >::max_upper_range; // used for
                                                                   // KDTREE

private:
    prec_t theta_1;  ///< Initial angle
    prec_t theta_2;  ///< Final angle
    prec_t radius_1; ///< Initial radius
    prec_t radius_2; ///< Final radius
    int Ntheta; ///< Number of nodes in theta
    int Nradius; ///< Number of nodes in radius
    random_t knots_dist;
    prec_t epru; ///< Used to control randomness in regular distr.

public:

/*!
 *  Default constructor.
 */
    SemiCircleKnots() : Knots< SemiCircleKnots<Tprec>, knot2D > () { }
/*!
 *  Constructor to generate knots in the centers of the "cells" of
 *  a gridded mesh, totally random knots or controlled random knots. 
 *  \param r1 radius 1 as shown in the next figure.
 *  \param r2 radius 2 as shown in the next figure.
 *  \param t1 angle in degrees 1 as shown in the next figure.
 *  \param t2 angle in degrees 2 as shown in the next figure.
 *  \param Nr Number of points in radius.
 *  \param Nt Number of points in theta.
 *  \param knd distribution of points (0 unif, 1 rand, 2 rand unif).
 * \image html  SemiCircle.png "Semicircle definition." width=5cm 
 * \image latex SemiCircle.eps "Semicircle definition." width=5cm
 */
    SemiCircleKnots(prec_t r1, prec_t r2, int Nr, 
		    prec_t t1, prec_t t2, int Nt, 
		    random_t knd = U)
      : Knots< SemiCircleKnots<Tprec>, knot2D > (),  
	theta_1(t1), theta_2(t2), radius_1(r1), radius_2(r2),
	Ntheta(Nt), Nradius(Nr), knots_dist(knd) {

        Dim = 2;            // A semicircle is defined in 2D.
 
	NB = 2 * (Ntheta + Nradius) - 4;
// Uniform or total random knots distribution
	if (knots_dist == U || knots_dist == R) 
	    NI = (Ntheta - 2) * (Nradius - 2);
// Random knots inside cells
	else if (knots_dist = RU)
	    NI = (Ntheta - 1) * (Nradius - 1);

	N = NI + NB;
	xyz = new knot2D[N];

// Transform degrees to radians
	prec_t dg_to_ra = PI / 180.0;	
	theta_1 = dg_to_ra * theta_1;
	theta_2 = dg_to_ra * theta_2;

	prec_t lt = 2 * PI * radius_1 * (theta_2 - theta_1);
	prec_t rd = radius_2 - radius_1;
// Set the max uppper range for the binary search of the RANGE in kdtree
	max_upper_range = (rd > lt) ? lt : rd;	

// Default randomnesss value:
	epru = 0.1;
    } 

    ~SemiCircleKnots() { }
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
	epru = (1 - e) * 0.5 ;}
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
bool SemiCircleKnots<prec_t>::print() 
{
    Knots<SemiCircleKnots, knot2D >::print();
    std::cout << "\n +-----[ SemiCircle information ]-----+ " 
	      << "\n | RBF : Radius 1 = " << radius_1
	      << "\n | RBF : Radius 2 = " << radius_2
	      << "\n | RBF : Angle 1 = " << theta_1
	      << "\n | RBF : Angle 2 = " << theta_2
	      << "\n | RBF : Boundary points in theta = " << Ntheta
	      << "\n | RBF : Boundary points in radius = " << Nradius
	      << "\n | RBF : Randomness (epsilon [0 to 1]) = " << 1 - 2 * epru
	      << "\n +-----+";
    return 0;
}

/*!
 *  This function generate nodes in a uniform grid (knd = 0), randomly
 *  distributed (knd = 1) and random controlled ( knd = 2, 0 < epru < 1).
 *  The random controlled is generated as follows: Taking the values of \c Nr 
 *  and \c Nt an "invisible grid" is defined inside the semicircle, this grid 
 * will contains (\c Nr - 1) \c X (\c Nt - 1) "cells", for each cell a knot is 
 *  generated and this knot is allowed to be in a random position inside its 
 *  cell, the randomness depends on the \c ep parameter. The order of the
 *  knots in the \c xyz array is as shown in the next figure:
 * \image html  SemiCircleKnots.png "Knots order" width=5cm
 * \image latex SemiCircleKnots.eps "Knots order" width=5cm
 */
template<typename prec_t>
bool SemiCircleKnots<prec_t>::constructKnots() 
{
    prec_t dif_ra = radius_2 - radius_1; 
    prec_t dif_th = theta_2 - theta_1;
    prec_t dr = dif_ra / (Nradius - 1 );
    prec_t dt = dif_th / (Ntheta - 1);
    prec_t R, T, x0, x1, y0, y1;
    prec_t ep = dr > dt ? dt * 0.1 : dr * 0.1;

// ----- Interior points.
    int ei = Nradius - 2, ej = Ntheta - 2;

    int index = 0; // The xyz is a C-array, so begins in 0.

    switch (knots_dist) 
    {
	case 0: // Uniform 
	    for (int j = 0; j < ej ; ++j)
		for (int i = 0; i < ei; ++i, ++index) {   
		    R = radius_1 + (i+1) * dr; 
		    T = theta_1  + (j+1) * dt;
		    xyz[index].coord[0] = R * cos(T) ;
		    xyz[index].coord[1] = R * sin(T) ;
		}
	    break;
	case 1: // Random	
	    srand( time(0) );
	    for (int j = 0; j < ej ; ++j)
		for (int i = 0; i < ei; ++i, ++index) { 
		    R = radius_1 + 
			(1-ep) * dif_ra * static_cast<prec_t>(rand())/RAND_MAX;
		    T = theta_1 + 
			(1-ep) * dif_th * static_cast<prec_t>(rand())/RAND_MAX;
		    xyz[index].coord[0] = R * cos(T) ;
		    xyz[index].coord[1] = R * sin(T) ; 
		}
	    break;
	case 2: // Random inside "cells"
	    srand( time(0) );
	    ++ei; ++ej;
	    for (int j = 0; j < ej ; ++j)
		for (int i = 0; i < ei; ++i, ++index) { 
		    x0 = (i + epru) * dr; x1 = (i + 1 - epru) * dr;// Mapping
		    y0 = (j + epru) * dt; y1 = (j + 1 - epru) * dt;// to a cell
		    R = radius_1 + x0 +
			( static_cast<prec_t>(rand())/RAND_MAX ) * (x1 - x0);
		    T = theta_1  + y0 +
			( static_cast<prec_t>(rand())/RAND_MAX ) * (y1 - y0);
		    xyz[index].coord[0] = R * cos(T) ;
		    xyz[index].coord[1] = R * sin(T) ;  
		}
	    break;
	default:
	    break;

    }

/* */
// ----- Boundary points

// Edge 1: Radius = constant
    for (int j = 1; j < Ntheta - 1; ++j, ++index) {
	T = theta_1 + j * dt;	    
	xyz[index].coord[0] = radius_1 * cos(T) ;
	xyz[index].coord[1] = radius_1 * sin(T) ;
    }
// Edge 2: Radius = constant
    for (int j = 1; j < Ntheta - 1; ++j, ++index) {
	T = theta_1 + j * dt;	    
	xyz[index].coord[0] = radius_2 * cos(T) ;
	xyz[index].coord[1] = radius_2 * sin(T) ;
    }

// Edge 3: Theta = constant
    for (int i = 1; i < Nradius - 1; ++i, ++index) {
	R = radius_1 + i * dr; 
	xyz[index].coord[0] = R * cos(theta_1) ;
	xyz[index].coord[1] = R * sin(theta_1) ;
    }

// Edge 3: Theta = constant
    for (int i = 1; i < Nradius - 1; ++i, ++index) {
	R = radius_1 + i * dr; 
	xyz[index].coord[0] = R * cos(theta_2) ;
	xyz[index].coord[1] = R * sin(theta_2) ;
    }

/* */
// Corners

// 1
    xyz[index].coord[0] = radius_1 * cos(theta_1) ;
    xyz[index].coord[1] = radius_1 * sin(theta_1) ;

    index++;

// 2
    xyz[index].coord[0] = radius_1 * cos(theta_2) ;
    xyz[index].coord[1] = radius_1 * sin(theta_2) ;

    index++;

// 3
    xyz[index].coord[0] = radius_2 * cos(theta_1) ;
    xyz[index].coord[1] = radius_2 * sin(theta_1) ;
 
    index++;

// 4
    xyz[index].coord[0] = radius_2 * cos(theta_2) ;
    xyz[index].coord[1] = radius_2 * sin(theta_2) ;

/* */
    return 0;
}

/*!
 *  Add special points located in the boundary of the semicircle and 
 *  in the geometric centroid.
 *  \param P  indices of neighbors knots.
 *  \param j  actual knot.
 *  \param nsp number of special points.
 */
//NOTE: Remember that the indices in the matrix P begin in 1, and the
//      indices in xyz[] begin in 0. 
template<typename prec_t>
inline bool SemiCircleKnots<prec_t>::addSpecialPoints(iMat &P, int j, int nsp)
{
    int sp[9], side_x = Nradius - 2, side_y = Ntheta - 2;
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


    int bi = 0; //P.firstRow();
    int ei = P.rows(); //P.lastRow();
    
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



template<typename prec_t>
bool SemiCircleKnots<prec_t>::gridKnots() 
{
    prec_t dr = (radius_2 - radius_1) / (Nradius - 1 );
    prec_t dt = (theta_2 - theta_1) / (Ntheta - 1);
    prec_t R, T;
    int index = 0; // The xyz is a C-array, so begins in 0.
    for (int j = 0; j < Ntheta ; ++j)
	for (int i = 0; i < Nradius; ++i, ++index) {   
	    R = radius_1 + i * dr; 
	    T = theta_1  + j * dt;
	    xyz[index].coord[0] = R * cos(T) ;
	    xyz[index].coord[1] = R * sin(T) ;
	}
    
    return 0;
}
/* */

#endif // _SEMICIRCLEKNOTS_HPP_

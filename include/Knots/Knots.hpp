
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

#ifndef _KNOTS_HPP_
#define _KNOTS_HPP_

#include <fstream>
#include <string>
#include <vector>
#include "Traits.hpp"
#include "kdtree++/kdtree.hpp"

//  This is a hack to deal with 1D, 2D and 3D, so for example, to use 
//  knot2D I need to add -D_2D_ to the compiler options (-D_3D_ for knot3D
//  and -D_1D_ for knot1D). Note that these options are mutually exclusive. 
//  I haven't found a fancy solution to this issue...
#ifdef _1D_
#include "Knots/knot1D.hpp"
#endif
#ifdef _2D_
#include "Knots/knot2D.hpp"
#endif
#ifdef _3D_
#include "Knots/knot3D.hpp"
#endif

/*!
 ***************************************************************************
 *  Base class for knots generation.
 *  This class contains the definition of general variables and methods to
 *  be used inside more particular classes. The template parameters are:
 *  \c Tdomain type of the domain of the derived class.
 *  \c Tknot type of the knot (knot2D or knot3D).
 *  \author  Luis M. de la Cruz [ Mon Dec 17 11:38:56 GMT 2007 ]
 *  \todo To find a way to generate knot for 1D, 2D or 3D without -D_1D_, 
 *  -D_2D_ and -D_3D_ compiler options.
 ***************************************************************************
 */
template<typename Tdomain, typename Tknot>
class Knots
{
protected:
    int N;      ///< Total number of knots.
    int NI;     ///< Iternal knots.
    int NB;     ///< knots on the boundary.
    int Dim;    ///< Dimension of the problem (1,2 or 3).
    Tknot *xyz; ///< Coordinates array.
    prec_t max_upper_range; ///< Every derived class must to define this bound.

public:
    Knots() : max_upper_range(1.0) { }
    ~Knots() {  }
/*! 
 *  Curiosly Recursive Template Pattern (CRTP).  
 *  Here we use the Curiosly Recursive Template Pattern (CRTP also known as
 *  the Barton-Nackman trick), to delegate responsabilities to derived 
 *  classes. The idea is to get the actual type of derived classes inside of
 *  the base class. With this trick, the compiler is able to define the correct
 *  function that must be called when an instance of the derived class calls
 *  to the general function, defined in the base class. This is also a way to
 *  construct static polymorphism, to avoid the dynamic one and get better 
 *  performance. See the C++ Templates - The Complete Guide by David 
 *  Vandevoorde and Nicolai M. Josuttis, Chapter 11, pag. 174. 
 */
    Tdomain& asDerived() { return static_cast<Tdomain&>(*this); }
/*!
 *  Delegate responsabilities to the derived classes.
 *  The responsability of knots construction is delegated to the subclasses.
 *  Every derived class must to define (reimplement) this function, and 
 *  every such class has the freedom of defining the way the knots are 
 *  constructed.
 */
    bool constructKnots() { return asDerived().constructKnots(); }    
/*!
 *  KDTree construction.
 *  The knots are stored in the \c tree variable
 *  defined in knot2D or knot3D, depending on the dimension of the problem.
 */
    bool calcKDTree() 
    {
// NOTE: The index range in the whole program is [1 to N], on the other
// hand, kdtree++ works with the range [0 to N-1], so I need to add 1.
	for (int i = 0; i < N; ++i) {
	    xyz[i].index = i + 1; 
	    tree.insert(xyz[i]); 
	}
    }
/*!
 *  Find the neighbors points of a defined knot inside a range (2D).
 *  \param x x-axis coordinate of the knot.
 *  \param y y-axis coordinate of the knot.
 *  \param range the size of the neighborhood where the neighbors live.
 *  \return \c neighbors_array will contain the list of neighbors 
 *          (it is a std::vector<Tknot>).
 */
// NOTE: The functions in kdtree++ don't look for the neighbors in spheres (3D)
// or circles (2D), but in squares(2D) or cubes(3D) of side defined by a range.
// The tree must be calculated before with the function calcKDTree().
    inline const std::vector<Tknot>& findNeighbors(const prec_t x, 
						   const prec_t y, 
						   const prec_t range) 
    {
	Tknot target;
	target.coord[0] = x;
	target.coord[1] = y;
	neighbors_array.clear();
	tree.find_within_range(target, range, 
			       std::back_inserter(neighbors_array));
	std::cout << "\n | " <<neighbors_array.size()<< " neighbors found ";

	return neighbors_array;
    }

/*!
 *  Find the neighbors points of a defined knot inside a range (1D).
 *  \param x x-axis coordinate of the knot.
 *  \param range the size of the neighborhood where the neighbors live.
 *  \return \c neighbors_array will contain the list of neighbors 
 *          (it is a std::vector<Tknot>).
 */
    inline const std::vector<Tknot>& findNeighbors(const prec_t x, 
						   const prec_t range) 
    {
	Tknot target;
	target.coord[0] = x;
	neighbors_array.clear();
	tree.find_within_range(target, range, 
			       std::back_inserter(neighbors_array));
	std::cout << "\n | " <<neighbors_array.size()<< " neighbors found ";
	return neighbors_array;
    }

/*!
 *  Find the neighbors points of a defined knot inside a range (3D).
 *  \param x x-axis coordinate of the knot.
 *  \param y y-axis coordinate of the knot.
 *  \param z z-axis coordinate of the knot.
 *  \param range the size of the neighborhood where the neighbors live.
 *  \return \c neighbors_array will contain the list of neighbors 
 *          (it is a std::vector<Tknot>).
 */
    inline const std::vector<Tknot>& findNeighbors(const prec_t x, 
						   const prec_t y,
						   const prec_t z,
						   const prec_t range) 
    {
	Tknot target;
	target.coord[0] = x;
	target.coord[1] = y;
	target.coord[2] = z;
	neighbors_array.clear();
	tree.find_within_range(target, range, 
			       std::back_inserter(neighbors_array));
	std::cout << "\n | " <<neighbors_array.size()<< " neighbors found ";
	return neighbors_array;
    }

/*!
 *  Find the neighbors to all knots in the domain, this information is used 
 *  by the ACBF preconditioner (acbfprecond). This works for 1D, 2D and 3D.
 *  \param j index of the working knot.
 *  \param P Integer Matrix containing the indices of neighbors points for
 *         every knot.
 *  \param max_knots Maximum number of knots allowed inside the range.
 *  \return 
 */
    inline bool findNeighbors(int j, iMat &P, int max_knots, int sp = 0) 
    {

// NOTE: The max_upper_range must be defined inside every derived class.
	prec_t upper_range = max_upper_range;
	prec_t lower_range = 0.0;
	prec_t range, rantol = 0.000001;
	int nneighbors; // number of neighbors.

// NOTE: the kdtree uses C-arrays style, so they begin in 0.
	int id = j - 1; 
	
// Binary search for the optimum range
	while ( true ) {
	    range = (upper_range + lower_range) * 0.5; 
	    nneighbors = tree.count_within_range(xyz[id], range);      
	    if ( nneighbors < max_knots )
		lower_range = range;
	    else if ( nneighbors > max_knots )  
		upper_range = range;
	    else 
		break;     
	    if((upper_range-lower_range < rantol) && (nneighbors > max_knots))
		break;
	}

// Find the neighbors and store them into an array	
	neighbors_array.clear();
	tree.find_within_range(xyz[id], range, 
			       std::back_inserter(neighbors_array));
	neighbors_i = neighbors_array.begin();
		
// NOTE: If the number of points is greater than the Max, the next code 
//       eliminate this exceedence.
	int dif = nneighbors - max_knots;
	int mid_max_knots = max_knots / 2;
	int position;
	if ( dif ) {
	    for (int n = 1; neighbors_i != neighbors_array.end(); 
		 ++neighbors_i, ++n)
		if (neighbors_i->index == j) { 
		    position = n; 
		    break; 
		}
	    
	    neighbors_i = neighbors_array.begin();		
	    
	    if ( position > mid_max_knots) {
		if ( nneighbors - position - mid_max_knots >= 0  ) 
		    neighbors_i += (position - mid_max_knots);	  
		else 
		    neighbors_i += (nneighbors - max_knots);	       
	    }
	}

// Store the global indices of max_knots neighbors to the i-esimo knot.
	for (int n = 1; n <= max_knots; ++neighbors_i, ++n) 
	    P(n, j) = neighbors_i->index;

// Next function must be defined in the derived classes.
	if ( sp )
	    asDerived().addSpecialPoints(P, j, sp);

/* *
	std::ofstream nefile ("NE.dat");
	std::ofstream tgfile ("TG.dat");
	for (int n = 1; n <= max_knots; ++neighbors_i, ++n) {
	    nefile << xyz[ P(n, j)-1 ].coord[0] << "\t"
		   << xyz[ P(n, j)-1 ].coord[1] << "\n";
	}
	tgfile << xyz[id].coord[0] << "\t"
	       << xyz[id].coord[1] << "\n";

	std::cout << "\n+- Knots:  j = " << j << " +--+\n";

	tgfile.close();
	nefile.close();
	    
/* */
    } 

/*!
 *  Set the total number of knots \c N.
 */
    inline bool setN(int n) { N = n; }
/*!
 *  Set the number of interior knots \c NI.
 */
    inline bool setNI(int ni) { NI = ni; }
/*!
 *  Set the number of knots on the boundary \c NB.
 */
    inline bool setNB(int nb) { NB = nb; }
/*!
 *  Get the Knots coordinates.
 *  \param n is of type axis_t, this argument can be X, Y or Z.
 *  \return coord wll contain the n-axis coordinates of all knots.
 */
    inline Vec getKnots(axis_t n) { 
      Vec coord(N); 
// NOTE: Vec begins in 1, and xyz in 0
/****
      for(int i = 1; i <= N; ++i) 
	  coord(i) = xyz[i-1].coord[n];
****/
// NOTE: The above code is not valid for EIGEN!!
      for(int i = 0; i < N; ++i) 
	  coord(i) = xyz[i].coord[n];
      return coord; 
    }
/*!
 *  Get the total number of knots.
 */
    inline int getTotalKnots() { return N; }
/*!
 *  Get the number of interior knots.
 */
    inline int getInteriorKnots() { return NI; }
/*!
 *  Get the number of knots on the boundary.
 */
    inline int getBoundaryKnots() { return NB; }

/*!
 *  Read the knots from a file, \c N must be defined before,
 *  actually it must be defined in the constructor of derived classes.
 */
    bool readFromFile(std::string name) {
	std::ifstream file( name.c_str() );
	if ( !file ) {
	    std::cout << "\n +-----[ ERROR ]-----+"
		      << "\n | Knots::readFromFile(): !Can't open file \""
		      << name << "\" "
		      << "\n +-----+\n";	   
	    exit(1);
	}
	if ( Dim == 1 )
	    for (int i = 0; i < N; ++i) // These are C-arrays 
  	        file >> xyz[i].coord[0];
	if ( Dim == 2 )
	    for (int i = 0; i < N; ++i) // These are C-arrays 
  	        file >> xyz[i].coord[0]
		     >> xyz[i].coord[1];  
	if ( Dim == 3 )
	    for (int i = 0; i < N; ++i) // These are C-arrays
	        file >> xyz[i].coord[0] 
		     >> xyz[i].coord[1] 
		     >> xyz[i].coord[2];
	file.close();
	return 0;
    }
/*!
 *  Write the knots to a file,  \c N must be defined before,
 *  actually it must be defined in the constructor of derived classes.
 */
    bool writeToFile(std::string name) {	
	std::ofstream file( name.c_str() );
	if ( !file ) {
	    std::cout << "\n +-----[ ERROR ]-----+"
		      << "\n | Knots::writeToFile(): !Can't open file \""
		      << name << "\" "
		      << "\n +-----+\n";
	    exit(1);
	}	

	if ( Dim == 1 )
	    for (int i = 0; i < N; ++i) // These are C-arrays
	      file << xyz[i].coord[0] << "\t" << 0.0 << "\n";
	if ( Dim == 2 ) {	   	    
	    for (int i = 0; i < N; ++i) // These are C-arrays
 	        file << xyz[i].coord[0] << "\t"
		     << xyz[i].coord[1] 
//		     << "\t" << 0 // For gnu plotting in 3D
		     << "\n";  
	}
	if ( Dim == 3 )
	    for (int i = 0; i < N; ++i) // These are C-arrays
  	        file << xyz[i].coord[0] << "\t"
		     << xyz[i].coord[1] << "\t"
		     << xyz[i].coord[2] << "\n";
	file.close();
	return 0;	
    }
/*!
 *  Write the knots to files defined by \c I and \c J coordinates of a
 *  Cartesian topology in the MPI context, \c N must be defined before,
 *  actually it must be defined in the constructor of derived classes.
 */
    bool writeToFile(std::string name, int I, int J) {

	std::string istr, jstr;
	std::ostringstream inum, jnum;
	inum.width(1); inum.fill('0'); inum << I; istr = inum.str();
	jnum.width(1); jnum.fill('0'); jnum << J; jstr = jnum.str();
	name += istr + "_" + jstr + ".dat" ;

	std::ofstream file( name.c_str() );
	if ( !file ) {
	    std::cout << "\n +-----[ ERROR ]-----+"
		      << "\n | Knots::writeToFile(): !Can't open file \""
		      << name << "\" "
		      << "\n +-----+\n";
	    exit(1);
	}	
	if ( Dim == 1 )
	    for (int i = 0; i < N; ++i) // These are C-arrays
	      file << xyz[i].coord[0] << "\n";
	if ( Dim == 2 )
	    for (int i = 0; i < N; ++i) // These are C-arrays
 	        file << xyz[i].coord[0] << "\t"
		     << xyz[i].coord[1] 
//		     << "\t" << 0 // For gnu plotting in 3D
		     << "\n";  
	if ( Dim == 3 )
	    for (int i = 0; i < N; ++i) // These are C-arrays
  	        file << xyz[i].coord[0] << "\t"
		     << xyz[i].coord[1] << "\t"
		     << xyz[i].coord[2] << "\n";
	file.close();
	return 0;	
    }

/*!
 *  Print info about the knots to the standard output (screen).
 */
    bool print() {
	std::cout << "\n +-----[ Knots information ]-----+ " 
		  << "\n | RBF : Interior knots = " << NI
		  << "\n | RBF : Boundary knots = " << NB
		  << "\n | RBF : Total knots    = " << N
		  << "\n +-----+";
	return 0;	
    }

}; // Class Knots


#endif // _KNOTS_HPP_

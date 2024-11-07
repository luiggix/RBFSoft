
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

#ifndef _CARTCOMM_HPP_
#define _CARTCOMM_HPP_

#include <iostream>
#include <mpi.h>

/*!
 ***************************************************************************
 * 2D Catersian communicator in the MPI environment. 
 * This class is used to construct a Cartesian topology to be used in simple
 * rectangular-like domains. It is supposed all subdomains that belongs to
 * this topology are rectangles as well. The \c comm variable is an 
 * MPI::Cartcomm communicator, and to be able to access all its functionality 
 * from outside the class in an easy way, this variable is left public. The 
 * numeration of the processors in the topology is as follows:
 * \verbatim    
   ^ y-axis (physical domain)
   |
   +---------+---------+---------+---------+
   |         |         |         |         |   Number of Processors
 ^ |    0    |    1    |    2    |    3    |   
 | |  (0,0)  |  (0,1)  |  (0,2)  |  (0,3)  |   x-axis: NP_J = 4
 | |         |         |         |         |   y-axis: NP_I = 3
 | +---------+---------+---------+---------+   
   |         |         |         |         |   I = 0,...,2; J = 0,...,3
   |    4    |    5    |    6    |    7    |   
 I |  (1,0)  |  (1,1)  |  (1,2)  |  (1,3)  |   NOTE: the I variable runs in
   |         |         |         |         |   y-axis, while J runs in x-axis
 | +---------+---------+---------+---------+   this is contrary to the usual
 | |         |         |         |         |   in numerical methods, and
 | |    8    |    9    |   10    |   11    |   similar to matrix numeration.
 v |  (2,0)  |  (2,1)  |  (2,2)  |  (2,3)  |   
   |         |         |         |         |   
   +---------+---------+---------+---------+---> x-axis (physical domain)
                 <---- J ---->
   \endverbatim   
 * \todo Extend the functionality of CartComm class to 1D and 3D (templates)
 ***************************************************************************
 *  \author  Luis M. de la Cruz [ Tue Mar 18 11:10:23 GMT 2008 ]
 ***************************************************************************
 */
class CartComm
{
public:
    CartComm() { };    
    CartComm(int argc, char **argv, int r);
    
    inline int getNumProc_I() { return num_proc[0]; }
    inline int getNumProc_J() { return num_proc[1]; }
    inline int get_I() { return coords[0]; }    
    inline int get_J() { return coords[1]; }
    inline int getNeighborProc(cartDir_t n) { return neighbors[n]; }
    inline const int* getNeighbors() { return neighbors; }

    void print () {
	std::cout << "\n Proc " << rank 
		  << " (" << coords[1] << "," << coords[1] << ") "
		  << " | UP = " << neighbors[UP]
		  << " | DOWN = " << neighbors[DOWN]
		  << " | LEFT = " << neighbors[LEFT]
		  << " | RIGHT = " << neighbors[RIGHT]
		  << " \n";	   
    }


    MPI::Cartcomm comm; ///< Communicator for a Cartesian topology.

private:
    int num_proc[2];
    int coords[2];
    int neighbors[4];
    bool periods[2];
    int rank;
};

/*!
 *  The Cartesian topology is constructed here. 
 */
CartComm::CartComm(int argc, char **argv, int r)
{
    rank = r;    
    if (argc >= 3) {
        num_proc[0] = std::atoi(argv[2]); // Because the Cartcomm numeration
        num_proc[1] = std::atoi(argv[1]); // is "upside down"
    } else if (rank == 0) {
	std::cerr << "\n\n DD::CartComm : usage: " 
		  << argv[0] << " ns_x ns_y" 
		  << "\n ns_x : number of subdomains in x-axis"
		  << "\n ns_y : number of subdomains in y-axis \n\n";
	std::exit(0);
    }

    // create mpi_cart:  2-dimensional, non-periodic
    periods[0] = false;
    periods[1] = false;
    comm = MPI::COMM_WORLD.Create_cart(2, num_proc, periods, true);
    
    // find your position in the cart
    comm.Get_coords(rank, 2, coords);

    // find your neighbors up, down, left and right,
    // NOTE: if there is not neigbor in some direction, the result is -1.
    comm.Shift(0, 1, neighbors[UP],  neighbors[DOWN]);
    comm.Shift(1, 1, neighbors[LEFT], neighbors[RIGHT]);
}

#endif // _CARTCOMM_HPP_

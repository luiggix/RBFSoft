
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

/*! 
 ***************************************************************************
 * \file Heat2DParallel.cpp
 *  Poisson equation in 2D in parallel Addtive Schwarz.
 *  The boundary problem to be solved is the same as in the example 02Heat2D. 
 *  In this case the problem is solved in parallel using the Additive Schwarz
 *  algorithm. 
 * \image html  topovirtual.png "(a) MPI Cartesian communicator, (b) Partition of the domain" width=2.5cm 
 * \image latex topovirtual.eps "(a) MPI Cartesian communicator, (b) Partition of the domain" width=5cm 
 *  \par Input
 * The \c input file contains the following initial data:
 * - \c hx, \c hy, \c Nx, \c Ny : Lenghts of the rectangle, and number of 
 *                                points in each axis.
 * - \c over_x, \c over_y : Overlapping of the subdomains.
 * - \c rtype, \c ep, \c layer : Type of knots distribution, randomness
 * - \c DDM Max. domain decomposition steps
 * - \c Ngx, \c Ngy, Grid for evaluation.
 * - \c max_knots : neighbors used to construct the ACBF preconditioner.
 * - \c tol_gmres, \c tol_ddm, tolerance for GMRES and DDM algo respectively.
 * \par Output
 * - \c xy_A_B.dat (x,y) coordinates of random points;
 * - \c sol_A_B.dat (x,y,u) Solution for the subdomain (A,B) see above fig
 * - \c err_A_B.dat (x,y,e) Error for the subdomain (A,B) see above fig
 * - \c exa_A_B.dat (x,y,u) Exact solution for the subdomain (A,B) see above fig
 ***************************************************************************
 *  \author  Luis M. de la Cruz [ Tue Feb 19 12:33:23 GMT 2008 ]
 ***************************************************************************
 */
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>

#include "Traits.hpp"
#include "Knots/OVRectangleKnots.hpp"
#include "DD/CartComm.hpp"
#include "RBF/RBF.hpp"
#include "RBF/MQ.hpp"
#include "Solvers/Krylov/diagonalprec.hpp"
#include "Solvers/Krylov/acbfprec.hpp"
#include "Solvers/Krylov/gmres.hpp"
#include "Solvers/Gauss.hpp" 
#include "GNUplot.hpp"

//==========================================================================
//                            FUNCTION PROTOTYPES
//
bool read_data_and_Bcast(CartComm&);

template<typename RBF, typename RBF_2DX, typename RBF_2DY>
void fillMatrices(Mat&, const Vec&, const Vec&, RBF, RBF_2DX, RBF_2DY);
prec_t evalSolution(Vec&, const Vec&, const Vec&, const Vec&, prec_t, 
		    string, string, string, prec_t&);

prec_t exchangeInfo(Vec &u, const Vec&, const Vec&, const Vec&, 
		    prec_t, CartComm&, 
		    OVRectangleKnots<prec_t> &);
prec_t exactSolution(prec_t, prec_t);

//==========================================================================
//                            GLOBAL DATA
//
// Data from "input" file
prec_t hx, hy, over_x, over_y, ep, tol_gmres, tol_ddm;
int Nx, Ny, rtype, max_num_steps, Ngx, Ngy, max_knots, rho;

// Other useful global information
prec_t sx, sy, c;
int N, NI;
int size, rank, I, J, NP_I, NP_J;
//GNUplot plotter;

using namespace std;
int main( int argc, char **argv)
{
    MPI::Init(argc, argv);
    rank = MPI::COMM_WORLD.Get_rank();       
    size = MPI::COMM_WORLD.Get_size();       

// ----- Construction of a Cartesian topology
    CartComm cart(argc, argv, rank);    
    I    = cart.get_I();
    J    = cart.get_J();
    NP_I = cart.getNumProc_I();
    NP_J = cart.getNumProc_J();

// ----- Reading and broadcasting info from the input file
    read_data_and_Bcast(cart);       

// ----- Point generation 
    random_t RT = static_cast<random_t>(rtype);
    OVRectangleKnots<prec_t> rect(hx, Nx, hy, Ny, over_x, over_y, cart, RT);
    rect.setRandomness(ep);
    rect.constructKnots();
    rect.constructOverlapping(cart);    
    rect.writeToFile("xyz_", I, J);
    rect.print();

    Vec x  = rect.getKnots(X);
    Vec y  = rect.getKnots(Y);
    c      = rect.getShapeParameter();
    N      = rect.getTotalKnots();
    NI     = rect.getInteriorKnots();
    Nx     = rect.getNx();            // Number of points in this subdomain.
    Ny     = rect.getNy();
    sx     = rect.getShift(X);        // Shift for this subdomain.
    sy     = rect.getShift(Y);
    over_x = rect.getOverlap(X);      // Overlapping for this subdomain.
    over_y = rect.getOverlap(Y);

// ----- Fill the matrices using MQ and its derivatives.  
    Mat G(N, N);
    Vec u(N), lambda(N), lambda_old(N);
    
    fillMatrices(G, x, y,
		 RBF::MQ<prec_t, 2>(c), 
		 RBF::MQ_2DX<prec_t, 2>(c), 
		 RBF::MQ_2DY<prec_t, 2>(c));

// ----- Boundary conditions.
    int bi, ei;
    int side_x = Nx - 2, side_y = Ny - 2;  

    if (I == NP_I - 1) {            // Just proc. in the bot of the rectangle. 
	bi = NI + 2 * side_y + 1; 
	ei = bi + side_x;
	for(int i = bi; i < ei; ++i) // Edge 4:
	    u(i) = 100.0;
	u(N-3) = 100.00;             // Corner 1:
	u(N-1) = 100.00;             // Corner 3:
    }

// ----- Construction of the Preconditioner.
//    DiagonalPrec<Mat> Jprec(G);    
    ACBFPrec<Mat, OVRectangleKnots<prec_t> > Aprec(G);
    Aprec.construct(rect, max_knots, 0);

// ----- DDM LOOP
    prec_t local_error, error = 1.0, error_lambda = 1.0, sum;
    int niter, ip = 1, ipi = 1, total_iter = 0;

    while ( error > tol_ddm && ip <= max_num_steps ) {
//    while ( error_lambda > tol_gmres && ip <= max_num_steps ) {

//	lambda = Solver::Gauss<prec_t>( G, u );
	niter = gmres (G, lambda, u, Aprec, N-1 , tol_gmres);	 
//	niter = gmres (G, lambda, u, Jprec, N-1 , tol_gmres);	 
//	niter = gmres (G, lambda, u, N-1 , tol_gmres); 

	local_error = exchangeInfo(u, x, y, lambda, c, cart, rect);
	cart.comm.Allreduce(&local_error, &error, 1, MPI::DOUBLE, MPI_MAX);

	sum = 0;
	for(int i = 1; i <= N; ++i) {
	    sum += fabs(lambda_old(i) - lambda(i));
	    lambda_old(i) = lambda(i);
	}
	cart.comm.Allreduce(&sum, &error_lambda, 1, MPI::DOUBLE, MPI_MAX);

	cout << "\n Iter = " << ip << " | Rank = " << rank
	     << " | GMRES = " << niter
  	     << " | E = " << error
	     << " | EL = " << error_lambda
	     << "\n"; 
	ip++;
	total_iter += niter;
    }

    std::string filename, istr, jstr;
    std::ostringstream inum, jnum;
    inum.width(1); inum.fill('0'); inum << I; istr = inum.str();
    jnum.width(1); jnum.fill('0'); jnum << J; jstr = jnum.str();    
    filename = "sol_";
    filename += istr + "_" + jstr + ".dat" ;
    string file2 = "exa_";
    string file3 = "err_";
    file2 += istr + "_" + jstr + ".dat" ;
    file3 += istr + "_" + jstr + ".dat" ;

    prec_t maximum;
    local_error = evalSolution(u, x, y, lambda, c, filename, file2, file3,
			       sum);

    cart.comm.Allreduce(&total_iter, &ipi, 1, MPI::INT, MPI_SUM);

    cart.comm.Allreduce(&local_error, &error, 1, MPI::DOUBLE, MPI_SUM);
    cart.comm.Allreduce(&sum, &maximum, 1, MPI::DOUBLE, MPI_MAX);

    if ( rank == 0)
	cout << "\n | GMRES total iter = " << ipi 
	     << "\t | Max error = " << maximum
	     << "\t | RMS error = " << sqrt( error )
	     << "\n"; 



    MPI::Finalize();
    
    return 0;
}

//==========================================================================
//                            FUNCTION DEFINITIONS
//--------------------------------------------------------------------------
//
bool read_data_and_Bcast(CartComm& cart) 
{
    int    dataBCast_Int[8];
    prec_t dataBCast_Flt[7];

    if (rank == 0) {
	std::cout << "\n\n"
		  << "          MASTER PROCESS rank = " << rank << "\n"
		  << " +----------------------------------------+\n"
		  << " |       TUNA::RBF FOR PDE SOLVING        |\n"
		  << " +----------------------------------------+\n"
		  << "\n";
	
// ----- Reading data from "input" file	
	ifstream input_file ("input");
	input_file >> dataBCast_Flt[0]; // length in x-axis
	input_file >> dataBCast_Flt[1]; // length in y-axis
	input_file >> dataBCast_Int[0]; // Num. of pts in x-axis
	input_file >> dataBCast_Int[1]; // Num. of pts in y-axis
	input_file >> dataBCast_Flt[2]; // Size of overlapping in x-axis
	input_file >> dataBCast_Flt[3]; // Size of overlapping in y-axis
	input_file >> dataBCast_Int[2]; // 0 - Unif, 1 - Rand, 2 - Rand Unif
	input_file >> dataBCast_Flt[4]; // Randomness (0-1.0)
	input_file >> dataBCast_Int[3]; // Maximum num of DDM steps
	input_file >> dataBCast_Int[4]; // Grid for evaluation: x-axis
	input_file >> dataBCast_Int[5]; // Grid for evaluation: y-axis
	input_file >> dataBCast_Int[6]; // Max knots for ACBF precond
	input_file >> dataBCast_Flt[5]; // Tolerance for GMRES
	input_file >> dataBCast_Flt[6]; // Tolerance for DDM
	input_file >> dataBCast_Int[7]; // Rho .....
	input_file.close();
    }
// ----- Broadcast the info to all processors
    cart.comm.Bcast(dataBCast_Int, 8, MPI::INT, 0);    
    cart.comm.Bcast(dataBCast_Flt, 7, MPI::DOUBLE, 0); 

// ----- Using global variables.
    hx        = dataBCast_Flt[0];
    hy        = dataBCast_Flt[1];
    over_x    = dataBCast_Flt[2];
    over_y    = dataBCast_Flt[3];
    ep        = dataBCast_Flt[4];
    tol_gmres = dataBCast_Flt[5];
    tol_ddm   = dataBCast_Flt[6];
    Nx            = dataBCast_Int[0];
    Ny            = dataBCast_Int[1];
    rtype         = dataBCast_Int[2];	
    max_num_steps = dataBCast_Int[3];
    Ngx           = dataBCast_Int[4];
    Ngy           = dataBCast_Int[5];
    max_knots     = dataBCast_Int[6];
    rho           = dataBCast_Int[7];    
    return 0;
}

//--------------------------------------------------------------------------
// Fill the matrix of the system
//
template<typename RBF, typename RBF_2DX, typename RBF_2DY>
void fillMatrices(Mat& G, const Vec& x, const Vec& y,
		  RBF rbf, RBF_2DX rbf_2dx, RBF_2DY rbf_2dy)
{
// ----- Wl matrix
    for(int j = 1; j <= N; ++j) 
	for(int i = 1; i <= NI; ++i) {
	    G(i,j) = rbf_2dx ( x(i), y(i), x(j), y(j) ) + 
		rbf_2dy ( x(i), y(i), x(j), y(j) );
	}
// ----- Wb matrix
    for(int j = 1; j <= N; ++j) 
	for(int i = NI + 1; i <= N; ++i)
	    G(i,j) = rbf ( x(i), y(i), x(j), y(j) ); 
}

//--------------------------------------------------------------------------
// Evaluate the numerical solution on a mesh of Ngx X Ngy
//
prec_t evalSolution(Vec& u, const Vec& x, const Vec& y, const Vec& lambda, 
		    prec_t c, string file_sol, string file_exa, 
		    string file_err, prec_t& maximum)
{
    prec_t sol, exa, diff, sum = 0.0, max = 0.0;

    std::ofstream sol_file (file_sol.c_str());
    std::ofstream exa_file (file_exa.c_str());
/* *
    for(int i = 1; i <= NI; ++i) {
	sol = RBF::eval(x(i), y(i), x, y, lambda, RBF::MQ<prec_t, 2>(c));
	sol_file << x(i) << "\t" << y(i) << "\t" << sol << std::endl;   
//	exa = exactSolution(x(i), y(i) );
//	exa_file << x(i) << "\t" << y(i) << "\t" << exa << std::endl;
	diff = fabs(sol - exa);
	sum += diff;	  	
    }
    
    for(int i = NI + 1; i <= N; ++i)
	sol_file << x(i) << "\t" << y(i) << "\t" << u(i) << std::endl;   


/* */
    std::ofstream err_file (file_err.c_str());

    prec_t xmi,ymi;
    prec_t dx_g = (hx / NP_J) / (Ngx - 1);
    prec_t dy_g = (hy / NP_I) / (Ngy - 1);

    for(int j = 0; j < Ngy ; ++j) {
      for(int i = 0; i < Ngx ; ++i) {
	xmi = dx_g * i + sx + over_x;
	ymi = dy_g * j + sy + over_y;
	
	sol = RBF::eval(xmi, ymi, x, y, lambda, RBF::MQ<prec_t, 2>(c));
	exa = exactSolution(xmi, ymi);

// ----- Dirichlet boundary conditions
	if (J == 0        && i == 0      ) { sol = 0.0; exa = 0.0; }
	if (J == NP_J - 1 && i == Ngx - 1) { sol = 0.0; exa = 0.0; }
	if (I == 0        && j == Ngy - 1) { sol = 0.0; exa = 0.0; }
	if (I == NP_I - 1 && j == 0      ) { sol = 100.0; exa = 100.0; }

	diff = fabs(sol - exa);
	if (diff > max) max = diff;

	sol_file << xmi << "\t" << ymi << "\t" << sol << std::endl;   
	exa_file << xmi << "\t" << ymi << "\t" << exa << std::endl;
	err_file << xmi << "\t" << ymi << "\t" << diff << std::endl;
	
	sum += diff * diff;	  
      }
      sol_file << std::endl;
      exa_file << std::endl;
      err_file << std::endl;
    }

    sum /= (Ngx * Ngy * NP_I * NP_J);

    err_file.close();
    exa_file.close();
    sol_file.close();

    maximum = max;

    return sum;
}

//--------------------------------------------------------------------------
// Evaluate the exact solution in the (x, y) point.
//
prec_t exactSolution(prec_t x, prec_t y)
{
    prec_t c1, s1, s2, s3;
    prec_t sum = 0.0, hxg = 1.0, hyg = 2.0;
    int Nmax = 20;

    for(int n = 1; n <= Nmax; n += 2) {
	c1 = 2.0 / (n * PI);
	s1 = sinh( n * PI * (hyg - y) / hxg );
	s2 = sinh( n * PI * hyg / hxg );
	s3 = sin( n * PI * x / hxg );
	sum += c1 * s1 * s3 / s2;
    }
    return 2.0 * 100.0 * sum ;
}


//--------------------------------------------------------------------------
// Get information in the overlapping from neigbors processors.
//
prec_t exchangeInfo(Vec &u, const Vec& x, const Vec& y, 
		    const Vec& lambda, prec_t c, 
		    CartComm& cart, 
		    OVRectangleKnots<prec_t>& rect)
{
    prec_t sol;
    int ei, bi;

// ----- Get the solution in the overlapping. An interpolation is done
//       using the lambda's calculated for the subdomain, to the overlapping
//       points stored in arrayL, arrayR, arrayU and arrayD from the class
//       OVRectangleKnots. The left, right, up and down vectors defined here
//       will contain the interpolation values at last.
    Vec left = rect.getArrayL();
    ei = left.length();
    for(int i = 1; i < ei; ++i) 
	left(i) = RBF::eval(left(ei), left(i), x, y, 
			lambda, RBF::MQ<prec_t, 2>(c));

    Vec right = rect.getArrayR();
    ei = right.length();
    for(int i = 1; i < ei; ++i) 
	right(i) = RBF::eval(right(ei), right(i), x, y, 
			     lambda,RBF::MQ<prec_t, 2>(c));

    Vec up = rect.getArrayU();
    ei = up.length();
    for(int i = 1; i < ei; ++i)
	up(i) = RBF::eval(up(i), up(ei), x, y, 
			lambda,RBF::MQ<prec_t, 2>(c));

    Vec down = rect.getArrayD();
    ei = down.length();
    for(int i = 1; i < ei; ++i)
	down(i) = RBF::eval(down(i), down(ei), x, y, 
			lambda,RBF::MQ<prec_t, 2>(c));

// Vectors used to receive the information from neigbor subdomains.
    Vec aux_l(Ny), aux_r(Ny), aux_u(Nx), aux_d(Ny);

    MPI::Request req[8];
    MPI::Status  stat[8];

    MPI::Datatype Edge_x = MPI::DOUBLE.Create_contiguous(Nx);
    MPI::Datatype Edge_y = MPI::DOUBLE.Create_contiguous(Ny);
    Edge_x.Commit();
    Edge_y.Commit();

    const int *neighbor = cart.getNeighbors();

    req[0] = cart.comm.Isend(&down(1) , 1, Edge_x, neighbor[DOWN],  0);
    req[1] = cart.comm.Isend(&up(1)   , 1, Edge_x, neighbor[UP],    0);
    req[2] = cart.comm.Isend(&right(1), 1, Edge_y, neighbor[RIGHT], 0);
    req[3] = cart.comm.Isend(&left(1) , 1, Edge_y, neighbor[LEFT],  0);

    req[4] = cart.comm.Irecv(&aux_u(1), 1, Edge_x, neighbor[UP],    0);
    req[5] = cart.comm.Irecv(&aux_d(1), 1, Edge_x, neighbor[DOWN],  0); 
    req[6] = cart.comm.Irecv(&aux_l(1), 1, Edge_y, neighbor[LEFT],  0);
    req[7] = cart.comm.Irecv(&aux_r(1), 1, Edge_y, neighbor[RIGHT], 0);

    MPI::Request::Waitall(8, req, stat);    

    prec_t error = 0.0;
    int side_x = Nx - 2, side_y = Ny - 2;
    if (neighbor[LEFT] != -1) {
	bi = NI + 1;
	ei = bi + side_y;
	for(int i = bi, ii = 2; i < ei; ++i, ++ii) {
	    error += fabs(u(i) - aux_l(ii));
	    u(i) = aux_l(ii);
	}	
	if (neighbor[DOWN] != -1) {  // Corner 1    //
	    error += fabs(u(N-3) - aux_l(1));       // Just calculate the error
	    u(N-3) = aux_l(1);                      // and update the solution
	}                                           // when you are not in a
	if (neighbor[UP] != -1)   {  // Corner 2    // real boundary, where 
	    error += fabs(u(N-2) - aux_l(Ny));      // a Dirichlet condition 
	    u(N-2) = aux_l(Ny);                     // is imposed.
	}                                           //
    }

    if (neighbor[RIGHT] != -1) {
	bi = NI + 1 + side_y;
	ei = bi + side_y;
	for(int i = bi, ii = 2; i < ei; ++i, ++ii) {
	    error += fabs(u(i) - aux_r(ii));
	    u(i) = aux_r(ii);
	}
	if (neighbor[DOWN] != -1) {  // Corner 3    //
	    error += fabs(u(N-1) - aux_r(1));       // Just calculate the error
	    u(N-1) = aux_r(1);                      // and update the solution 
	}                                           // when you are not in a
	if (neighbor[UP] != -1)   {  // Corner 4    // real boundary, where 
	    error += fabs(u(N  ) - aux_r(Ny));      // a Dirichlet condition
	    u(N  ) = aux_r(Ny);                     // is imposed.
	}                                           //
    }

    if (neighbor[DOWN] != -1) {
	bi = NI + 1 + 2 * side_y;
	ei = bi + side_x;
	for(int i = bi, ii = 2; i < ei; ++i, ++ii) {
	    error += fabs(u(i) - aux_d(ii));
	    u(i) = aux_d(ii);
	}
    }

    if (neighbor[UP] != -1) {
	bi = NI + 1 + 2 * side_y + side_x;
	ei = bi + side_x;
	for(int i = bi, ii = 2; i < ei; ++i, ++ii) {
	    error += fabs(u(i) - aux_u(ii));
	    u(i) = aux_u(ii);
	}
    }

    return error;
}

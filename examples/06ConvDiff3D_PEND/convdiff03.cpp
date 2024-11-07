/*! 
 ***************************************************************************
 *
 *  \file convdiff03.cpp
 *  Forced Convection in 3D.
 *  In this example the time-dependent convection-diffusion equation is 
 *  solved. The equations is written as follows:
 *  \f[ \frac{\partial T}{\partial t} + 
 *  u \frac{\partial T}{\partial x} + 
 *  v \frac{\partial T}{\partial y} +
 *  w \frac{\partial T}{\partial z} =
 *  \Gamma\left( \frac{\partial^2 T}{\partial x^2} +
 *               \frac{\partial^2 T}{\partial z^2} +
 *               \frac{\partial^2 T}{\partial y^2} \right)\f]
 *  where \f$(u, v, w)\f$ is a prescribed velocity field that fulfills the 
 *  continuity equation and is given by the next formula:
 *  \f{eqnarray*}{
 *    u(x,y) & = & -A\cos(\pi y) \sin(\pi \lambda x / l_x ) \mbox{and} \\
 *    v(x,y) & = & \frac{A \lambda}{l_x} \sin(\pi y) \cos(\pi \lambda x / l_x)\\
 *    w(x,y) & = & 0.
 *  \f}
 *  The initial condition is \f$T=0\f$ and the boundary conditions are as
 *  shown in the next figure:
 * \image html  Xforcedconv.png "Geometry and boundary conditions" width=5cm 
 * \image latex  Xforcedconv.eps "Geometry and boundary conditions" width=5cm 
 *  \par Compiling and running
 *  Modify the variables FLENS, ATLAS and RBF in the file \c rules.in, 
 *  to the path where these libraries are installed. Then type the next 
 *  \verbatim
    % make
    % ./convdiff03 \endverbatim   
 * The \c input file contains the initial data to setup 
 * the problem.
 * \par Output
 * \c xy_knots.dat coordinates of random points;
 * \c solution.dat (x,y,u) evaluation of the solution in a grid.
 ***************************************************************************
 *  \author  Luis M. de la Cruz Wed [ Dec 19 14:49:56 GMT 2007 ]
 ***************************************************************************
 */

#include <ctime>    
#include <cstdlib>  
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

#include "Traits.hpp"
#include "Knots/BoxKnots.hpp"
#include "Solvers/Gauss.hpp" 
#include "Solvers/Krylov/diagonalprec.hpp"
#include "Solvers/Krylov/identityprec.hpp"
#include "Solvers/Krylov/acbfprec.hpp"
#include "Solvers/Krylov/gmres.hpp"
#include "RBF/RBF.hpp"
#include "RBF/MQ.hpp"

//==========================================================================
//                            FUNCTION PROTOTYPES
//
void initialVelocity(Vec&, Vec&, const Vec&, const Vec&, 
		     prec_t, prec_t);

template<typename RBF, 
	 typename RBF_1DX, typename RBF_1DY, typename RBF_1DZ,
	 typename RBF_2DX, typename RBF_2DY, typename RBF_2DZ>
void fillMatrices(Mat&, const Vec &, const Vec&, const Vec&, 
		  RBF, 
		  RBF_1DX, RBF_1DY, RBF_1DZ, 
		  RBF_2DX, RBF_2DY, RBF_2DZ,
		  Vec&, Vec&, prec_t, prec_t);


bool evalSolution(Vec&, const Vec &, const Vec&, const Vec&, const Vec&, 
		  prec_t);

//==========================================================================
//                            GLOBAL DATA
//
int Nx, Ny, Nz, N, NI, NB;
prec_t hx, hy;

//==========================================================================
//                            MAIN FUNCTION
//
int main( int argc, char * argv[])
{
    timer time;
    prec_t total_time = 0;

    std::cout << "\n\n"
	      << " +----------------------------------------------------+\n"
	      << " |       RADIAL BASIS FUNCTION FOR PDE SOLVING        |\n"
	      << " +----------------------------------------------------+\n"
	      << " | Author: Luis M. de la Cruz                         |\n"
	      << " | Date  : Wed Oct 31 16:55:58 GMT 2007               |\n"
	      << " +----------------------------------------------------+\n"
	      << "\n";

    int max_iter, max_knots;
    prec_t ep, hz, c, dt, A;

    std::ifstream input_file ("input");	
    input_file >> hx;   // length in x-axis
    input_file >> hy;   // length in y-axis
    input_file >> hz;   // length in z-axis
    input_file >> Nx; // Boundary points in x-axis
    input_file >> Ny; // Boundary points in y-axis
    input_file >> Nz; // Boundary points in y-axis
    input_file >> ep;   // Randomness of the interior points (0-0.5)
    input_file >> max_iter; // Number of time steps
    input_file >> dt;   // Time step
    input_file >> max_knots; 
    input_file >> A; // amplitude;
    input_file.close();

 // ----- Point generation 

    std::cout << "\n +-----+  "
	      << "\n | Calculating knots ... ";  
   
    BoxKnots<prec_t> box(hx, Nx, hy, Ny, hz, Nz);
    N  = box.getTotalKnots();
    NI = box.getInteriorKnots();
    NB = box.getBoundaryKnots();
    box.print();
    time.tic();
    box.constructKnots();
    box.writeToFile("xy_knots.dat");

    Vec x = box.getKnots(X);
    Vec y = box.getKnots(Y);
    Vec z = box.getKnots(Z);

    std::cout << "\n | Knots generation elapsed time = " << time.toc()
	      << "\n +-----+  ";

    c = 1.0 / sqrt (1.0 * N);

// ----- Initial condition
    Vec uvel(N), vvel(N), wvel(N); 
    prec_t rolls = 1.0;
    initialVelocity(uvel, vvel, x, y, A, rolls);

// ----- Fill the matrices using MQ and its derivatives.      

    std::cout << "\n +-----+  ";
    std::cout << "\n | Filling the linear system ...  ";
    
    Mat G(N, N);
    Vec u(N), lambda(N);
    
    time.tic();
    fillMatrices(G, x, y, z,
		 RBF::MQ<prec_t, 3>(c), 
		 RBF::MQ_1DX<prec_t, 3>(c), 
		 RBF::MQ_1DY<prec_t, 3>(c), 
		 RBF::MQ_1DZ<prec_t, 3>(c), 
		 RBF::MQ_2DX<prec_t, 3>(c), 
		 RBF::MQ_2DY<prec_t, 3>(c), 
		 RBF::MQ_2DZ<prec_t, 3>(c), 
		 uvel, vvel, dt, 1.0);

    std::cout << "\n | Elapsed time = " << time.toc();

    time.tic();
    std::cout << "\n | Construction of the preconditioner....";

//    ACBFPrec<Mat, BoxKnots<prec_t> > precond(G);
//    precond.construct(box, max_knots);
//    IdentityPrec<Mat> precond(G);
    DiagonalPrec<Mat> precond(G);
    std::cout << "\n | Elapsed time = " << time.toc();
    std::cout << "\n +-----+  ";

// ----- Boundary conditions.

    int pausa;
    int bi, ei;

// Dirichlet conditions in vertical walls.
    int Nxy = (Nx - 2) * (Ny - 2);
    int Nxz = (Nx - 2) * (Nz - 2);
    int Nyz = (Ny - 2) * (Nz - 2);

// Face 1: Dirichlet Temp = 0.5
    bi = NI + 1;
    ei = bi + Nyz;
    for(int i = bi; i < ei; ++i) {
	u(i) = 0.5;
    }

// Face 2: Dirichlet Temp = -0.5
    bi = ei ;
    ei = bi + Nyz;
    for(int i = bi; i < ei; ++i) {
	u(i) = -0.5;
    }

// Edges 1 and 2: Dirichlet Temp = 0.5
    bi = ei + 2 * Nxz +  2 * Nxy;
    ei = bi + 2 * (Nz - 2);
    for(int i = bi; i < ei; ++i) {
	u(i) = 0.5;
    }

// Edges 3 and 4: Dirichlet Temp = -0.5
    bi = ei;
    ei = bi + 2 * (Nx - 2);
    for(int i = bi; i < ei; ++i) {
	u(i) = -0.5;
    }

// Edges 5 and 6: Dirichlet Temp = 0.5
    bi = ei;
    ei = bi + 2 * (Ny - 2);
    for(int i = bi; i < ei; ++i) {
	u(i) = 0.5;
    }

// Edges 7 and 8: Dirichlet Temp = -0.5
    bi = ei;
    ei = bi + 2 * (Ny - 2);
    for(int i = bi; i < ei; ++i) {
	u(i) = -0.5;
    }

// Corners 1,2,3 and 4: Dirichlet Temp = 0.5
    for (int i = N - 8, count = 1; count <= 4; ++i, ++count)
	u(i) = 0.5;

// Corners 5,6,7 and 8: Dirichlet Temp = -0.5
    for (int i = N, count = 1; count <= 4; --i, ++count)
	u(i) = -0.5;

// ----- Main Loop

    int niter;
    prec_t tol = 1.0e-8;
    for(int t = 1; t <= max_iter; ++t) 
    {			
	time.tic();
/* *
// Gauss Solver

	lambda = Solver::Gauss<prec_t>( G, u );
	std::cout << "\n | " << t ;
/* */
// GMRES Solver 

	niter = gmres (G, lambda, u, precond, N-1 , tol);    
	std::cout << "\n | " << t <<"\t GMRES Iterations : " << niter;
/* */
	total_time += time.toc();

	evalSolution(u, x, y, z, lambda, c);
//	std::cin >> pausa;
    }
    std::cout << "\n Total time ( Ax = b ) = " << total_time;

    std::cout << "\n +-----+  ";
    std::cout << "\n + Happy finish :-) ";
    std::cout << "\n +-----+ \n\n";    

    return 0;
}

//==========================================================================
//                            FUNCTION DEFINITIONS
//--------------------------------------------------------------------------
// Calculate the solenoidal velocity
// 
void initialVelocity(Vec& u, Vec& v, const Vec& x, const Vec& y,
		     prec_t A, prec_t rolls)
{
    prec_t A2 = A * rolls / hx;
    std::ofstream uvel_file ("uvel.dat");
    std::ofstream vvel_file ("vvel.dat");
    int cnt = 1;
    for(int k = 2; k < Nz; ++k) {
	for(int j = 2; j < Ny; ++j) {
	    for(int i = 2; i < Nx; ++i, ++cnt) {
		u(cnt) = -A * cos(PI * y(cnt)) * sin(PI * rolls * x(cnt) /hx);
		v(cnt) = A2 * sin(PI * y(cnt)) * cos(PI * rolls * x(cnt) /hx);
		uvel_file << x(cnt) << "\t" << y(cnt) << "\t" << u(cnt) <<"\n";
		vvel_file << x(cnt) << "\t" << y(cnt) << "\t" << v(cnt) <<"\n";
	    }	    
	    uvel_file << "\n";
	    vvel_file << "\n";	
	}
	uvel_file << "\n";
	vvel_file << "\n";
    }
    uvel_file.close();
    vvel_file.close();
}
//--------------------------------------------------------------------------
// Fill the matrix of the system
//
template<typename RBF, 
	 typename RBF_1DX, typename RBF_1DY, typename RBF_1DZ,
	 typename RBF_2DX, typename RBF_2DY, typename RBF_2DZ>
void fillMatrices(Mat& G, const Vec& x, const Vec& y, const Vec& z,
		  RBF rbf, 
		  RBF_1DX rbf_1dx,
		  RBF_1DY rbf_1dy,
		  RBF_1DZ rbf_1dz,
		  RBF_2DX rbf_2dx, 
		  RBF_2DY rbf_2dy, 
		  RBF_2DZ rbf_2dz, 
		  Vec& uvel, Vec& vvel, prec_t dt, prec_t D)
{
    int N = G.numRows(); //NI + NB;

// ----- Wl matrix

    for(int j = 1; j <= N; ++j) 
	for(int i = 1; i <= NI; ++i) {
	    G(i,j) = rbf ( x(i), y(i), z(i), x(j), y(j), z(j) ) +
		dt * ( uvel(i) * rbf_1dx (x(i), y(i), z(i), x(j), y(j), z(j)) +
		       vvel(i) * rbf_1dy (x(i), y(i), z(i), x(j), y(j), z(j)) -
		       D * ( rbf_2dx ( x(i), y(i), z(i), x(j), y(j), z(j) ) + 
			     rbf_2dy ( x(i), y(i), z(i), x(j), y(j), z(j) ) +
			     rbf_2dz ( x(i), y(i), z(i), x(j), y(j), z(j) )
			   ) 
		    );
	}

    int Nxy = (Nx - 2) * (Ny - 2);
    int Nxz = (Nx - 2) * (Nz - 2);
    int Nyz = (Ny - 2) * (Nz - 2);

// ----- Wb matrix
    int bi, ei;
    for(int j = 1; j <= N; ++j) {

// Face 1 and 2 : Dirichlet
	bi = NI + 1;
	ei = bi + 2 * Nyz;
	for(int i = bi; i < ei; ++i) {
	    G(i,j) = rbf ( x(i), y(i), z(i), x(j), y(j), z(j) );   
	}

// Face 3 and 4 : Neumann
	bi = ei;      
	ei = bi + 2 * Nxz;    
	for(int i = bi; i < ei; ++i) {
	    G(i,j) = rbf_1dy ( x(i), y(i), z(i), x(j), y(j), z(j) ); 
	}

// Face 5 and 6 : Neumann
	bi = ei;      
	ei = bi + 2 * Nxy;    
	for(int i = bi; i < ei; ++i) {
	    G(i,j) = rbf_1dz ( x(i), y(i), z(i), x(j), y(j), z(j) ); 
	}

// Edges 1,2,3 and 4 : Dirichlet as in faces 1 and 2
	bi = ei;
	ei = bi + 4 * (Nz - 2);
	for(int i = bi; i < ei; ++i) {
	    G(i,j) = rbf ( x(i), y(i), z(i), x(j), y(j), z(j) ); 
	}

// Edges 5,6,7 and 8 : Dirichlet as in faces 1 and 2
	bi = ei;
	ei = bi + 4 * (Ny - 2); 
	for(int i = bi; i < ei; ++i) {
	    G(i,j) = rbf ( x(i), y(i), z(i), x(j), y(j), z(j) ); 
	}

// Edges 9,10,11 and 12 : Neumann as in faces 5 and 6
	bi = ei;
	ei = bi + 4 * (Nx - 2);
	for(int i = bi; i < ei; ++i) {
	    G(i,j) = rbf_1dy ( x(i), y(i), z(i), x(j), y(j), z(j) ); 
	}

// Corners : Dirichlet as in faces 3 and 4
	bi = ei;
	ei = bi + 8;
	for(int i = bi; i < ei; ++i) {
	    G(i,j) = rbf ( x(i), y(i), z(i), x(j), y(j), z(j) );
	}
    }
}
//--------------------------------------------------------------------------
// Evaluate the numerical solution on a mesh of Nx by Ny
//

bool evalSolution(Vec& u, const Vec& x, const Vec& y, const Vec& z, 
		  const Vec& lambda, prec_t c)
{
    int N = Nx * Ny * Nz;
    int NI = (Nx - 2) * (Ny - 2) * (Nz - 2);

    for(int i = 1; i <= N; ++i) {	  
	u(i) = RBF::eval(x(i), y(i), z(i), x, y, z, 
			 lambda, RBF::MQ<prec_t,3>(c),0);
    }

    int Nxy = (Nx - 2) * (Ny - 2);
    int Nxz = (Nx - 2) * (Nz - 2);
    int Nyz = (Ny - 2) * (Nz - 2);
    int bi, ei;   

/*
// Face 3: 
    bi = NI + 2 * Nyz + 1;
    ei = bi + Nxz;
    for(int i = 1; i <= Nx - 2; ++i)
	u(i) = ;
    
    bi = NI + Nx + Ny;
    ei = bi + Nx - 3;
    for(int i = bi, ii = 0; i <= ei; ++i, ++i)
	u(i) = u(NI - ii);
    
/* */

    std::ofstream pos_sol ("solution.dat");
    int entero = (Nz / 2) - 1;
    bi = entero * Nxy + 1;
    ei = bi + Nxy;
    for(int i = bi; i < ei; ++i) {
	pos_sol << x(i) <<" "<< y(i) <<"  " << u(i) << std::endl;
	if ( !(i % (Nx - 2)) ) pos_sol << std::endl;
    }

//    std::cout << "\n | RMS error = " << sqrt (sum /(Nx * Ny));
//    std::cout << "\n | Max error = " << max;

    return 0;
} 


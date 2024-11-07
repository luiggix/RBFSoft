/*! 
 ***************************************************************************
 *  \file    convdiff01.cpp
 *  Non-stat adv-diff equation in 1D.
 *  In this example the time-dependent advection-diffusion equation in 1D 
 *  is solved. The equation is written as follows:
 *  \f[ \frac{\partial f}{\partial t} + u \frac{\partial f}{\partial x} 
 *  - D \frac{\partial^2 f}{\partial x^2} = 0 \f]
 *  where \f$u\f$ represents a convective velocity. The initial and 
 *  boundary conditions are:
 *  \f{eqnarray*}{
 *  f(x,0) & = & 0 \textrm{ for } 0 \le x < \infty, \\
 *  f(0,t) & = & 1 \textrm{ for } t \ge 0,  \\
 *  f(L,t) & = & 0 \textrm{ for } t > 0, L \rightarrow \infty,
 *  \f}
 * \image html  soladvdiff1D.png "MQ-RBF solution for D = 0.0001 and 50 points" width=5cm 
 * \image latex soladvdiff1D.eps "MQ-RBF solution for D = 0.0001 and 50 points" width=5cm 
 *  \par Input
 * The \c input file contains the initial data to setup the problem:
 * - \c h lenght of the domain
 * - \c N number of points
 * - \c max_iter max time iterations
 * - \c dt time step
 * - \c D Peclet number
 * - \c max_knots neighbor knots for the ACBF precond.
 * \par Output
 * - \c xy_knots.dat x-coordinates of the points.
 * - \c solution.dat evaluation of the final solution.
 * - \c error.dat error against the exact solution.
 * - \c exact.dat exact solution.
 ***************************************************************************
 * \author  Luis M. de la Cruz [ Thu Sep  6 14:35:41 BST 2007 ]
 ***************************************************************************
 */

#include <cmath>
#include <iostream>
#include <fstream>
#include "Traits.hpp"
#include "Knots/LineKnots.hpp"
#include "Solvers/Gauss.hpp" 
#include "Solvers/Krylov/identityprec.hpp"
#include "Solvers/Krylov/diagonalprec.hpp"
#include "Solvers/Krylov/acbfprec.hpp"
#include "Solvers/Krylov/gmres.hpp"
#include "RBF/RBF.hpp"
#include "RBF/MQ.hpp"


#include "GNUplot.hpp"

//==========================================================================
//                            FUNCTION PROTOTYPES
//
template<typename RBF, typename RBF_1DX, typename RBF_2DX>
void fillMatrices(Mat&, const Vec&, RBF, RBF_1DX, RBF_2DX, 
		  prec_t, prec_t, prec_t);
prec_t evalSolution(Vec&, const Vec&, const Vec&, prec_t, 
		  prec_t, prec_t, prec_t);

//==========================================================================
//                            GLOBAL DATA
//
int N, NI, NB; // Number of points involved in the calculation
double errorL1, errorL2, errorLmax;
double nL1, nL2, nLmax;
double ReL1, ReL2, ReLmax;

//==========================================================================
//                            MAIN FUNCTION
//
int main( int argc, char * argv[])
{
    timer time;

    std::cout << "\n\n"
	      << " +----------------------------------------+\n"
	      << " |       TUNA::RBF FOR PDE SOLVING        |\n"
	      << " +----------------------------------------+\n"
	      << "\n";

    int max_iter, max_knots;
    prec_t h, c, dt, D, vel = 1.0;

    std::ifstream input_file ("input");
    input_file >> h   // length in x-axis
	       >> N   // Points in x-axis for the grid (to evaluate sol.)
	       >> D   // Peclet
	       >> dt   // Time step
	       >> max_iter // Number of time steps
	       >> max_knots // for ACBF precond
               >> c;
    input_file.close();

// ----- Point generation 


    LineKnots<prec_t> line(h, N);
    line.constructKnots();    
    N  = line.getTotalKnots();
    NI = line.getInteriorKnots();
    NB = line.getBoundaryKnots();
    line.print();
    line.constructKnots();
    line.writeToFile("xy_knots.dat");

    Vec x = line.getKnots(X);

//    c = 1.0 / sqrt(1.0 * N);

    prec_t dx = h / (N - 1.0);
    for(int i = 1; i <= N-2; ++i) 
	x(i) = i * dx;
    x(NI + 1) = 0.0;
    x(N) = h;

// ----- Fill the matrices using MQ and its derivatives.  

    Mat G(N, N);
    Vec u(N), lambda(N);

    fillMatrices(G, x,
		 RBF::MQ<prec_t, 1>(c), 
		 RBF::MQ_1DX<prec_t, 1>(c),
		 RBF::MQ_2DX<prec_t, 1>(c), 
		 dt, vel, D);

    double time_precond;
    time.tic();
    //    IdentityPrec<Mat> precond(G);
//    DiagonalPrec<Mat> precond(G);
    ACBFPrec<Mat, LineKnots<prec_t> > precond(G);
    precond.construct(line, max_knots, 2);
    time_precond = time.toc();
    

// ----- Boundary conditions.
    u(NI + 1) = 1.0;
    u(N) = 0.0;

    prec_t tol = 1.0e-6;
    prec_t total_time = 0;

#ifdef WITH_GNUPLOT
    GNUplot plotter;
    plotter("set grid");
#endif

    prec_t error;
    std::ofstream err_file ("error.dat");

    int niter = 0;
    double dniter = 0;

    int Restart = N - 1;
//    Restart /= 2;
// ----- Main loop
    for(int t = 1; t < max_iter; ++t) 
    {

// ----- Solve the Linear System  G * lambda = u       
	time.tic();

/* *
// Gauss Solver

	lambda = Solver::Gauss<prec_t>( G, u );
	total_time += time.toc();
/* */
// GMRES Solver
	niter = gmres (G, lambda, u, precond, Restart, tol); 
	total_time += time.toc();

	dniter += niter;

/* *
	std::cout << "\n |  it = "<< t << "\t GMRES Iterations : " << niter;
/* */


	error = evalSolution(u, x, lambda, c, dt * t, vel, D);
//	err_file << t << "\t" << error << "\n";

/* *
	cout << "\n Error L1 = " << errorL1
	     << "\t Error L2 = " << errorL2 
	     << "\n Error R L1 = " << ReL1
	     << "\t Error R L2 = " << ReL2
	     << "\t Error R Lmax = " << ReLmax 
	     << "\t Error Lmax = " << errorLmax 
	     << "\t Total time = " << total_time 
	     << "\n Time preconditioner = " << time_precond
	     << endl;
/* */
#ifdef WITH_GNUPLOT
	plotter("p [0:2.5][-0.1:1.1] \"solution.dat\" w lp,\"exact.dat\" w l");
#endif

    }

    dniter /= max_iter;
    cout << "\n " <<  c << " "
	 << errorL1 << " " 
	 << errorL2 << " "
	 << ReL1 << " "
	 << ReL2 << " "
	 << ReLmax << " "
	 << errorLmax << " "
	 << dniter << " "
	 << total_time << endl;

#ifdef WITH_GNUPLOT
    int pausa; 
    std::cout << "\n\n >---> Press any key and then <enter> to finish " ;
    cin >> pausa;
#endif

    std::cout << "\n +-----+  ";
    std::cout << "\n + Happy finish :-) ";
    std::cout << "\n +-----+ \n\n";    

    return 0;
}

//==========================================================================
//                            FUNCTION DEFINITIONS
//--------------------------------------------------------------------------
// Fill the matrix of the system
//
template<typename RBF, typename RBF_1DX, typename RBF_2DX>
void fillMatrices(Mat& G, const Vec& x, 
		  RBF rbf, 
		  RBF_1DX rbf_1dx, 
		  RBF_2DX rbf_2dx, 
		  prec_t dt, prec_t vel, prec_t D)
{
//    int N = G.numRows(); //NI + NB;

// ----- Wl matrix
    for(int j = 1; j <= N; ++j) 
	for(int i = 1; i <= NI; ++i)
	    G(i,j) = rbf ( x(i), x(j) ) +
		dt * ( vel * rbf_1dx ( x(i), x(j) ) -
		       D   * rbf_2dx ( x(i), x(j) ) ) ; 


// ----- Wb matrix
    for(int j = 1; j <= N; ++j) 
	for(int i = NI+1; i <= N; ++i)
	  G(i,j) = rbf ( x(i), x(j) );   
}

//--------------------------------------------------------------------------
// Evaluate the numerical solution
//
prec_t evalSolution(Vec& u, const Vec& x, const Vec& lambda, 
		  prec_t c, prec_t t, prec_t vel, prec_t D)
{
    int N = u.length();
    std::ofstream sol_file ("solution.dat"), exa_file ("exact.dat");
    prec_t ua;
    prec_t diff = 0.0, max = 0.0;

    nL1 = 0.0, nL2 = 0.0, nLmax = 0.0;
    errorL1 = errorL2 = errorLmax = 0.0;

    sol_file << x(NI+1) <<"\t" << u(NI+1) << std::endl;
    exa_file << x(NI+1) <<"\t" << u(NI+1) << std::endl;
    for(int i = 1; i <= NI; ++i) {
	u(i) = RBF::eval(x(i), x, lambda, RBF::MQ<prec_t, 1>(c));
	sol_file << x(i) <<"\t" << u(i) << std::endl; 

// Exact solution
	ua = 0.5 * ( erfc( 0.5 * (x(i) - vel * t) / sqrt(D * t) ) + 
		     exp(vel * x(i)) * exp(-D) * 
		     erfc( 0.5 * (x(i) + vel * t) / sqrt(D * t)) );
	exa_file << x(i) <<"\t" << ua << std::endl; 

	diff = fabs( u(i) - ua );
	errorL1 += diff;
	errorL2 += diff * diff;
	if ( max < diff ) errorLmax = diff;

	nL1 += fabs( ua );
	nL2 += ua * ua;
	if ( nLmax < nL1 ) nLmax = nL1;

    }
    sol_file << x(N) <<"\t" << u(N) << std::endl;
    exa_file << x(N) <<"\t" << u(N) << std::endl;

    sol_file.close();
    exa_file.close();

    ReL1 = errorL1 / nL1;
    ReL2 = sqrt (errorL2) / sqrt(nL2);
    ReLmax = errorLmax / nLmax;

    errorL1 /= N;
    errorL2 = sqrt( errorL2 / N );
    

    return errorL2;
} 

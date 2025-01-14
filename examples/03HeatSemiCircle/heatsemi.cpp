/*! 
 ***************************************************************************
 *  \file    heatsemi.cpp
 *  Laplace equation in a semicircle domain.
 *  The PDE to solve is: 
 *  \f[ \frac{\partial^2 T}{\partial x^2}+\frac{\partial^2 T}{\partial y^2}=0 \f]
 *  and the boundary conditions are (see the figure):
 * - \f$ T = \sin(\theta) / R_{0A} \f$  for segment AC;
 * - \f$ T = \sin(\theta) / R_{0B} \f$  for segment BD;
 * - \f$ T = 0 \f$ for segment AB;
 * - \f$ T = 1 / r \f$ for segment CD.
 * \image html  heatsemi.png "Geometry, boundary conditions and solution" width=5cm 
 * \image latex heatsemi.eps "Geometry, boundary conditions and solution" width=5cm 
 * 
 * \par Input
 * The \c input file contains the initial data to setup the problem.
 * - \c hx, \c hy, \c Nx, \c Ny : Lenghts of the rectangle, and number of 
 *                                points in each axis.
 * - \c rtype, \c ep, \c layer : Type of knots distribution, randomness, layer
 *                                around the domain.
 * - \c c : Shape parameter for MQ-RBF kernel (\f$c < 0\f$ implies \f$ c = 1/\sqrt{N}\f$).
 * - \c fsol : choose a method to solve the system.
 * - \c max_knots : neighbors used to construct the ACBF preconditioner. 
 * .
 * \par Output
 * - \c xy_knots.dat coordinates of random points;
 * - \c solrand.dat (x,y,u) evaluation of the solution in the collocation pts.
 * - \c errrand.dat (x,y,|e-u|) error distribution in the collocation pts.
 * - \c exarand.dat (x,y,e) exact solution in the collocation pts.
 * - \c solgrid.dat (x,y,u) evaluation of the solution in a grid;
 * - \c errgrid.dat (x,y,|e - u|) error distribution in a grid;
 * - \c exagrid.dat (x,y,e) exact solution in a grid;
 * .
 ***************************************************************************
 *  \author  Luis M. de la Cruz [ Thu Sep  6 14:35:41 BST 2007 ]
 ***************************************************************************
 */

#include <iostream>
#include <fstream>
#include <string>

#include "Traits.hpp"
#include "Knots/SemiCircleKnots.hpp"
#include "Solvers/Gauss.hpp" 
#include "Solvers/Krylov/diagonalprec.hpp"
#include "Solvers/Krylov/identityprec.hpp"
#include "Solvers/Krylov/acbfprec.hpp"
#include "Solvers/Krylov/gmres.hpp"
#include "RBF/RBF.hpp"
#include "RBF/MQ.hpp"
#include "GNUplot.hpp"

//==========================================================================
//                            FUNCTION PROTOTYPES
//
template<typename RBF, typename RBF_2DX, typename RBF_2DY>
void fillMatrices(Mat&, const Vec&, const Vec&, RBF, RBF_2DX, RBF_2DY);
bool evalSolution(const Vec&, const Vec&, const Vec&, const Vec&, prec_t);

//==========================================================================
//                            GLOBAL DATA
//
int NX, NY, N, NI, NB, Nt, Nr;
prec_t r1, r2, t1, t2, ep;

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

    int max_knots, fsol, rtype;

// ----- Reading data from "input" file

    std::ifstream input_file ("input");
    input_file >> r1 >> r2 >> t1 >> t2 >> Nr >> Nt
	       >> rtype >> ep >> fsol >> max_knots >> NX >> NY;
    input_file.close();
    random_t RT = static_cast<random_t>(rtype);

    prec_t c;

// ----- Point generation 
    std::cout << "\n +-----+  "
	      << "\n | Calculating knots ... ";   

    SemiCircleKnots<prec_t> semi(r1, r2, Nr, t1, t2, Nt, RT);
    N  = semi.getTotalKnots();
    NI = semi.getInteriorKnots();
    NB = semi.getBoundaryKnots();
    semi.setRandomness(ep);
    semi.print();
    time.tic();
    semi.constructKnots();
    Vec x = semi.getKnots(X);
    Vec y = semi.getKnots(Y);

    std::cout << "\n | Knots generation elapsed time = " << time.toc()
	      << "\n +-----+  "; 

    semi.writeToFile("xy_knots.dat");

    c = 1.0 / sqrt (1.0 * N);

// ----- Fill the matrices using MQ and its derivatives.  
    std::cout << "\n +-----+  ";
    std::cout << "\n | Filling the linear system ...  ";

    Mat G(N, N);
    Vec u(N), lambda(N);

    time.tic();
    fillMatrices(G, x, y,
		 RBF::MQ<prec_t, 2>(c), 
		 RBF::MQ_2DX<prec_t, 2>(c), 
		 RBF::MQ_2DY<prec_t, 2>(c));

// ----- Boundary conditions.

// Edge 1: Dirichlet: sin(theta) / r1
    int bi, ei;
    int side_x = Nr - 2, side_y = Nt - 2;    
    bi = NI + 1;
    ei = bi + side_y;
    for(int i = bi; i < ei; ++i) 
	u(i) = y(i)/ ( r1 * r1 );  // sin(theta) / r1
    
// Edge 2: Dirichlet: sin(theta) / r2
    bi = ei;
    ei = bi + side_y;
    for(int i = bi; i < ei; ++i) 
	u(i) = y(i)/ ( r2 * r2 );  // sin(theta) / r2

// Edge 3: Dirichlet: T = 0
// Nothing to do

// Edge 4: Dirichlet: 1 / r
    bi = ei + side_x;
    ei = bi + side_x;
    for(int i = bi; i < ei; ++i) 
	u(i) = 1.0 / y(i);  // 1 / r

// Corner 2: Dirichlet T = 1 / r
    u(N-2) = 1.0 / r1;
// Corner 4: Dirichlet T = 1 / r
    u(N) = 1.0 / r2;

    std::cout << "\n | Elapsed time = " << time.toc()
	      << "\n +-----+  ";   

// ----- Solve the Linear System  G * lambda = u
    std::cout << "\n +-----+  "
	      << "\n | Solving the linear system : ";

    prec_t tol = 1.0e-8;
    int niter;

    IdentityPrec<Mat> Iprec(G);
    DiagonalPrec<Mat> Jprec(G);
    ACBFPrec<Mat, SemiCircleKnots<prec_t> > Aprec(G);

    switch (fsol) 
    {
	case 1: // Gauss Solver
	    time.tic();
	    std::cout << " Gauss Elimination ....";	    
	    lambda = Solver::Gauss<prec_t>( G, u );
	    break;
	case 2: // GMRES No preconditioned
	    time.tic();
	    std::cout << " GMRES NO PREC ....";
	    niter = gmres (G, lambda, u, Iprec, N-1 , tol);    
	    break;
	case 3: // GMRES Jacobi preconditioner
	    time.tic();
	    std::cout << " GMRES JACOBI ....";
	    niter = gmres (G, lambda, u, Jprec, N-1 , tol); 	    
	    break;
	case 4: // GMRES ACBF preconditioner
	    time.tic();
	    std::cout << " GMRES ACBF ....";
	    Aprec.construct(semi, max_knots, 8, true);
	    std::cout << "\n | Elapsed time - ACBF precond construction: " 
		      << time.toc();
	    time.tic();
	    niter = gmres (G, lambda, u, Aprec, N-1 , tol); 	    
	    break;
	default:
	    break;
    }
    if (fsol != 1) {
	std::cout << "\n | GMRES Iterations : " << niter;
    }
    std::cout << "\n | Elapsed time for solving Ax = b : "
	      << time.toc() << " seconds "
	      << "\n +-----+  ";

// ----- Evaluation of the solution.

    std::cout << "\n +-----+  ";
    std::cout << "\n | Evaluating the solution on the grid ... ";
    time.tic();
    evalSolution(x, y, lambda, u, c);
    std::cout << "\n | Elapsed time = " << time.toc() ;   
    std::cout << "\n +-----+  ";

#ifdef WITH_GNUPLOT
    int pausa;
    GNUplot plotter;
    plotter("set view 60,15");
    plotter("set contour");
    plotter("set cntrparam level 10");
    plotter("splot \"solrand.dat\" w p, \"solgrid.dat\" w l");
    std::cout << "\n\n >---> Press any key and then <enter> to finish " ;
    std::cin >> pausa;
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
template<typename RBF, typename RBF_2DX, typename RBF_2DY>
void fillMatrices(Mat& G, const Vec& x, const Vec& y,
		  RBF rbf, 
		  RBF_2DX rbf_2dx, 
		  RBF_2DY rbf_2dy)
{
// ----- Wl matrix
    for(int j = 1; j <= N; ++j) 
	for(int i = 1; i <= NI; ++i) {
	  G(i,j) = rbf_2dx ( x(i), y(i), x(j), y(j) ) + 
	           rbf_2dy ( x(i), y(i), x(j), y(j) );
	}
// ----- Wb matrix
    for(int j = 1; j <= N; ++j) 
	for(int i = NI+1; i <= N; ++i)
	  G(i,j) = rbf ( x(i), y(i), x(j), y(j) );   

}

//--------------------------------------------------------------------------
// Evaluate the numerical solution on a mesh of Nx by Ny
//
bool evalSolution(const Vec& x, const Vec& y, const Vec& lambda, 
		  const Vec& u, prec_t c)
{

    std::ofstream pos_sol ("solrand.dat");
    std::ofstream exa_sol ("exarand.dat");
    std::ofstream err_sol ("errrand.dat");
    prec_t diff, sol, exa, sum = 0.0, max = 0.0;

    for(int i = 1; i <= NI; ++i) {	  
	sol = RBF::eval(x(i), y(i), x, y, lambda, RBF::MQ<prec_t,2>(c),0);
	pos_sol << x(i) <<" "<< y(i) <<"  "<< sol << std::endl;   
	if ( !(i % (Nr-2)) ) pos_sol << "\n";

	exa = y(i) / ( x(i)*x(i) + y(i)*y(i) ); // sin(theta)/r
	exa_sol << x(i) <<" "<< y(i) <<"  " << exa << std::endl;   
	if ( !(i % (Nr-2)) ) exa_sol << "\n";

	diff = fabs ( sol - exa );
	err_sol << x(i) <<" "<< y(i) <<"  " << diff << std::endl;   
	if ( !(i % (Nr-2)) ) err_sol << "\n";

	if (diff > max) max = diff;
	sum += (diff * diff);
    }   
    std::ofstream pos_solb ("solutionB.dat");    
    for(int i = NI+1; i <= N; ++i) {
	pos_solb << x(i) <<" "<< y(i) <<"  "<< u(i) << std::endl;   
    }
    pos_sol.close();
    pos_solb.close();    
    exa_sol.close();
    err_sol.close();
    std::cout << "\n | RMS error (rand) = " << sqrt ( sum /(NI) );
    std::cout << "\n | Max error (rand) = " << max;    
    
    sum = 0; max = 0;
    int ii = 1;
    prec_t solg;
    std::ofstream pos_grd ("solgrid.dat");
    std::ofstream exa_grd ("exagrid.dat");
    std::ofstream err_grd ("errgrid.dat");

    SemiCircleKnots<double> aux(r1, r2, NX, t1, t2, NY);
    aux.gridKnots();
    Vec xx = aux.getKnots(X);
    Vec yy = aux.getKnots(Y);
    for(int j = 1; j <= NY; ++j) {
	for(int i = 1; i <= NX; ++i, ++ii) {
	    solg = RBF::eval(xx(ii), yy(ii), 
			     x, y, lambda, RBF::MQ<prec_t,2>(c),0);
	    pos_grd << xx(ii) <<" "<< yy(ii) <<"  "<< solg << std::endl;   
	    
	    sol = yy(ii) / ( xx(ii)*xx(ii) + yy(ii)*yy(ii) ); // sin(theta)/r
	    exa_grd << xx(ii) <<" "<< yy(ii) <<"  " << sol << std::endl;   
	    
	    diff = fabs ( sol - solg );
	    err_grd << xx(ii) <<" "<< yy(ii) <<"  " << diff << std::endl;   
	    
	    if (diff > max) max = diff;
	    sum += (diff * diff);	    
	}
	pos_grd << "\n";
	exa_grd << "\n";
	err_grd << "\n";
    }
    pos_grd.close();
    exa_grd.close(); 
    err_grd.close(); 

    std::cout << "\n | RMS error (grid) = " << sqrt ( sum /(NX*NY) );
    std::cout << "\n | Max error (grid) = " << max;    
   
    return 0;
} 



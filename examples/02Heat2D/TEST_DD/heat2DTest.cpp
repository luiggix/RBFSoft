/*! 
 ***************************************************************************
 *  \file    heat2DTest.cpp
 *  Poisson equation in 2D.
 *  The equation is: 
 *  \f[ \frac{\partial T}{\partial^2 x}+\frac{\partial T}{\partial^2 x}=0 \f]
 *  and the boundary conditions are as shown in the next figure:
 * \image html  Xrectx1_y2.png "Geometry and boundary conditions" width=2.5cm 
 * \image latex Xrectx1_y2.eps "Geometry and boundary conditions" width=5cm 
 *  \par Compiling and running
 *  Modify the variables FLENS, ATLAS and RBF in the file \c rules.in, 
 *  according to your installation, then type the next commands:
 *  \verbatim
    % make
    % ./heat2D_1 \endverbatim   
 * The \c input file contains the initial data to setup the problem.
 * \par Output
 * \c xy_knots.dat (x,y) coordinates of random points;
 * \c solution.dat (x,y,u) evaluation of the solution in a grid;
 * \c errordis.dat (x,y,e) error distribution in a grid;
 * \c exactsol.dat (x,y,u) exact solution in a grid;
 ***************************************************************************
 *  \author  Luis M. de la Cruz [ Thu Sep  6 14:35:41 BST 2007 ]
 ***************************************************************************
 */
#include <iostream>
#include <fstream>
#include <string>

#include "Traits.hpp"
#include "Knots/RectangleKnots.hpp"
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
bool evalSolution(const Vec&, const Vec&, const Vec&, prec_t);
prec_t exactSolution(prec_t, prec_t);

//==========================================================================
//                            GLOBAL DATA
//
const int Ngx = 21;    
const int Ngy = 41;    
int Nx, Ny, N, NI, NB; 
prec_t hx, hy;         

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

    prec_t c, ep;
    int max_knots, rtype, fsol, layer;

// ----- Reading data from "input" file

    std::ifstream input_file ("input");
    input_file >> hx;   // length in x-axis
    input_file >> hy;   // length in y-axis
    input_file >> Nx;   // Boundary points in x-axis
    input_file >> Ny;   // Boundary points in y-axis
    input_file >> rtype;// 0 - Uniform, 1 - Random, 2 - Random "uniform"
    input_file >> ep;   // Randomness of the interior points (0-1.0)
    input_file >> layer;// Ghost (1) or not Ghost (0)
    input_file >> c;    // Shape parameter
    input_file >> fsol; // 1 Gauss, 2 GMRES, 3 GMRES-Jacobi, 4 GMRES-ACBF
    input_file >> max_knots; // Max. number of support knots
    input_file.close();
    random_t RT = static_cast<random_t>(rtype);

// ----- Point generation 

    std::cout << "\n +-----+  "
	      << "\n | Calculating knots ... ";    

    RectangleKnots<prec_t> rect(hx, Nx, hy, Ny, RT, layer);
    N  = rect.getTotalKnots();
    NI = rect.getInteriorKnots();
    NB = rect.getBoundaryKnots();
    rect.setRandomness(ep);
    rect.print();
    time.tic();
    rect.constructKnots();
    rect.writeToFile("xy_knots.dat");

    Vec x = rect.getKnots(X);
    Vec y = rect.getKnots(Y);

    std::cout << "\n | Knots generation elapsed time = " << time.toc()
	      << "\n +-----+  ";  
  
    if (c < 0) c = 1.0 / sqrt (1.0 * N);
    std::cout << "\n +-----+"
	      << "\n | c = " << c
	      << "\n +-----+";

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

// Edge 3: Dirichlet u = 100
    int bi, ei;
    int side_x = Nx - 2, side_y = Ny - 2;    
    bi = NI + 2 * side_y + 1;
    ei = bi + side_x;
    for(int i = bi; i < ei; ++i)
	u(i) = 100.0;

// Corners 1 and 2: Dirichlet u = 100
    u(N-3) = 100;
    u(N-1) = 100;

    std::cout << "\n | Elapsed time = " << time.toc()
	      << "\n +-----+  ";    

// ----- Solve the Linear System  G * lambda = u.

    std::cout << "\n +-----+  "
	      << "\n | Solving the linear system : ";

    prec_t tol = 1.0e-8;
    int niter;

    IdentityPrec<Mat> Iprec(G);
    DiagonalPrec<Mat> Jprec(G);
    ACBFPrec<Mat, RectangleKnots<prec_t> > Aprec(G);

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
	    Aprec.construct(rect, max_knots, 0, true);
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
    evalSolution(x, y, lambda, c);
    std::cout << "\n | Elapsed time = " << time.toc() ;   
    std::cout << "\n +-----+  ";

#ifdef WITH_GNUPLOT
    int pausa;
    GNUplot plotter;
    plotter("set contour");
    plotter("set cntrparam level discrete 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95");
    plotter("splot \"solution.dat\" w l");
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
	for(int i = NI+1; i <= N; ++i)
	  G(i,j) = rbf ( x(i), y(i), x(j), y(j) );   


}

//--------------------------------------------------------------------------
// Evaluate the numerical solution on a mesh of Ngx X Ngy
//
bool evalSolution(const Vec& x, const Vec& y, const Vec& lambda, prec_t c)
{
    std::ofstream sol_file ("solution.dat");
    std::ofstream exa_file ("exactsol.dat");
    std::ofstream err_file ("errordis.dat");
    prec_t xmi,ymi;
    prec_t dx = hx / (Ngx - 1);
    prec_t dy = hy / (Ngy - 1);
    prec_t diff, sol, exa, sum = 0.0, max = 0.0;

    for(int j = 0; j < Ngy; ++j) {
      for(int i = 0; i < Ngx; ++i) {
	xmi = dx * i;
	ymi = dy * j;
	
	sol = RBF::eval(xmi, ymi, x, y, lambda, RBF::MQ<prec_t, 2>(c));
	exa = exactSolution(xmi, ymi);
	
// ----- Dirichlet boundary conditions
	if ( i == 0 || i == Ngx - 1) { sol = 0.0; exa = 0.0; }
	if ( j == Ngy - 1) { sol = 0.0; exa = 0.0; }
	if ( j == 0) { sol = 100.0; exa = 100.0; }
	
	diff = fabs(sol - exa);
	
	if (diff > max) max = diff;

	sol_file << xmi << "\t" << ymi << "\t" << sol << std::endl;   
	exa_file << xmi << "\t" << ymi << "\t" << exa << std::endl;
	err_file << xmi << "\t" << ymi << "\t" << diff << std::endl;
	
	sum += (diff * diff);
	  
      }
      sol_file << std::endl;
      exa_file << std::endl;
      err_file << std::endl;
    }
    std::cout << "\n | RMS error = " << sqrt (sum /(Nx * Ny));
    std::cout << "\n | Max error = " << max;
    sol_file.close();
    exa_file.close();
    err_file.close();

    return 0;
}
 
//--------------------------------------------------------------------------
// Evaluate the exact solution in the (x, y) point.
//
prec_t exactSolution(prec_t x, prec_t y)
{
    prec_t c1, s1, s2, s3;
    prec_t sum = 0.0, hy = 2.0, hx = 1.0;
    int Nmax = 20;

    for(int n = 1; n <= Nmax; n += 2) {
	c1 = 2.0 / (n * PI);
	s1 = sinh( n * PI * (hy - y) / hx );
	s2 = sinh( n * PI * hy / hx );
	s3 = sin( n * PI * x / hx );
	sum += c1 * s1 * s3 / s2;
    }
    return 2.0 * 100.0 * sum ;
}

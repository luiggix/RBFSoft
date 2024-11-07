/*! 
 ***************************************************************************
 *  \file    heat2D_blocks.cpp
 *  Poisson equation in 2D (blocks).
 *  Same as heat2D.cpp but using RectBlocksKnots<> class.
 ***************************************************************************
 *  \author  Luis M. de la Cruz [ Wed Jan 30 16:21:16 GMT 2008 ]
 ***************************************************************************
 */

#include <ctime>    
#include <cstdlib>  
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

#include "Traits.hpp"
#include "Knots/RectBlocksKnots.hpp"
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
const int Ngx = 21;    // Size of the grid
const int Ngy = 41;    // to evaluate the solution
int Nx, Ny, N, NI, NB; // Number of points involved in the calculation
prec_t hx, hy;         // Size of the domain

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

    int max_knots, rtype, layer;
    int NBL, NBR, NBS, NBN, NIx, NIy;
    prec_t ep, c;

// ----- Reading data from "input" file

    std::ifstream input_file ("input_blocks");
    input_file >> hx;   // length in x-axis
    input_file >> hy;   // length in y-axis
    input_file >> Nx; // Boundary points in x-axis
    input_file >> Ny; // Boundary points in y-axis
    input_file >> NBL >> NBR >> NBS >> NBN >> NIx >> NIy;
    input_file >> rtype;
    input_file >> ep;   // Randomness of the interior points (0-1.0)
    input_file >> layer;// Ghost (1) or not Ghost (0)
    input_file >> c;    // Shape parameter
    input_file >> max_knots; // Max. number of support knots
    input_file.close();

    random_t RT = static_cast<random_t>(rtype);
    
// ----- Point generation 
    std::cout << "\n +-----+  "
	      << "\n | Calculating knots ... ";
    
    RectBlocksKnots<prec_t> rect(hx, Nx, hy, Ny, RT, layer,
				 NBL, NBR, NBS, NBN, 
				 NIx, NIy);
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

// ----- Solve the Linear System  G * lambda = u
    std::cout << "\n +-----+  "
	      << "\n | Solving the linear system : ";

    prec_t tol = 1.0e-8;
    int niter;

    time.tic();
    std::cout << " GMRES ....";
    ACBFPrec<Mat, RectBlocksKnots<prec_t> > precond(G);
    precond.construct(rect, max_knots, 8, true);    
    std::cout << "\n | Elapsed time - preconditioner : " 
	      << time.toc();
    time.tic();

    niter = gmres (G, lambda, u, precond, N-1 , tol);    
    std::cout << "\n | GMRES Iterations : " << niter;

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
//    plotter("set contour");
//    plotter("set cntrparam level discrete 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95");
//    plotter("splot \"solution.dat\" w l");
    plotter("p \"xy_knots.dat\" w p");
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
// Evaluate the numerical solution on a mesh of Ngx by Ngy
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


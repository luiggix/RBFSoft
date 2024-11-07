/*! 
 ***************************************************************************
 *
 *  \file convdiff2D.cpp
 *  Forced Convection in 2D.
 *  In this example the time-dependent convection-diffusion equation is 
 *  solved. The equation is written as follows:
 *  \f[ \frac{\partial T}{\partial t} + 
 *  u \frac{\partial T}{\partial x}  + v \frac{\partial T}{\partial y} =
 *  \Gamma\left( \frac{\partial^2 T}{\partial x^2} +
 *               \frac{\partial^2 T}{\partial y^2} \right)\f]
 *  where \f$(u, v)\f$ is a prescribed velocity field that fulfills the 
 *  continuity equation and is given by the next formula:
 *  \f{eqnarray*}{
 *    u(x,y) & = & -A\cos(\pi y) \sin(\pi \lambda x / l_x ) \mbox{and} \\
 *    v(x,y) & = & \frac{A \lambda}{l_x} \sin(\pi y) \cos(\pi \lambda x / l_x).
 *  \f}
 *  The initial condition is \f$T=0\f$ and the boundary conditions are as
 *  shown in the next figure:
 * \image html  solconv2D.png "Geometry, boundary conditions and solution" width=5cm 
 * \image latex solconv2D.eps "Geometry, boundary conditions and solution" width=5cm 
 * \par Input
 * The \c input file contains the initial data to setup the problem.
 * - \c hx, \c hy, \c Nx, \c Ny : Lenghts of the rectangle, and number of 
 *                                points in each axis.
 * - \c rtype, \c ep, \c layer : Type of knots distribution, randomness, layer
 *                                around the domain.
 * - \c max_iter time iterations, \c dt Step time step
 * - \c max_knots : neighbors used to construct the ACBF preconditioner.
 * - \c c : Shape parameter for MQ-RBF kernel (\f$c < 0\f$ implies \f$ c = 1/\sqrt{N}\f$).
 * - \c cboundary : Shape parameter for MQ-RBF kernel on the boundary (\f$c < 0\f$ implies \f$ c = 1/\sqrt{N}\f$).
 * - \c A max magnitude of the velocity.
 * - \c Nprofile number of points to make the velocity profiles.
 * \par Output
 * - \c xy_knots.dat coordinates of random points;
 * - \c solution.dat (x,y,u) evaluation of the solution in a grid.
 * - \c profBF_x.dat, \c profRBF_y.dat velocity profiles through the middle of the cavity.
 ***************************************************************************
 *  \author  Luis M. de la Cruz Wed [ Dec 19 14:49:56 GMT 2007 ]
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
void initialVelocity(Vec&, Vec& , const Vec&, const Vec&, prec_t, prec_t);

template<typename RBF, 
	 typename RBF_1DX, typename RBF_1DY, 
	 typename RBF_2DX, typename RBF_2DY>
void fillMatrices(Mat&, const Vec&, const Vec&,
		  RBF rbf, 
		  RBF_1DX rbf_1dx, RBF_1DY rbf_1dy,
		  RBF_2DX rbf_2dx, RBF_2DY rbf_2dy, 
		  const Vec&, const Vec&, prec_t, prec_t,
		  RBF, RBF_1DY);

prec_t evalSolution(Vec&, const Vec&, const Vec&, const Vec&, prec_t, prec_t);

//==========================================================================
//                            GLOBAL DATA
//
int Nx, Ny, N, NI, NB, NG;
prec_t hx, hy;

int stop;
#ifdef WITH_GNUPLOT
    GNUplot plotter;
#endif

//==========================================================================
//                            MAIN FUNCTION
//
int main( int argc, char * argv[])
{
    timer time;
    prec_t total_time = 0;

    std::cout << "\n\n"
	      << " +----------------------------------------+\n"
	      << " |       TUNA::RBF FOR PDE SOLVING        |\n"
	      << " +----------------------------------------+\n"
	      << "\n";

    int max_iter, max_knots, rtype, layer, Nprofile;
    prec_t ep, c, dt, A, cboundary;
    prec_t rolls;

    std::ifstream input_file ("input");
    input_file >> hx   // length in x-axis
	       >> hy   // length in y-axis
	       >> Nx   // Boundary points in x-axis
	       >> Ny   // Boundary points in y-axis
	       >> rtype
	       >> ep   // Randomness of the interior points (0-0.5)
	       >> layer 
	       >> max_iter // Number of time steps
	       >> dt   // Time step
	       >> max_knots 
	       >> c
	       >> cboundary
	       >> A   // amplitude
	       >> Nprofile
	       >> rolls;
    input_file.close();
    random_t RT = static_cast<random_t>(rtype);

    std::cout << "\n >--> hx = " << hx << "\t hy = " << hy
	      << "\t Nx = " << Nx << "\t Ny = " << Ny
	      << "\n >--> ep = " << ep << "\t layer = " << layer
	      << "\t max = " << max_iter << "\t dt = " << dt
	      << "\t maxk = " << max_knots << "\t A = " << A
	      << "\n";

 // ----- Point generation 

    std::cout << "\n +-----+  "
	      << "\n | Calculating knots ... ";  
   
    RectangleKnots<prec_t> rect(hx, Nx, hy, Ny, RT, layer);

    bool corners = true;
    rect.setCorners(corners);

    N  = rect.getTotalKnots();
    NI = rect.getInteriorKnots();
    NB = rect.getBoundaryKnots();
    NG = rect.getGhostKnots();
    rect.setRandomness(ep);
    rect.print();
    time.tic();
    rect.constructKnots();      
    rect.writeToFile("xy_knots.dat");

    Vec x = rect.getKnots(X);
    Vec y = rect.getKnots(Y);

    std::cout << "\n | Knots generation elapsed time = " << time.toc()
	      << "\n +-----+  ";

   if (c < 0 ) c = 1.0 / sqrt (1.0 * N);
   if (cboundary < 0 ) cboundary = 1.0 / sqrt (1.0 * N);

   std::cout << "\n +-----+  "
	     << "\n | c = " << c
	     << "\n | cboundary = " << cboundary
	     << "\n +-----+  ";

   /*
#ifdef WITH_GNUPLOT
    plotter("set size 0.73, 1.0");
    plotter("set grid");
    //    plotter("set view 0,0,1.5,1");
    plotter("set view 0,0");
    plotter("plot \"xy_knots.dat\" w p");
    plotter("unset surface");
    plotter("unset clabel");
    std::cout << "\n >---> Press <1> and then <enter> to continue ";
    std::cin >> stop;           
#endif
   */

// ----- Initial condition 
    Vec uvel(N), vvel(N); 
    initialVelocity(uvel, vvel, x, y, A, rolls);

// ----- Fill the matrices using MQ and its derivatives.  

    std::cout << "\n +-----+  ";
    std::cout << "\n | Filling the linear system ...  ";
    
    Mat G(N, N);
    Vec u(N), lambda(N);

    time.tic();
    fillMatrices(G, x, y,
		 RBF::MQ<prec_t, 2>(c), 
		 RBF::MQ_1DX<prec_t, 2>(c), 
		 RBF::MQ_1DY<prec_t, 2>(c), 
		 RBF::MQ_2DX<prec_t, 2>(c), 
		 RBF::MQ_2DY<prec_t, 2>(c), 
		 uvel, vvel, dt, 1.0,
		 RBF::MQ<prec_t, 2>(cboundary),
		 RBF::MQ_1DY<prec_t, 2>(cboundary));

    std::cout << "\n | Elapsed time = " << time.toc();

    time.tic();
// ----- Preconditioner construction
    std::cout << "\n | Construction of the preconditioner....";

    ACBFPrec<Mat, RectangleKnots<prec_t> > precond(G);
    precond.construct(rect, max_knots, 8, true);
//    IdentityPrec<Mat> precond(G); // This prec adds some spurious oscilla...
//    DiagonalPrec<Mat> precond(G);

    std::cout << "\n | Elapsed time = " << time.toc();
    std::cout << "\n +-----+  ";

// ----- Boundary conditions.
    int bi, ei;
    int side_x = Nx - 2, side_y = Ny - 2;

// Edge 1: Dirichlet Temp = 0.5;
    bi = NI + 1;  // MQ
    ei = bi + side_y;
    for(int i = bi; i < ei; ++i)
	u(i) = 0.5;
	
// Edge 2: Dirichlet Temp = -0.5;
    bi = ei;
    ei = bi + side_y;
    for(int i = bi; i < ei; ++i)
	u(i) = -0.5;

    bi = ei;
    ei = bi + side_x;
    for(int i = bi; i < ei; ++i)
	u(i) = 0.0;

    if (corners) {
// Corners:
	u(N-3) = 0.5;
	u(N-2) = 0.5;
	u(N-1) = -0.5;
	u(N) = -0.5;
    }
// ----- Main Loop

    int niter, NtotalIter;
    double dniter = 0;
    prec_t tol = 1.0e-8, tolerance = 1.0e-8, error = 0; 

#ifdef WITH_GNUPLOT
    plotter("set grid");
    plotter("set contour");
    plotter("set cntrparam level discrete -.4, -.3, -.2, -.1, 0, .1, .2, .3, .4");
    plotter("set view 0,0");
    //    plotter("set view 0,0,1.5,1");
    plotter("set size 0.73, 1.0");
    plotter("unset surface");
    plotter("unset clabel");
#endif

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
	total_time += time.toc();
	std::cout << "\n | " << t <<"\t GMRES Iterations : " << niter;
/* */
	dniter += niter;
	error = evalSolution(u, x, y, lambda, c, dt);
	std::cout << "\t error = " << error;

#ifdef WITH_GNUPLOT
//	plotter("plot \"error.dat\" w lp");
//	plotter("splot \"solgrid.dat\" notitle w l, \"temp.05001\" w l");
	plotter("splot \"solgrid.dat\" notitle w l");
//	plotter("splot \"solution.dat\" w l");
#endif
	if (error <= tolerance) break;

	NtotalIter = t;
//	std::cin >> stop;       
    }

    std::cin >> stop;           

    prec_t t_x, t_y, dl = 1.0 / (Nprofile-1);
    std::ofstream profile_x("profRBF_x.dat"), profile_y("profRBF_y.dat");
    for(int i = 0; i < Nprofile; ++i) {	
	t_x = RBF::eval(0.5, dl*i, x, y, lambda, RBF::MQ<prec_t,2>(c), 0);
	t_y = RBF::eval(dl*i, 0.5, x, y, lambda, RBF::MQ<prec_t,2>(c), 0);
	profile_y << t_x << "\t" << dl * i << "\n";
	profile_x << dl * i << "\t" << t_y << "\n";
    }
    profile_x.close();
    profile_y.close();

#ifdef WITH_GNUPLOT
    plotter("plot \"profRBF_x.dat\" w lp");
    std::cout << "\n | Press <1> and then <enter> "; std::cin >> stop; 
    plotter("plot \"profRBF_y.dat\" w lp");
    std::cout << "\n | Press <1> and then <enter> "; std::cin >> stop; 
#endif

    std::cout << "\n Total time ( Ax = b ) = " << total_time;
    std::cout << "\n Average GMRES iter = " << dniter / NtotalIter;
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
//    std::ofstream uvel_file ("uvel.dat");
//    std::ofstream vvel_file ("vvel.dat");
    int cnt = 1;
    for(int j = 2; j < Ny; ++j) {
	for(int i = 2; i < Nx; ++i, ++cnt) {
	    u(cnt) = -A * cos(PI * y(cnt)) * sin(PI * rolls * x(cnt) / hx); 
	    v(cnt) = A2 * sin(PI * y(cnt)) * cos(PI * rolls * x(cnt) / hx);
//	    uvel_file << x(cnt) << "\t" << y(cnt) << "\t" << u(cnt) << "\n";
//	    vvel_file << x(cnt) << "\t" << y(cnt) << "\t" << v(cnt) << "\n";
	}	    
//	uvel_file << "\n";
//	vvel_file << "\n";	
    }
//    uvel_file.close();
//    vvel_file.close();
}
//--------------------------------------------------------------------------
// Fill the matrix of the system
//
template<typename RBF, 
	 typename RBF_1DX, typename RBF_1DY, 
	 typename RBF_2DX, typename RBF_2DY>
void fillMatrices(Mat& G, const Vec& x, const Vec& y, 
		  RBF rbf, 
		  RBF_1DX rbf_1dx,
		  RBF_1DY rbf_1dy,
		  RBF_2DX rbf_2dx, 
		  RBF_2DY rbf_2dy, 
		  const Vec& uvel, const Vec& vvel, prec_t dt, prec_t D,
		  RBF rbf_b, RBF_1DY rbf_1dy_b)
{
// ----- Wl matrix

    for(int j = 1; j <= N; ++j) 
	for(int i = 1; i <= NI; ++i) {  
	    G(i,j) = rbf ( x(i), y(i), x(j), y(j) ) +
		dt * ( uvel(i) * rbf_1dx ( x(i), y(i), x(j), y(j) ) +
		       vvel(i) * rbf_1dy ( x(i), y(i), x(j), y(j) ) -
		       D * ( rbf_2dx ( x(i), y(i), x(j), y(j) ) + 
			     rbf_2dy ( x(i), y(i), x(j), y(j) ) ) );
	}

// ----- Wb matrix
    int side_x = Nx - 2, side_y = Ny - 2;

    int bi, ei;
    for(int j = 1; j <= N; ++j) {
// Edges 1 and 2: Dirichlet
	bi = NI + 1;
	ei = bi + 2 * side_y;
	for(int i = bi; i < ei; ++i) {
	    G(i,j) = rbf_b ( x(i), y(i), x(j), y(j) );   
//	    G(i,j) = rbf ( x(i), y(i), x(j), y(j) );   
	}

// Edges 3 and 4: Neumman
	bi = ei;
	ei = bi + 2 * side_x;
	for(int i = bi; i < ei; ++i) {
	    G(i,j) = rbf_1dy_b ( x(i), y(i), x(j), y(j) );
//	    G(i,j) = rbf_1dy ( x(i), y(i), x(j), y(j) );
	}
/* */
// Corners: Dirichlet
	G(ei    ,j) = rbf_b ( x(ei    ), y(ei    ), x(j), y(j) ); // Corner 1
	G(ei + 1,j) = rbf_b ( x(ei + 1), y(ei + 1), x(j), y(j) ); // Corner 2
	G(ei + 2,j) = rbf_b ( x(ei + 2), y(ei + 2), x(j), y(j) ); // Corner 3
	G(ei + 3,j) = rbf_b ( x(ei + 3), y(ei + 3), x(j), y(j) ); // Corner 4

/* */

    }
}
//--------------------------------------------------------------------------
// Evaluate the numerical solution on a mesh of Nx by Ny
//

prec_t evalSolution(Vec& u, const Vec& x, const Vec& y, 
		    const Vec& lambda, prec_t c, prec_t dt)
{
// Look, I just need to Update the values of interior nodes, because these
// are the "source" for the next iteration. The boundary values are prescribed
// and do not change during the calculation.
    prec_t uold, unew, error = 0, beta = 1.0;
    for(int i = 1; i <= NI; ++i) { 
	uold = u(i);
	u(i) = RBF::eval(x(i), y(i), x, y, lambda, RBF::MQ<prec_t,2>(c),0);
//	u(i) = beta * unew + (1 - beta) * uold;
	error += fabs( uold - u(i) );
    }

    int side_x = Nx - 2, side_y = Ny - 2;
    int bi, ei;

/* *
// Does it is really needed?????

// Edge 3 : Neumman dT/dy = 0
    bi = NI + 2 * side_y + 1;
    ei = bi + side_x;
    for(int i = bi, ii = 1; i < ei; ++i, ++ii)
//	u(i) = u(ii);
	u(i) = RBF::eval(x(i), y(i), x, y, lambda, RBF::MQ<prec_t,2>(c),0);


// Edge 4 : Neumman dT/dy = 0    
    bi = ei;
    ei = bi + side_x;
    int bii = NI - side_x + 1; 
    for(int i = bi, ii = bii; i < ei; ++i, ++ii) 
//	u(i) = u(ii);
	u(i) = RBF::eval(x(i), y(i), x, y, lambda, RBF::MQ<prec_t,2>(c),0);

/* */

// Solution in the interior collocation nodes.
    std::ofstream sol_file ("solution.dat");    
    for(int i = 1; i <= NI; ++i)
	sol_file << x(i) <<" "<< y(i) <<"  "<< u(i) << std::endl;   

/* *
// Prescribed values on the boundary.
    std::ofstream b_file ("solutionB.dat");    
    for(int i = NI+1; i <= N; ++i)
	b_file << x(i) <<" "<< y(i) <<"  "<< u(i) << std::endl;   
    b_file.close();
/* */

    sol_file.close();


// Evaluation of the solution in a grid of 20 by 20
    int Ngrid = 20;
    std::ofstream grid_file ("solgrid.dat");    
    prec_t ddx = 1. / Ngrid;
    prec_t ddy = 1. / Ngrid;
    for(int i = 0; i <= Ngrid; ++i) {
	for(int j = 0; j <= Ngrid; ++j) 
/* */
	    if (i == 0) {
		grid_file << i * ddx << "\t" << j * ddy << "\t"<< 0.5 << "\n";
	    } else if (i == Ngrid) {
		grid_file << i * ddx << "\t" << j * ddy << "\t"<< -0.5 << "\n";
/* *
	    } else if (j == 0) {
		grid_file << i * ddx << "\t" << j * ddy << "\t"
			  << RBF::eval( i * ddx, ddy, 
					x, y, lambda, 
					RBF::MQ<prec_t,2>(c)) << "\n";
	    } else if (j == Ngrid) {
		grid_file << i * ddx << "\t" << j * ddy << "\t"
			  << RBF::eval( i * ddx, (j-1) * ddy, 
					x, y, lambda, 
					RBF::MQ<prec_t,2>(c)) << "\n";
/* */
	    } else {
		
		grid_file << i * ddx << "\t" << j * ddy << "\t" 
			  << RBF::eval( i * ddx, j * ddy, 
					x, y, lambda, 
					RBF::MQ<prec_t,2>(c)) << "\n";
/* */
	    }
/* */
	grid_file << "\n";
    }
    grid_file.close();

    return error / N;
} 


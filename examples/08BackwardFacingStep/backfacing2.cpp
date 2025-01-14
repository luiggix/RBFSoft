/*! 
 ***************************************************************************
 *  \file backfacing2.cpp
 *  Backward-acing Step in 2D.
 *  In this example the Backward-facing step is solved in the vorticity -
 *  stream function formulation. The equations to solve are:
 *  \f[ \frac{\partial^2 \psi}{\partial x^2} + 
 *      \frac{\partial^2 \psi}{\partial y^2} = w \f]
 * and
 *  \f[ \frac{\partial^2 w}{\partial x^2} + 
 *      \frac{\partial^2 w}{\partial y^2} = 
 *      Re \left( u \frac{\partial w}{\partial x} + 
 *                v \frac{\partial w}{\partial y} \right) \f]
 *  where \f$ \psi \f$ is the stream function and the vorticity and the 
 *  velocity are defined as follows
 *  \f[ (u, v) = (\frac{\partial \psi}{\partial y}, 
                     -\frac{\partial \psi}{\partial x}) \\
     w = \frac{\partial v}{\partial x} - \frac{\partial u}{\partial y} \f]
 *  The boundary conditions are the typical for the backward-facing step, 
 *  see for example [].
 * \image html solbackfacing.png "Streamfunction and Vorticity" width=5cm 
 * \image latex solbackfacing.eps "Streamfunction and Vorticity" width=5cm 
 *  \par Input
 * The \c input file contains the initial data to setup the problem.
 * - \c hx, hy, Nx, Ny Size of the cavity and number of points on each axis.
 * - \c rtype, \c ep point distribution, randomness
 * - \c dt, \c max_iter, \c max_knots time step, time iterations, neighbors
 * - \c c shape parameter for MQ-RBF kernel (\f$c < 0\f$ implies \f$ c = 1/\sqrt{N}\f$).
 * - \c Re Reydolds number.
 * - \c beta_s, beta_v, beta_vb under-relaxation parameters for streamfunction, vorticity and boundary conditions for vorticity.
 * - \c gmres_tol, \c tolerance Tolerance for GMRES, global tolerance
 * - \c scale_x, scale factor for x-coordinates.
 * \par Output
 * - \c xy_knots.dat coordinates of random points;
 * - \c solStream.dat (x,y,u) Streamfunction.
 * - \c solVort.dat (x,y,u) Vorticity.
 ***************************************************************************
 *  \author  Luis M. de la Cruz Wed [ Mon Apr  7 09:58:44 BST 2008 ]
 ***************************************************************************
 */

#include <ctime>    
#include <cstdlib>  
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

#include "Traits.hpp"
#include "Knots/RectangleKnots.hpp"

#include "Solvers/Gauss.hpp" 
#include "Solvers/Krylov/diagonalprec.hpp"
//#include "Solvers/Krylov/identityprec.hpp"
#include "Solvers/Krylov/acbfprec.hpp"
#include "Solvers/Krylov/gmres.hpp"
#include "RBF/RBF.hpp"
#include "RBF/MQ.hpp"
//#include "RBF/TPS.hpp"
#include "GNUplot.hpp"

//==========================================================================
//                            FUNCTION PROTOTYPES
//
template<typename RBF, 
	 typename RBF_2DX, 
	 typename RBF_2DY>
void streamLaplace(Mat&, const Vec&, const Vec&, RBF, RBF_2DX, RBF_2DY, 
		   prec_t);

template<typename RBF, 
	 typename RBF_1DX, typename RBF_1DY,
	 typename RBF_2DX, typename RBF_2DY>
void vortEquation(Mat&, const Vec&, const Vec&, RBF,
		  RBF_1DX, RBF_1DY, RBF_2DX, RBF_2DY, 
		  const Vec& u, const Vec& v, prec_t, prec_t);

void evalVelocity(Vec&, Vec&, const Vec&, const Vec&, const Vec&, prec_t);

prec_t evalSolution(Mat&, const Vec&, const Vec&, const Vec&, prec_t, prec_t, string);

void vorticityBC(Vec&, Mat&, const Mat&, int*, prec_t, Vec&, prec_t);

void writeSolution(const Vec&, const Vec&, const Vec&, prec_t, string);

bool writeMatrix(const Mat&, string, int, int, int, int);

//==========================================================================
//                            GLOBAL DATA
//
int Nx, Ny, N, NI, NB, NG;
prec_t hx, hy, dx, dy;

prec_t top_wall_value = 0.5;
int NgridX = 300, NgridY = 30;

#ifdef WITH_GNUPLOT
GNUplot plotStream, plotVort, plotter;
int stop;
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

//    std::cout << "\n >---> top_wall_value = "; std::cin >> top_wall_value;

    int max_iter, max_knots, rtype;
    prec_t ep, c, Re, beta_s, beta_v, beta_vb, dt, scale_x;
    prec_t gmres_tol, tolerance, error_s = 1, error_v = 1; 

    std::ifstream input_file ("input2");
    input_file >> hx        // length in x-axis
	       >> hy        // length in y-axis
	       >> Nx        // Boundary points in x-axis
	       >> Ny        // Boundary points in y-axis
	       >> rtype     // Type of knots distribution
	       >> ep        // Randomness of the interior points (0-1)
	       >> dt        // Pseudo-time step
	       >> max_iter  // Number of time steps
	       >> max_knots // Neighbor knots for ACBF precond
	       >> c         // Shape parameter
	       >> Re        // Reynolds number
	       >> beta_s    // Under-relaxation coefficient for streamfunction
	       >> beta_v // Under-relaxation coefficient for vorticity
	       >> beta_vb // Under-relaxation coefficient for boundary vort
	       >> gmres_tol // Tolerance for GMRES
	       >> tolerance// Global tolerance
	       >> scale_x; // Scale for x;
    input_file.close();
    random_t RT = static_cast<random_t>(rtype);

    dx = hx / (Nx - 1);
    dy = hy / (Ny - 1);

 // ----- Point generation 

    std::cout << "\n +-----+  "
	      << "\n | Calculating knots ... ";  
   
    RectangleKnots<prec_t> rect(hx, Nx, hy, Ny, RT);
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

    Vec x = rect.getKnots(X);  //  Coordinates of the
    Vec y = rect.getKnots(Y);  //  "physical" domain
    Vec xs = x / scale_x;   //  Scaled coordinates
    Vec ys = y / hy;   //  for RBF calculations

    std::cout << "\n | Knots generation elapsed time = " << time.toc()
	      << "\n +-----+  ";

    prec_t beta_c;
    if (c < 0)
	c = 1.0 / sqrt (1.0 * N);
    else {
	beta_c = c;
	c = dx < dy ? beta_c * dx : beta_c * dy;
    }

    std::cout << "\n +-----+  "
	      << "\n | c = " << c
	      << "\n | dx = " << dx << "\t dy = " << dy
	      << "\n | hx = " << hx << "\t hy = " << hy
	      << "\n +-----+  ";

// ----- Declaration of matrices and vectors

    Mat Gstream(N, N), Gvort(N,N);
    Vec b(N), lstream(N), lvort(N), uvel(N), vvel(N);
    Mat sgrid(Nx, Ny), vgrid(Nx, Ny);      

// ----- Calculate the Gram Matrix for streamfunction

    std::cout << "\n +-----+  "
	      << "\n | Filling the linear system: streamfunction ...  ";
    
    time.tic();
// ----- MQ 
/* */
    streamLaplace( Gstream, xs, ys, 
		   RBF::MQ<prec_t,2>(c),
		   RBF::MQ_2DX<prec_t,2>(c), 
		   RBF::MQ_2DY<prec_t,2>(c), dt );
// ----- TPS
/* *
    streamLaplace( Gstream, xs, ys, 
		   RBF::TPS<prec_t,2>(),
		   RBF::TPS_2DX<prec_t,2>(), 
		   RBF::TPS_2DY<prec_t,2>(), dt );
/* */

    std::cout << "\n | Elapsed time = " << time.toc()
	      << "\n +-----+  ";

// ----- Preconditioner construction

    std::cout << "\n +-----+  "
	      << "\n | Preconditioner for streamfunction ...";
    time.tic();
//    ACBFPrec<Mat, RectangleKnots<prec_t> > preStream(Gstream);
//    preStream.construct(rect, max_knots, 8, true);
    std::cout << "\n | Elapsed time = " << time.toc();
    std::cout << "\n +-----+  ";

    int bi, ei, niter, t = 0, id;
    int side_x = Nx - 2, side_y = Ny - 2;
    int point[9];
//
//           6                               7
//           |                               |
//           v                               v   
//       *---*---*---*---*---*---*---*---*---*---*
//       |                                       |
// 1 --> *                                       * <-- 3
//       |                                       |
//       *                                       *
//       |                                       |
// 8 --> *                                       * 
//       |                                       |
//       *                                       *
//       |                                       |
// 0 --> *                                       * <-- 2
//       |                                       |
//       *---*---*---*---*---*---*---*---*---*---*
//           ^                               ^ 
//           |                               |
//           4                               5
//
    point[0] = NI + 1;       point[1] = point[0] + side_y - 1;
    point[2] = point[1] + 1; point[3] = point[2] + side_y - 1;
    point[4] = point[3] + 1; point[5] = point[4] + side_x - 1;
    point[6] = point[5] + 1; point[7] = point[6] + side_x - 1;
    point[8] = point[0] + side_y / 2;

// ----- Inlet
    prec_t y_coord;
    for(int i = point[8], ii = Ny / 2 + 1; i <= point[1]; ++i, ++ii) {
	y_coord = y(i) - 0.5;
	sgrid(1, ii) = 2 * y_coord * y_coord * (3  - 4 * y_coord);
    }
    
// ----- Outlet
    for(int i = point[2], ii=2; i <= point[3]; ++i, ++ii) {
	y_coord = y(i) - 0.5 ;
	sgrid(Nx, ii) = 0.75 * y_coord - y_coord * y_coord * y_coord + 0.25;  
    }

//    top_wall_value = Re / 2 ;
    for(int i = point[6]; i <= point[7]; ++i)
	sgrid(i,Ny) = top_wall_value; 
    sgrid(1 ,Ny) = top_wall_value;
    sgrid(Nx,Ny) = top_wall_value;

#ifdef WITH_GNUPLOT
plotStream("set ylabel \"Streamfunction\"");
plotStream("set contour");
plotStream("unset surface");
plotStream("set grid");
plotStream("set view 0,0");
plotStream("set cntrparam level discrete 0.5, 0.4, 0.3, 0.2, 0.1, 0.01, 0.001, 0.0001, 0.0, -0.001, -0.01, -0.1");
//plotStream("set cntrparam level 20");

plotVort("set ylabel \"Vorticity\"");
plotVort("set contour");
plotVort("unset surface");
plotVort("set grid");
plotVort("set view 0,0");
//plotVort("set cntrparam level discrete -5, -4, -3, -2, -1, -0.5, 0, 0.5, 1, 2, 3");
plotVort("set cntrparam level 20");
#endif

// ----- Main Loop : BEGIN

    while ( (error_v > tolerance) && (t < max_iter) ) {
	t++;

// ----- Calculate the Gram Matrix for the vorticity

// ----- MQ
/* */
	vortEquation(Gvort, xs, ys, 
		     RBF::MQ<prec_t,2>(c),
		     RBF::MQ_1DX<prec_t,2>(c),
		     RBF::MQ_1DY<prec_t,2>(c),
		     RBF::MQ_2DX<prec_t,2>(c), 
		     RBF::MQ_2DY<prec_t,2>(c),
		     uvel, vvel, Re, dt);
// ----- TPS
/* *
	vortEquation(Gvort, xs, ys, 
		     RBF::TPS<prec_t,2>(),
		     RBF::TPS_1DX<prec_t,2>(),
		     RBF::TPS_1DY<prec_t,2>(),
		     RBF::TPS_2DX<prec_t,2>(), 
		     RBF::TPS_2DY<prec_t,2>(),
		     uvel, vvel, Re, dt);
/* */

// ----- Boundary conditions for the Vorticity

	vorticityBC(b, vgrid, sgrid, point, beta_vb, y, Re);

// ----- Solving for the vorticity
	id = 1;
	for(int j = 2; j < Ny; ++j)
	    for(int i = 2; i < Nx; ++i, ++id) 
		b(id) = vgrid(i,j);	
	
//	lvort = Solver::Gauss<prec_t>( Gvort, b );
	niter = gmres (Gvort, lvort, b, N-1 , gmres_tol);
	std::cout << "\n | " << t <<"\t Vorticity : GMRES Iterations : " 
		  << niter;

	error_v = evalSolution(vgrid, xs, ys, lvort, c, beta_v, "solVort.dat");
	std::cout << "\n | Vorticity : error = " << error_v;

#ifdef WITH_GNUPLOT
	plotVort("splot [][][] \"solVort.dat\" w l");
//	std::cout << "\n | Press <1> and then <enter> "; std::cin >> stop;  
#endif

	b = 0;
	id = 1;
	for(int j = 2; j < Ny; ++j)
	    for(int i = 2; i < Nx; ++i, ++id) 
		b(id) = sgrid(i,j) + dt * vgrid(i,j);

	for(int i = point[8], ii = Ny / 2 + 1; i <= point[1]; ++i, ++ii) 
	    b(i) = sgrid(1, ii);
    
	for(int i = point[2], ii = 2; i <= point[3]; ++i, ++ii) 
	    b(i) = sgrid(Nx, ii);
	
	for(int i = point[6], ii = 2; i <= point[7]; ++i, ++ii)
	    b(i) = sgrid(ii,Ny) = top_wall_value; 
	b(N-2) = top_wall_value;
	b(N) = top_wall_value;

//	lstream = Solver::Gauss<prec_t>( Gstream, b );
//	niter = gmres (Gstream, lstream, b, preStream, N-1 , gmres_tol);
	niter = gmres (Gstream, lstream, b, N-1 , gmres_tol);
	std::cout << "\n | " << t <<"\t Streamfuntion : GMRES Iterations : " 
		  << niter;

	error_s = evalSolution(sgrid, xs, ys, lstream, c, beta_s, "solStream.dat");
	std::cout << "\n | Streamfunction : error = " << error_s; 

	evalVelocity(uvel, vvel, xs, ys, lstream, c);

#ifdef WITH_GNUPLOT
//	plotVort("p [][] \"u_profiles.dat\" w p");
	plotStream("splot [][][] \"solStream.dat\" w l");
//	std::cout << "\n | Press <1> and then <enter> "; std::cin >> stop; 
#endif

    }    
#ifdef WITH_GNUPLOT
    std::cout << "\n | Press <1> and then <enter> "; std::cin >> stop; 
#endif

// Main Loop : END

// ----- Write profiles of the velocity
    int Ngrid = 50;
    prec_t ul, dl = 1.0 / Ngrid;
    std::ofstream ufile14 ("u_line_14.dat"), ufile30 ("u_line_30.dat");
//    ufile << 0.0 << "\t" << 0.0 << "\n";
    for(int i = 0; i <= Ngrid; ++i) {	
	ul = RBF::eval(7.0, dl*i, xs, ys, lstream,RBF::MQ_1DY<prec_t,2>(c),0);
	ufile14 << ul << "\t" << dl * i << "\n";
	ul = RBF::eval(15.0, dl*i, xs, ys, lstream,RBF::MQ_1DY<prec_t,2>(c),0);
	ufile30 << ul << "\t" << dl * i << "\n";
    }
//    ufile << 1.0 << "\t" << 1.0 << "\n";

    std::cout << "\n Total time ( Ax = b ) = " << total_time;
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
void streamLaplace(Mat& G, const Vec& x, const Vec& y, 
		   RBF rbf, RBF_2DX rbf_2dx, RBF_2DY rbf_2dy, prec_t dt)
{
// ----- Wl matrix
    for(int j = 1; j <= N; ++j) 
	for(int i = 1; i <= NI; ++i) {
	    G(i,j) =  rbf( x(i), y(i), x(j), y(j) ) -
		dt * ( rbf_2dx ( x(i), y(i), x(j), y(j) ) +
		       rbf_2dy ( x(i), y(i), x(j), y(j) ) );
	}
    
// ----- Wb matrix
    for(int j = 1; j <= N; ++j) 
	for(int i = NI+1; i <= N; ++i)
	    G(i,j) = rbf ( x(i), y(i), x(j), y(j) );   
}
//--------------------------------------------------------------------------
// Fill the matrix of the system
//
template<typename RBF, 
	 typename RBF_1DX, typename RBF_1DY,
	 typename RBF_2DX, typename RBF_2DY>
void vortEquation(Mat& G, const Vec& x, const Vec& y,
		  RBF rbf, 
		  RBF_1DX rbf_1dx, RBF_1DY rbf_1dy, 
		  RBF_2DX rbf_2dx, RBF_2DY rbf_2dy, 
		  const Vec& u, const Vec& v, prec_t Re, prec_t dt)
{
// ----- WL matrix
    std::cout << "\n >---> Re = " << Re
	      << "\t dt = " << dt;
    for(int j = 1; j <= N; ++j) 
	for(int i = 1; i <= NI; ++i) {
	    G(i,j) = rbf( x(i), y(i), x(j), y(j) ) -
		dt * ( rbf_2dx ( x(i), y(i), x(j), y(j) ) + 
		       rbf_2dy ( x(i), y(i), x(j), y(j) ) -  
		       Re * ( u(i) * rbf_1dx( x(i), y(i), x(j), y(j) ) +
			      v(i) * rbf_1dy( x(i), y(i), x(j), y(j) ) ) );
	}

// ----- WB matrix
    for(int j = 1; j <= N; ++j) 
	for(int i = NI+1; i <= N; ++i)
	    G(i,j) = rbf ( x(i), y(i), x(j), y(j) );   
}


//--------------------------------------------------------------------------
// Evaluate the velocity
//
void evalVelocity(Vec& u, Vec& v, const Vec& x, const Vec& y,
		  const Vec& lstream, prec_t c)
{
    for(int i = 1; i <= NI; ++i) {
/* */
	u(i) = RBF::eval(x(i), y(i), x, y, lstream,RBF::MQ_1DY<prec_t,2>(c),0);
	v(i) =-RBF::eval(x(i), y(i), x, y, lstream,RBF::MQ_1DX<prec_t,2>(c),0);
/* *
	u(i) = RBF::eval(x(i), y(i), x, y, lstream,RBF::TPS_1DY<prec_t,2>(),0);
	v(i) =-RBF::eval(x(i), y(i), x, y, lstream,RBF::TPS_1DX<prec_t,2>(),0);
/* */

    }
}

//--------------------------------------------------------------------------
// Evaluate the numerical solution on a mesh of Nx by Ny
//

prec_t evalSolution(Mat& fgrid, const Vec& x, const Vec& y, const Vec& lambda,
		    prec_t c, prec_t beta, string name)
{

    prec_t fnew, error = 0;
    int id = 1;
    for(int j = 2; j < Ny; ++j)
	for(int i = 2; i < Nx; ++i, ++id) {
	    fnew = RBF::eval(x(id), y(id), x, y, lambda, RBF::MQ<prec_t,2>(c));
//	    fnew = RBF::eval(x(id), y(id), x, y, lambda, RBF::TPS<prec_t,2>());
	    fnew = beta * fnew + (1 - beta) * fgrid(i,j);
	    error += fabs( fnew - fgrid(i, j) );
	    fgrid(i, j) = fnew;
	}

/* */
    std::ofstream grid_file (name.c_str());    
    prec_t xd, yd;
    prec_t ddx = hx / (Nx - 1);
    prec_t ddy = hy / (Ny - 1);
    id = 1;
    for(int j = 1; j <= Ny; ++j) { 
	for(int i = 1; i <= Nx; ++i, ++id) {
	    xd = (i-1) * ddx;
	    yd = (j-1) * ddy;
	    grid_file << xd <<"\t"<< yd << "\t" << fgrid(i,j) <<"\n";
	}
	grid_file << "\n";	    
    }
    grid_file.close();

/* *
    std::ofstream b_file ("solRand.dat");    
    prec_t fnew, error = 0;
    for(int i = 1; i <= NI; ++i) {	  
	fnew = RBF::eval(x(i), y(i), x, y, lambda, RBF::MQ<prec_t,2>(c),0);
	error += fabs( fnew - f(i) );
	f(i) = fnew;
	b_file << x(i) << "\t" << y(i) << "\t" << f(i) << "\n";
    }
    b_file.close();
/* */

    return error / N;
}

void writeSolution(const Vec& x, const Vec& y, const Vec& lambda,
		    prec_t c, string name)
{
    prec_t xd, yd;
    prec_t ddx = hx / (NgridX - 1.0);
    prec_t ddy = hy / (NgridY - 1.0);
    std::ofstream grid_file (name.c_str());    
    for(int j = 1; j <= NgridY; ++j) {
	for(int i = 1; i <= NgridX; ++i) {
	    xd = (i-1) * ddx;
	    yd = (j-1) * ddy;
	    grid_file << xd <<"\t"<< yd << "\t" 
		      << RBF::eval( xd / hx, yd/ hy, x, y, lambda, RBF::MQ<prec_t,2>(c))
//		      << RBF::eval( xd / hx, yd/ hy, x, y, lambda, RBF::TPS<prec_t,2>())
		      <<"\n";
	}
	grid_file << "\n";	    
    }
    grid_file.close();
} 


void vorticityBC(Vec& b, Mat& vgrid, const Mat& sgrid, int* point, 
		 prec_t beta, Vec& y, prec_t Re)
{
    b = 0;  // Source term

// Edges 1: Backward facing step - Inlet

/* */
// Inlet
    prec_t y_coord;
    for(int i = point[8], k = 2; i <= point[1]; ++i, ++k) {
	y_coord = y(i) - 0.5;
	b(i) = -12 * ( 1 - 4 * y_coord );

	b(i) += ( sgrid(1,k) - 2 * sgrid(2,k) + sgrid(3,k) ) / (dx * dx);
    }
/* */

    for(int i = point[0], k = 2; i < point[8]; ++i, ++k) {
//    for(int i = point[0], k = 2; i <= point[1]; ++i, ++k) {
// Second order approx. (Jensen formula):
	b(i) = -(4 * sgrid(2,k) - 0.5 * sgrid(3,k) - 3.5 * sgrid(1,k)) 
	    / (dx * dx);
    }
    
// Edge 2: Outlet
    for(int i = point[2], k = 2; i <= point[3]; ++i, ++k) {
/* */
	y_coord = y(i) - 0.5;
	b(i) = 6 * y_coord;
/* *
// Second order approx. (Jensen formula):
	b(i) = -(4 * sgrid(Nx-1,k) - 0.5 * sgrid(Nx-2,k) - 3.5 * sgrid(Nx,k)) 
	    / (dx * dx);
/* */
    }	

// Edges 3: Dirichlet
    for(int i = point[4], k = 2; i <= point[5]; ++i, ++k) {
// Second order approx. (Jensen formula):
	b(i) = -(4 * sgrid(k, 2) - 0.5 * sgrid(k, 3) - 3.5 * sgrid(k,1) ) 
	    / (dy * dy);
    }
    
// Edges 4: Dirichlet
    for(int i = point[6], k = 2; i <= point[7]; ++i, ++k) {
// Second order approx. (Jensen formula):
	b(i) = -(4 * sgrid(k, Ny-1) - 0.5 * sgrid(k, Ny-2) - 3.5 * sgrid(k,Ny))
	    / (dy*dy);
    }
// Corners:
// Average between wall neighbors:
/* */
   b(N-3) = 0.5 * ( b(point[0]) + b(point[4]) );
/* *
   b(N-2) = 0.5 * ( b(point[1]) + b(point[6]) );
   b(N-1) = 0.5 * ( b(point[2]) + b(point[5]) );
   b(N  ) = 0.5 * ( b(point[3]) + b(point[7]) );
/* */
   b(N-1) = -3;
   b(N  ) =  3;
   b(N-2) = 12;
/* */

// ----- Under-relaxation on the vorticity boundary conditions
    for(int i = point[0], k = 2; i <= point[1]; ++i, ++k) {
	b(i) = beta * b(i) + (1 - beta) * vgrid(1, k);
	vgrid(1, k) = b(i);
    }
    for(int i = point[2], k = 2; i <= point[3]; ++i, ++k) {
	b(i) = beta * b(i) + (1 - beta) * vgrid(Nx, k);
	vgrid(Nx, k) = b(i);
    }
    for(int i = point[4], k = 2; i <= point[5]; ++i, ++k) {
	b(i) = beta * b(i) + (1 - beta) * vgrid(k, 1);
	vgrid(k, 1) = b(i);
    }
    for(int i = point[6], k = 2; i <= point[7]; ++i, ++k) {
	b(i) = beta * b(i) + (1 - beta) * vgrid(k, Ny);
	vgrid(k, Ny) = b(i);	
    }
// Corners: under-relaxation
    b(N-3) = beta * b(N-3) + (1 - beta) * vgrid(1, 1); 
/* *
    b(N-2) = beta * b(N-2) + (1 - beta) * vgrid(1,Ny); 
    b(N-1) = beta * b(N-1) + (1 - beta) * vgrid(Nx,1); 
    b(N  ) = beta * b(N  ) + (1 - beta) * vgrid(Nx,Ny);
/* */
    vgrid(1 , 1 ) = b(N-3);
    vgrid(1 , Ny) = b(N-2);
    vgrid(Nx, 1 ) = b(N-1);
    vgrid(Nx, Ny) = b(N  );
}

bool writeMatrix(const Mat& G, string name, int bi, int ei, int bj, int ej)
{    
    std::ofstream Gram_file (name.c_str());
    for (int j = bj; j <= ej; ++j){
	for (int i = bi; i <= ei; ++i) 
	    Gram_file << i << " " << j << " " << G(i,j) << std::endl;
	Gram_file << std::endl;
    }

    return 0;
}


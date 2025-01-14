/*! 
 ***************************************************************************
 *  \file backfacing.cpp
 *  Backward-Facing Step in 2D (static).
 *  In this example the Backward-Facing Step problem is solved in the 
 *  vorticity - stream function formulation. The equations to solve are:
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
 *  The boundary conditions are:
 *  COMING SOON....(add a figure).
 *  \par Compiling and running
 *  Modify the variables FLENS, ATLAS and RBF in the file \c rules.in, 
 *  according to your installation, then type the next commands:
 *  \verbatim
    % make
    % ./lid_driven \endverbatim   
 * The \c input file contains the initial data to setup the problem.
 * \par Output
 * \c xy_knots.dat coordinates of random points;
 * \c solution.dat (x,y,u) evaluation of the solution in a grid.
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

//#include "Solvers/Gauss.hpp" 
#include "Solvers/Krylov/diagonalprec.hpp"
//#include "Solvers/Krylov/identityprec.hpp"
#include "Solvers/Krylov/acbfprec.hpp"
#include "Solvers/Krylov/gmres.hpp"
#include "RBF/RBF.hpp"
#include "RBF/MQ.hpp"
#include "GNUplot.hpp"

//==========================================================================
//                            FUNCTION PROTOTYPES
//
template<typename RBF, 
	 typename RBF_2DX, 
	 typename RBF_2DY>
void streamLaplace(Mat&, const Vec&, const Vec&, RBF, RBF_2DX, RBF_2DY);

template<typename RBF, 
	 typename RBF_1DX, typename RBF_1DY,
	 typename RBF_2DX, typename RBF_2DY>
void vortEquation(Mat&, const Vec&, const Vec&, RBF,
		  RBF_1DX, RBF_1DY, RBF_2DX, RBF_2DY, 
		  const Vec& u, const Vec& v, prec_t);

void evalVelocity(Vec&, Vec&, const Vec&, const Vec&, const Vec&, prec_t);

prec_t evalSolution(Mat&, const Vec&, const Vec&, const Vec&, prec_t, prec_t, string);

void vorticityBC(Vec&, Mat&, const Mat&, int*, prec_t, Vec&);

void writeSolution(const Vec&, const Vec&, const Vec&, prec_t, string);

bool writeMatrix(const Mat&, string);

//==========================================================================
//                            GLOBAL DATA
//
int Nx, Ny, N, NI, NB, NG;
prec_t hx, hy, dx, dy;

prec_t top_wall_value = 0.0;
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

    std::cout << "\n >---> top_wall_value = "; std::cin >> top_wall_value;

    int max_iter, max_knots, rtype;
    prec_t ep, c, Re, beta_s, beta_v, beta_vb;
    prec_t gmres_tol, tolerance, error_s = 1, error_v = 1; 

    std::ifstream input_file ("input2");
    input_file >> hx        // length in x-axis
	       >> hy        // length in y-axis
	       >> Nx        // Boundary points in x-axis
	       >> Ny        // Boundary points in y-axis
	       >> rtype     // Type of knots distribution
	       >> ep        // Randomness of the interior points (0-1)
	       >> max_iter  // Number of time steps
	       >> max_knots // Neighbor knots for ACBF precond
	       >> c         // Shape parameter
	       >> Re        // Reynolds number
	       >> beta_s    // Under-relaxation coefficient for streamfunction
	       >> beta_v // Under-relaxation coefficient for vorticity
	       >> beta_vb // Under-relaxation coefficient for boundary vort
	       >> gmres_tol // Tolerance for GMRES
	       >> tolerance;// Global tolerance
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

    Vec x = rect.getKnots(X);
    Vec y = rect.getKnots(Y);
    Vec xs = x;
    Vec ys = y;

    prec_t xscale = 1. / hx;
    prec_t yscale = 1. / hy;
    xs = x * xscale;
    ys = y * yscale;

    std::cout << "\n | Knots generation elapsed time = " << time.toc()
	      << "\n +-----+  ";

    if (c < 0 ) c = 1.0 / sqrt (1.0 * N);

    std::cout << "\n +-----+  "
	      << "\n | c = " << c
	      << "\n +-----+  ";

/* *
#ifdef WITH_GNUPLOT
    plotter("set grid");
    plotter("p \"xy_knots.dat\" w p");
    std::cout << "\n >---> Press <1> and then <enter> to continue \n";
    std::cin >> stop;
#endif
/* */

// ----- Declaration of matrices and vectors

    Mat Gstream(N, N), Gvort(N,N);
    Vec b(N), lstream(N), lvort(N), uvel(N), vvel(N);
    Mat sgrid(Nx, Ny), vgrid(Nx, Ny);      

// ----- Calculate the Gram Matrix for streamfunction

    std::cout << "\n +-----+  "
	      << "\n | Filling the linear system: streamfunction ...  ";
    
    time.tic();
    streamLaplace( Gstream, xs, ys, 
		   RBF::MQ<prec_t,2>(c),
		   RBF::MQ_2DX<prec_t,2>(c), 
		   RBF::MQ_2DY<prec_t,2>(c) );

    std::cout << "\n | Elapsed time = " << time.toc()
	      << "\n +-----+  ";

/* *
#ifdef WITH_GNUPLOT
    writeMatrix(Gstream, "Gstream.dat");
    plotStream("set grid");
    plotStream("splot \"Gstream.dat\" w l");
    std::cout << "\n >---> Pres <1> and then <enter> to continue \n";
    std::cin >> stop;
#endif
/* */
    

// ----- Preconditioner construction

    std::cout << "\n +-----+  "
	      << "\n | Preconditioner for streamfunction ...";
    time.tic();
    ACBFPrec<Mat, RectangleKnots<prec_t> > preStream(Gstream);
    preStream.construct(rect, max_knots, 8, true);
//    DiagonalPrec<Mat> preStream(Gstream);
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

    std::cout << " \n middle = " << point[8] 
	      << " \t y = " << y(point[8])
	      << "\n";
//    int kkk; std::cin >> kkk;

    prec_t y_coord;
    for(int i = point[8], ii = Ny / 2 + 1; i <= point[1]; ++i, ++ii) {
	y_coord = y(i) - 0.5;
	sgrid(1, ii) = 2 * y_coord * y_coord * (3  - 4 * y_coord);
    }
    
    for(int i = point[2], ii=2; i <= point[3]; ++i, ++ii) {
	y_coord = y(i) - 0.5 ;
	sgrid(Nx, ii) = 0.75 * y_coord - y_coord * y_coord * y_coord + 0.25;  
    }

    for(int i = point[6]; i <= point[7]; ++i)
	sgrid(i,Ny) = top_wall_value; 
    sgrid(1,Ny) = top_wall_value;
    sgrid(Nx,Ny) = top_wall_value;

/* *
    std::ofstream sgrid_file ("solStream.dat");    
    prec_t xd, yd;
    prec_t ddx = hx / (Nx - 1);
    prec_t ddy = hy / (Ny - 1);   
    for(int j = 1; j <= Ny; ++j) { 
	for(int i = 1; i <= Nx; ++i) {
	    xd = (i-1) * ddx;
	    yd = (j-1) * ddy;
	    sgrid_file << xd <<"\t"<< yd << "\t" << sgrid(i,j) <<"\n";
	}
	sgrid_file << "\n";	    
    }    
    sgrid_file.close();    

#ifdef WITH_GNUPLOT
	plotStream("splot [][][] \"solStream.dat\" w l");
	std::cout << "\n | Press <1> and then <enter> "; std::cin >> stop; 
#endif
/* */

#ifdef WITH_GNUPLOT
plotStream("set ylabel \"Streamfunction\"");
plotStream("set contour");
plotStream("unset surface");
plotStream("set grid");
plotStream("set view 0,0");
//plotStream("set cntrparam level discrete -0.11, -0.1, -0.09, -0.07, -0.05, -0.03, -0.01, -0.001, -0.0001, -0.00001, 0.0, 0.00001, 0.0001, 0.001, 0.01");
plotStream("set cntrparam level 50");

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

	vortEquation(Gvort, xs, ys, 
		     RBF::MQ<prec_t,2>(c),
		     RBF::MQ_1DX<prec_t,2>(c),
		     RBF::MQ_1DY<prec_t,2>(c),
		     RBF::MQ_2DX<prec_t,2>(c), 
		     RBF::MQ_2DY<prec_t,2>(c),
		     uvel, vvel, Re);

//	ACBFPrec<Mat, RectangleKnots<prec_t> > preVort(Gvort);
//	preVort.construct(rect, max_knots, 8);
	DiagonalPrec<Mat> preVort(Gvort);

/* *
#ifdef WITH_GNUPLOT
	writeMatrix(Gvort, "Gvort.dat");
	plotVort("set grid");
	plotVort("splot \"Gvort.dat\" w l");
	std::cout << "\n >---> Pres <1> and then <enter> to continue \n";
	std::cin >> stop;
#endif
/* */

// ----- Boundary conditions for the Vorticity

	vorticityBC(b, vgrid, sgrid, point, beta_vb, y);

/* */
#ifdef WITH_GNUPLOT
	std::ofstream b_file ("b0.dat");    
	for(int i = NI+1; i <= N; ++i)
	    b_file << x(i) << "\t" << y(i) << "\t" << b(i) << "\n";
	b_file.close();
/*
	plotter("set grid");
	plotter("set view 80,50");
	plotter("splot \"b0.dat\" w p");
	std::cout << "\n >---> Pres <1> and then <enter> to continue \n";
	std::cin >> stop;
*/
#endif
/* */

// ----- Solving for the vorticity
	niter = gmres (Gvort, lvort, b, preVort, N-1 , gmres_tol);
	std::cout << "\n | " << t <<"\t Vorticity : GMRES Iterations : " 
		  << niter;

	error_v = evalSolution(vgrid, xs, ys, lvort, c, beta_v, "solVort.dat");
	std::cout << "\n | Vorticity : error = " << error_v;

#ifdef WITH_GNUPLOT
//	plotVort("splot [][][] \"solVort.dat\" w l");
//	std::cout << "\n | Press <1> and then <enter> "; std::cin >> stop;  
#endif

	b = 0;
	id = 1;
	for(int j = 2; j < Ny; ++j)
	    for(int i = 2; i < Nx; ++i, ++id) 
		b(id) = -vgrid(i,j);

	for(int i = point[8], ii = Ny / 2 + 1; i <= point[1]; ++i, ++ii) 
	    b(i) = sgrid(1, ii);
    
	for(int i = point[2], ii = 2; i <= point[3]; ++i, ++ii) 
	    b(i) = sgrid(Nx, ii);
	
	for(int i = point[6], ii = 2; i <= point[7]; ++i, ++ii)
	    b(i) = sgrid(ii,Ny) = top_wall_value; 
	b(N-2) = top_wall_value;
	b(N) = top_wall_value;


	niter = gmres (Gstream, lstream, b, preStream, N-1 , gmres_tol);
	std::cout << "\n | " << t <<"\t Streamfuntion : GMRES Iterations : " 
		  << niter;
	
	error_s = evalSolution(sgrid, xs, ys, lstream, c, beta_s, "solStream.dat");
	std::cout << "\n | Streamfunction : error = " << error_s; 

	evalVelocity(uvel, vvel, xs, ys, lstream, c);

	std::ofstream ufile ("u_profiles.dat");
	for(int i = 1; i <= NI; i += 2)
	    ufile << x(i) + uvel(i) << "\t" << y(i) << "\n";
	ufile.close();

	writeSolution(xs, ys, lstream, c, "solStream.dat");
	    
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
// Fill the matrix of the system
//
template<typename RBF, 
	 typename RBF_1DX, typename RBF_1DY,
	 typename RBF_2DX, typename RBF_2DY>
void vortEquation(Mat& G, const Vec& x, const Vec& y,
		  RBF rbf, 
		  RBF_1DX rbf_1dx, RBF_1DY rbf_1dy, 
		  RBF_2DX rbf_2dx, RBF_2DY rbf_2dy, 
		  const Vec& u, const Vec& v, prec_t Re)
{
// ----- WL matrix
    for(int j = 1; j <= N; ++j) 
	for(int i = 1; i <= NI; ++i) {
	    G(i,j) = rbf_2dx ( x(i), y(i), x(j), y(j) ) + 
		rbf_2dy ( x(i), y(i), x(j), y(j) ) -
		Re * ( u(i) * rbf_1dx( x(i), y(i), x(j), y(j) ) +
		       v(i) * rbf_1dy( x(i), y(i), x(j), y(j) ) );
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
	u(i) = RBF::eval(x(i), y(i), x, y, lstream,RBF::MQ_1DY<prec_t,2>(c),0);
	v(i) =-RBF::eval(x(i), y(i), x, y, lstream,RBF::MQ_1DX<prec_t,2>(c),0);
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
	    fnew = beta * fnew + (1 - beta) * fgrid(i,j);
	    error += fabs( fnew - fgrid(i, j) );
	    fgrid(i, j) = fnew;
	}

/*
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
*/

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
		      <<"\n";
	}
	grid_file << "\n";	    
    }
    grid_file.close();
} 


void vorticityBC(Vec& b, Mat& vgrid, const Mat& sgrid, int* point, prec_t beta, Vec& y)
{
    b = 0;  // Source term

// Edges 1: Backward facing step - Inlet

/**
    prec_t y_coord;
    for(int i = point[8], k = 2; i <= point[1]; ++i, ++k) {
	y_coord = y(i) - 0.5;
	b(i) = -12 * ( 1 - 4 * y_coord );
    }
**/

    for(int i = point[0], k = 2; i <= point[1]; ++i, ++k) {
// Second order approx. (Jensen formula):
	b(i) = -(4 * sgrid(2,k) - 0.5 * sgrid(3,k) - 3.5 * sgrid(1,k)) 
	    / (dx * dx);
    }
    
// Edge 2: Outlet
    for(int i = point[2], k = 2; i <= point[3]; ++i, ++k) {
/**
	y_coord = y(i) - 0.5;
	b(i) = 6 * y_coord;
**/
// Second order approx. (Jensen formula):
	b(i) = -(4 * sgrid(Nx-1,k) - 0.5 * sgrid(Nx-2,k) - 3.5 * sgrid(Nx,k)) 
	    / (dx * dx);

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
/* */
   b(N-2) = 0.5 * ( b(point[1]) + b(point[6]) );
   b(N-1) = 0.5 * ( b(point[2]) + b(point[5]) );
   b(N  ) = 0.5 * ( b(point[3]) + b(point[7]) );
/* *
   b(N-2) = 12;
   b(N-1) = -3;
   b(N  ) =  3;
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
/* */
    b(N-2) = beta * b(N-2) + (1 - beta) * vgrid(1,Ny); 
    b(N-1) = beta * b(N-1) + (1 - beta) * vgrid(Nx,1); 
    b(N  ) = beta * b(N  ) + (1 - beta) * vgrid(Nx,Ny);
/* */
    vgrid(1 , 1 ) = b(N-3);
    vgrid(1 , Ny) = b(N-2);
    vgrid(Nx, 1 ) = b(N-1);
    vgrid(Nx, Ny) = b(N  );
}

bool writeMatrix(const Mat& G, string name)
{    
    std::ofstream Gram_file (name.c_str());
    for (int j = 1; j <= N; ++j){
	for (int i = 1; i <= N; ++i) 
	    Gram_file << i << " " << j << " " << G(i,j) << std::endl;
	Gram_file << std::endl;
    }

    return 0;
}


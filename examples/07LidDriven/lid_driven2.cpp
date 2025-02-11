/*! 
 ***************************************************************************
 *  \file lid_driven2.cpp
 *  Lid Driven Cavity in 2D.
 *  In this example the lid-driven cavity is solved in the vorticity -
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
 *  The boundary conditions are the typical for the lid-driven cavity, 
 *  see for example [].
 * \image html solliddriven.png "Streamfunction and Vorticity" width=5cm 
 * \image latex solliddriven.eps "Streamfunction and Vorticity" width=5cm 
 *  \par Input
 * The \c input file contains the initial data to setup the problem:
 * - \c hx, hy, Nx, Ny Size of the cavity and number of points on each axis.
 * - \c rtype, \c ep point distribution, randomness
 * - \c dt, \c max_iter, \c max_knots time step, time iterations, neighbors
 * - \c c shape parameter for MQ-RBF kernel (\f$c < 0\f$ implies \f$ c = 1/\sqrt{N}\f$).
 * - \c Re Reydolds number.
 * - \c beta_s, beta_v, beta_vb under-relaxation parameters for streamfunction, vorticity and boundary conditions for vorticity.
 * - \c gmres_tol, \c tolerance Tolerance for GMRES, global tolerance
 * - \c Ngrid points for plotting the final solution.
 * - \c ORDER order for the vorticity BC.
 * \par Output
 * - \c xy_knots.dat coordinates of random points;
 * - \c solStream.dat (x,y,u) Streamfunction.
 * - \c solVort.dat (x,y,u) Vorticity.
 * - \c u_line.dat, \c v_line.dat velocity profiles.
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

void vorticityBC(Vec&, Mat&, const Mat&, int*, prec_t);

void writeSolution(const Vec&, const Vec&, const Vec&, prec_t, string, int,
		   Mat &);

void writeToFile_VTK(const Mat&, std::string, std::string);

bool writeMatrix(const Mat&, string);

//==========================================================================
//                            GLOBAL DATA
//
int Nx, Ny, N, NI, NB, NG;
prec_t hx, hy, dx, dy;

int ORDER;

#ifdef WITH_GNUPLOT
GNUplot plotStream, plotVort;
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

    int max_iter, max_knots, rtype, Ngrid;
    prec_t ep, c, Re, beta_s, beta_v, beta_vb, dt;
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
	       >> beta_v    // Under-relaxation coefficient for vorticity
	       >> beta_vb   // Under-relaxation coefficient for boundary vort
	       >> gmres_tol // Tolerance for GMRES
	       >> tolerance // Global tolerance
	       >> Ngrid    // Grid for final results
	       >> ORDER;   // Order for vorticity boundary
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

    std::cout << "\n | Knots generation elapsed time = " << time.toc()
	      << "\n +-----+  ";

    if (c < 0 ) c = 1.0 / sqrt (1.0 * N);

    std::cout << "\n +-----+  "
	      << "\n | c = " << c
	      << "\n +-----+  ";

// ----- Declaration of matrices and vectors
    Mat Gstream(N, N), Gvort(N,N);
    Vec b(N), lstream(N), lvort(N), uvel(N), vvel(N);
    Mat sgrid(Nx, Ny), vgrid(Nx, Ny);      

// ----- Calculate the Gram Matrix for streamfunction

    std::cout << "\n +-----+  "
	      << "\n | Filling the linear system: streamfunction ...  ";
    
    time.tic();
    streamLaplace( Gstream, x, y, 
		   RBF::MQ<prec_t,2>(c),
		   RBF::MQ_2DX<prec_t,2>(c), 
		   RBF::MQ_2DY<prec_t,2>(c), dt);

    std::cout << "\n | Elapsed time = " << time.toc()
	      << "\n +-----+  ";

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
    int point[8];
//
//           6       7
//           |       |
//           v       v   
//       *---*---*---*---*
//       |         *     |
// 1 --> *   *   *       * <-- 3
//       |           *   |
//       *  *    *       *
//       |            *  |
// 0 --> *   *   *  *    * <-- 2
//       |               |
//       *---*---*---*---*
//           ^       ^ 
//           |       |
//           4       5
//
    point[0] = NI + 1;       point[1] = point[0] + side_y - 1;
    point[2] = point[1] + 1; point[3] = point[2] + side_y - 1;
    point[4] = point[3] + 1; point[5] = point[4] + side_x - 1;
    point[6] = point[5] + 1; point[7] = point[6] + side_x - 1;

#ifdef WITH_GNUPLOT
plotStream("set contour");
plotStream("unset surface");
plotStream("set grid");
plotStream("set view 0,0,1.5,1");
plotStream("set size 0.73,1.0");
plotStream("set cntrparam level discrete -0.11, -0.1, -0.09, -0.07, -0.05, -0.03, -0.01, -0.001, -0.0001, -0.00001, 0.0, 0.00001, 0.0001, 0.001, 0.01");
plotVort("set contour");
plotVort("unset surface");
plotVort("set grid");
plotVort("set view 0,0,1.5,1");
plotVort("set size 0.73,1.0");
plotVort("set cntrparam level discrete -5, -4, -3, -2, -1, -0.5, 0, 0.5, 1, 2, 3");
#endif

// ----- Main Loop : BEGIN
    double total_Gvort = 0, total_Svort = 0, total_Sstream = 0;

    while ( (error_v > tolerance) && (t < max_iter) ) {
	t++;

// ----- Calculate the Gram Matrix for the vorticity
	time.tic();
	vortEquation(Gvort, x, y, 
		     RBF::MQ<prec_t,2>(c),
		     RBF::MQ_1DX<prec_t,2>(c),
		     RBF::MQ_1DY<prec_t,2>(c),
		     RBF::MQ_2DX<prec_t,2>(c), 
		     RBF::MQ_2DY<prec_t,2>(c),
		     uvel, vvel, Re, dt);
	total_Gvort += time.toc();

	DiagonalPrec<Mat> preVort(Gvort);

// ----- Boundary conditions for the Vorticity
	vorticityBC(b, vgrid, sgrid, point, beta_vb);

// ----- Solving for the vorticity

	id = 1;
	for(int j = 2; j < Ny; ++j)
	    for(int i = 2; i < Nx; ++i, ++id) 
		b(id) = vgrid(i,j);

	time.tic();
	niter = gmres (Gvort, lvort, b, preVort, N-1 , gmres_tol);
	total_Svort += time.toc();

	std::cout << "\n | " << t <<"\t Vorticity : GMRES Iterations : " 
		  << niter;

	error_v = evalSolution(vgrid, x, y, lvort, c, beta_v, "solVort.dat");
	std::cout << "\n | Vorticity : error = " << error_v;

#ifdef WITH_GNUPLOT
	plotVort("splot [-0.1:1.1][-0.1:1.1][] \"solVort.dat\" w l");
//	std::cout << "\n | Press <1> and then <enter> "; std::cin >> stop;  
#endif

	b = 0;
	id = 1;
	for(int j = 2; j < Ny; ++j)
	    for(int i = 2; i < Nx; ++i, ++id) 
		b(id) =  sgrid(i,j) + dt * vgrid(i,j);

	time.tic();
	niter = gmres (Gstream, lstream, b, preStream, N-1 , gmres_tol);
	total_Sstream += time.toc();

	std::cout << "\n | " << t <<"\t Streamfuntion : GMRES Iterations : " 
		  << niter;
	
	error_s = evalSolution(sgrid, x, y, lstream, c, beta_s, "solStream.dat");
	std::cout << "\n | Streamfunction : error = " << error_s; 

	evalVelocity(uvel, vvel, x, y, lstream, c);

#ifdef WITH_GNUPLOT
	plotStream("splot [-0.1:1.1][-0.1:1.1][] \"solStream.dat\" w l");
//	std::cout << "\n | Press <1> and then <enter> "; std::cin >> stop; 
#endif

    }    
#ifdef WITH_GNUPLOT
    std::cout << "\n | Press <1> and then <enter> "; std::cin >> stop; 
#endif

// Main Loop : END

// ----- Write profiles of the velocity
    prec_t ul, vl, dl = 1.0 / Ngrid;
    prec_t umin = 1, vmin = 1, vmax = 0;
    std::ofstream ufile ("u_line.dat"), vfile ("v_line.dat");
//    ufile << 0.0 << "\t" << 0.0 << "\n";
//    vfile << 0.0 << "\t" << 0.0 << "\n";
    for(int i = 0; i <= Ngrid; ++i) {	
	ul = RBF::eval(0.5, dl*i, x, y, lstream,RBF::MQ_1DY<prec_t,2>(c),0);
	vl =-RBF::eval(dl*i, 0.5, x, y, lstream,RBF::MQ_1DX<prec_t,2>(c),0);
	ufile << ul << "\t" << dl * i << "\n";
	vfile << dl * i << "\t" << vl << "\n";
	if (ul < umin) umin = ul;
	if (vl < vmin) vmin = vl;
	if (vl > vmax) vmax = vl;
    }
//    ufile << 1.0 << "\t" << 1.0 << "\n";
//    vfile << 1.0 << "\t" << 0.0 << "\n";

    ufile.close();
    vfile.close();
    
    std::cout << "\n | Umin (x=0.5) | Vmin (y=0.5) | Vmax (y=0.5) | ";    
    std::cout << "\n | " << umin << " " << vmin << " " << vmax;

// ----- Write the final solution
/* *
    Mat solPrint(Ngrid, Ngrid);
    writeSolution(x, y, lstream, c, "solStream.dat", Ngrid, solPrint);
    writeToFile_VTK(solPrint, "stre.vtk", "streamfunction");
    writeSolution(x, y, lvort, c, "solVort.dat", Ngrid, solPrint);
    writeToFile_VTK(solPrint, "vort.vtk", "vorticity");
/* */
#ifdef WITH_GNUPLOT    
//    plotVort("splot [-0.1:1.1][-0.1:1.1][] \"solVort.dat\" w l");
//    plotStream("splot [-0.1:1.1][-0.1:1.1][] \"solStream.dat\" w l");
//    std::cout << "\n | Press <1> and then <enter> "; std::cin >> stop; 
    plotStream("plot \"u_line.dat\" w lp, \"./Ghia_Bench/u_100_Ghia.dat\" w p");
    plotVort("plot \"v_line.dat\" w lp, \"./Ghia_Bench/v_100_Ghia.dat\" w p");
    std::cout << "\n | Press <1> and then <enter> "; std::cin >> stop; 
#endif

    std::cout << "\n Total time ( Ax = b ) = " << total_time;
    std::cout << "\n +-----+  ";
    std::cout << "\n + Happy finish :-) ";
    std::cout << "\n +-----+ \n\n";    

    cout << "\n\n Gvort = " << total_Gvort
	 << "\n Svort = " << total_Svort 
	 << "\n Sstream = " << total_Sstream
	 << "\n Total  = " << total_Gvort + total_Svort + total_Sstream
	 << "\n";

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
	    G(i,j) = rbf( x(i), y(i), x(j), y(j) ) -
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
    for(int j = 1; j <= N; ++j) 
	for(int i = 1; i <= NI; ++i) {
	    G(i,j) = rbf( x(i), y(i), x(j), y(j) ) -
		dt * ( rbf_2dx ( x(i), y(i), x(j), y(j) ) + 
		       rbf_2dy ( x(i), y(i), x(j), y(j) ) ) +
		dt * Re * ( u(i) * rbf_1dx( x(i), y(i), x(j), y(j) ) +
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

    std::ofstream grid_file (name.c_str());    
    prec_t xd, yd;
    prec_t ddx = 1. / (Nx - 1);
    prec_t ddy = 1. / (Ny - 1);
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
		    prec_t c, string name, int Ngrid, Mat& sol)
{
    prec_t xd, yd;
    prec_t dh = 1. / (Ngrid - 1);
    std::ofstream grid_file (name.c_str());    
    for(int j = 1; j <= Ngrid; ++j) {
	for(int i = 1; i <= Ngrid; ++i) {
	    xd = (i-1) * dh;
	    yd = (j-1) * dh;
	    sol(i,j) = RBF::eval( xd, yd, x, y, lambda, RBF::MQ<prec_t,2>(c));
	    grid_file << xd <<"\t"<< yd << "\t" << sol(i,j) << "\n";
	}
	grid_file << "\n";	    
    }
    grid_file.close();
} 


void vorticityBC(Vec& b, Mat& vgrid, const Mat& sgrid, int* point, prec_t beta)
{
    b = 0;  // Source term

    if (ORDER == 1 ) {
// First order approximation:
    for(int i = point[0], k = 2; i <= point[1]; ++i, ++k) // Edge 1: Dirichlet
	b(i) = -2 * sgrid(2, k) / (dx * dx);

    for(int i = point[2], k = 2; i <= point[3]; ++i, ++k) // Edge 2: Dirichlet
	b(i) = -2 * sgrid(Nx-1, k) / (dx * dx);

    for(int i = point[4], k = 2; i <= point[5]; ++i, ++k) // Edge 3: Dirichlet
	b(i) = -2 * sgrid(k, 2) / (dy * dy);

    for(int i = point[6], k = 2; i <= point[7]; ++i, ++k) // Edge 4: Dirichlet 
	b(i) = -2 * sgrid(k, Ny-1) / (dy * dy) - 2 / dy;  // (sliding wall)

    }

    if (ORDER == 2 ) {    
// Second order approx. (Jensen formula):
    for(int i = point[0], k = 2; i <= point[1]; ++i, ++k) // Edge 1: Dirichlet
	b(i) = -(4 * sgrid(2,k) - 0.5 * sgrid(3,k)) / (dx * dx);
    
    for(int i = point[2], k = 2; i <= point[3]; ++i, ++k) // Edge 2: Dirichlet
	b(i) = -(4 * sgrid(Nx-1, k) - 0.5 * sgrid(Nx-2, k)) / (dx * dx);
    
    for(int i = point[4], k = 2; i <= point[5]; ++i, ++k) // Edge 3: Dirichlet
	b(i) = -(4 * sgrid(k, 2) - 0.5 * sgrid(k, 3)) / (dy * dy);
    
    for(int i = point[6], k = 2; i <= point[7]; ++i, ++k) // Edge 4: Dirichlet 
	b(i) = -(4 * sgrid(k, Ny-1) - 0.5 * sgrid(k, Ny-2)) / (dy*dy) - 3/dy;
    }

    if (ORDER == 3 ) {
// High order approx. see Eurtuk[], Stortkuhl[];
    for(int i = point[0], k = 2; i <= point[1]; ++i, ++k) // Edge 1: Dirichlet
	b(i) = ( -( sgrid(2,k-1) + sgrid(2,k) + sgrid(2,k+1) ) / (3*dx*dx) 
		 -( 0.5 * vgrid(1,k-1) + 0.5 * vgrid(1,k+1) + 
		    0.25 * vgrid(2,k-1) + vgrid(2,k) + 
		    0.25 * vgrid(2,k+1) ) / 9 ) * 9 / 2;
    
    for(int i = point[2], k = 2; i <= point[3]; ++i, ++k) // Edge 2: Dirichlet
	b(i) = ( -( sgrid(Nx-1,k-1) + sgrid(Nx-1,k) + sgrid(Nx-1,k+1) ) 
		 / (3*dx*dx)
		 -( 0.5 * vgrid(Nx,k-1) + 0.5 * vgrid(Nx,k+1) + 
		    0.25 * vgrid(Nx-1,k-1) + vgrid(Nx-1,k) + 
		    0.25 * vgrid(Nx-1,k+1) ) / 9 ) * 9 / 2;
    
    for(int i = point[4], k = 2; i <= point[5]; ++i, ++k) // Edge 3: Dirichlet
	b(i) = ( -( sgrid(k-1,2) + sgrid(k,2) + sgrid(k+1,2) ) / (3*dy*dy) 
		 -( 0.5 * vgrid(k-1,1) + 0.5 * vgrid(k+1,1) + 
		    0.25 * vgrid(k-1,2) + vgrid(k,2) + 
		    0.25 * vgrid(k+1,2) ) / 9 ) * 9 / 2;
    
    for(int i = point[6], k = 2; i <= point[7]; ++i, ++k) // Edge 4: Dirichlet 
	b(i) = ( -1.0 / dy                                // (sliding wall)
		 -( sgrid(k-1,Ny-1) + sgrid(k,Ny-1) + sgrid(k+1,Ny-1) ) 
		 / (3*dy*dy)
		 -( 0.5 * vgrid(k-1,Ny) + 0.5 * vgrid(k+1,Ny) +
		    0.25 * vgrid(k-1,Ny-1) + vgrid(k,Ny-1) + 
		    0.25 * vgrid(k+1,Ny-1) ) / 9 ) * 9 / 2;
    }

// Corners:
// Average between wall neighbors:

    if (ORDER !=3) {
    b(N-3) = 0.5 * ( b(point[0]) + b(point[4]) );    // Bottom-left
    b(N-1) = 0.5 * ( b(point[2]) + b(point[5]) );    // Bottom-right
    b(N-2) = 0.5 * ( b(point[1]) + b(point[6]) );    // Upper-left
    b(N  ) = 0.5 * ( b(point[3]) + b(point[7]) );    // Upper-right
    } else {
// High order approx. see Eurtuk[], Stortkuhl[];
    b(N-3) = -3 * sgrid(2,2) / (dx * dx) 
	- ( 0.5 * vgrid(2,1) + 0.5 * vgrid(1,2) + 
	    0.25 * vgrid(2,2) );
    
    b(N-1) = -3 * sgrid(Nx-1,2) / (dx * dx) 
	- ( 0.5 * vgrid(Nx-1,1) + 0.5 * vgrid(Nx,2) + 
	    0.25 * vgrid(Nx-1,2) );
    
    b(N-2) = -4.5 / dx - 3 * sgrid(2,Ny-1) / (dx*dx)
	- ( 0.5 * vgrid(2,Ny) + 
	    0.5 * vgrid(1,Ny-1) + 
	    0.25 * vgrid(2,Ny-1) );
    
    b(N  ) = -4.5 / dx - 3 * sgrid(Nx-1,Ny-1) / (dx*dx) 
	- ( 0.5 * vgrid(Nx-1,Ny) + 
	    0.5 * vgrid(Nx,Ny-1) + 
	    0.25 * vgrid(Nx-1,Ny-1) );
    }
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
    b(N-2) = beta * b(N-2) + (1 - beta) * vgrid(1,Ny); 
    b(N-1) = beta * b(N-1) + (1 - beta) * vgrid(Nx,1); 
    b(N  ) = beta * b(N  ) + (1 - beta) * vgrid(Nx,Ny);
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

void writeToFile_VTK(const Mat& data, std::string name, std::string variable)
{
    int Nvtk_x = data.numRows();
    int Nvtk_y = data.numCols();
    int Nvtk_z = 1;
    prec_t dx = 1. / (Nvtk_x - 1);
    prec_t dy = 1. / (Nvtk_y - 1);
    prec_t dz = 1.;
    int Ntotal = Nvtk_x * Nvtk_y * Nvtk_z;

    std::ofstream file (name.c_str());        
    file << "# vtk DataFile Version 2.0 \n"
	 << "Lid_driven streamfuntion-vorticity \n"
	 << "ASCII \n"
	 << "DATASET RECTILINEAR_GRID \n" 
	 << "DIMENSIONS " 
	 << Nvtk_x << " " << Nvtk_y << " " << Nvtk_z;

    file << "\nX_COORDINATES " << Nvtk_x << " float \n";
    for(int i = 0; i < Nvtk_x; ++i) {
	file << i * dx << " ";
	if ( !(i % 6) ) file << "\n";
    }

    file << "\nY_COORDINATES " << Nvtk_y << " float \n"; 
    for(int i = 0; i < Nvtk_y; ++i) {
	file << i * dy << " ";
	if ( !(i % 6) ) file << "\n";
    }

    file << "\nZ_COORDINATES " << Nvtk_z << " float \n"; 
    for(int i = 0; i < Nvtk_z; ++i) {
	file << i * dz << " ";
	if ( !(i % 6) ) file << "\n";
    }
	
    file << "\nPOINT_DATA " << Ntotal
	 << "\nSCALARS " << variable << " float"
	 << "\nLOOKUP_TABLE default \n";
    for(int k = 1; k <= Nvtk_z ; ++k) {
	for(int j = 1; j <= Nvtk_y ; ++j) {
	    for(int i = 1; i <= Nvtk_x; ++i) {
		file << data(i,j) << " ";
		if ( !(i % 6) ) file << "\n";
	    }
	}
    }
        
    file.close();

}

/*! 
 ***************************************************************************
 *
 *  \file lid_driven1.cpp
 *  Lid Driven Cavity in 2D.
 *  In this example the time-dependent convection-diffusion equation is 
 *  solved. The equations is written as follows:
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
 * \image html  Xforcedconv.png "Geometry and boundary conditions" width=5cm 
 * \image latex  Xforcedconv.eps "Geometry and boundary conditions" width=5cm 
 *  \par Compiling and running
 *  Modify the variables FLENS, ATLAS and RBF in the file \c rules.in, 
 *  to the path where these libraries are installed. Then type the next 
 *  \verbatim
    % make
    % ./convdiff02 \endverbatim   
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
#include "Knots/RectangleKnots.hpp"

//#include "Solvers/Gauss.hpp" 
#include "Solvers/Krylov/diagonalprec.hpp"
//#include "Solvers/Krylov/identityprec.hpp"
//#include "Solvers/Krylov/acbfprec.hpp"
#include "Solvers/Krylov/gmres.hpp"
#include "RBF/RBF.hpp"
#include "RBF/MQ.hpp"
#include "GNUplot.hpp"

//==========================================================================
//                            FUNCTION PROTOTYPES
//
template<typename RBF, typename RBF_1DY, typename RBF_1DX, 
	 typename RBF_2DX, typename RBF_2DY>
void streamLaplace(Mat&, const Vec&, const Vec&, 
		   RBF_2DX, RBF_2DY, RBF, RBF_1DY, RBF_1DX);

template<typename RBF, typename RBF_1DY, typename RBF_1DX,
	 typename RBF_2DX, typename RBF_2DY>
void vortEquation(Mat&, const Vec&, const Vec&,
		  RBF_2DX, RBF_2DY, RBF, RBF_1DY, RBF_1DX ,
		  const Vec& u, const Vec& v, prec_t);

prec_t evalSolution(Vec&, const Vec&, const Vec&, const Vec&, prec_t, 
		    prec_t, Vec&, string);
void evalVelocity(Vec&, Vec&, const Vec&, const Vec&, 
		  const Vec&, const Vec&, const Vec&, prec_t);

//// FUNTION TO DEBUG, MUST BE ELIMINATED IN THE FINAL VERSION
bool dumpFiles(const Mat&, const Vec&, const Vec&, const Vec&);


//==========================================================================
//                            GLOBAL DATA
//
int Nx, Ny, N, NI, NB, NG;
int Ngrid = 40;
prec_t hx, hy;

GNUplot plotStream, plotVort;
int pausa;

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

    int max_iter, max_knots, rtype;
    prec_t ep, c, dx, dy, Re, beta;

    std::ifstream input_file ("input");
    input_file >> hx   // length in x-axis
	       >> hy   // length in y-axis
	       >> Nx // Boundary points in x-axis
	       >> Ny // Boundary points in y-axis
	       >> rtype
	       >> ep   // Randomness of the interior points (0-0.5)
	       >> max_iter // Number of time steps
	       >> max_knots 
	       >> c
	       >> Re // Reynolds number
	       >> beta;
    input_file.close();
    random_t RT = static_cast<random_t>(rtype);

    dx = hx / (Nx - 1);
    dy = hy / (Ny - 1);

 // ----- Point generation 

    std::cout << "\n +-----+  "
	      << "\n | Calculating knots ... ";  
   
    RectangleKnots<prec_t> rect(hx, Nx, hy, Ny, RT);
//    rect.constructKnotsDriven();
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
    
    RectangleKnots<prec_t> rect_grid(hx, Nx, hy, Ny);
    rect_grid.gridKnots();
    rect_grid.writeToFile("xy_grid.dat");
    Vec xgrid = rect_grid.getKnots(X);
    Vec ygrid = rect_grid.getKnots(Y);
    
// ----- Fill the matrices using MQ and its derivatives.  

    std::cout << "\n +-----+  ";
    std::cout << "\n | Filling the linear system ...  ";
    
    Mat Gstream(N, N), Gvort(N,N);
    Vec b(N);
    Vec stream(N), lstream(N), sgrid(Ngrid * Ngrid);
    Vec vort(N), lvort(N), vgrid(Ngrid * Ngrid);
    Vec uvel(N), vvel(N);
    
    time.tic();
    streamLaplace(Gstream, x, y, 
		  RBF::MQ_2DX<prec_t,2>(c), 
		  RBF::MQ_2DY<prec_t,2>(c),
		  RBF::MQ<prec_t,2>(c),
		  RBF::MQ_1DY<prec_t,2>(c),
		  RBF::MQ_1DX<prec_t,2>(c));
    
    std::cout << "\n | Elapsed time = " << time.toc();

    time.tic();
// ----- Preconditioner construction
    std::cout << "\n | Construction of the preconditioner....";
    DiagonalPrec<Mat> preStream(Gstream);
    std::cout << "\n | Elapsed time = " << time.toc();
    std::cout << "\n +-----+  ";

// ----- Boundary conditions.

    int bi, ei;
    int side_x = Nx - 2, side_y = Ny - 2;
    int point[8];
    point[0] = NI + 1;       point[1] = NI + side_y;
    point[2] = point[1] + 1; point[3] = point[2] + side_y - 1;
    point[4] = point[3] + 1; point[5] = point[4] + side_x - 1;
    point[6] = point[5] + 1; point[7] = point[6] + side_x - 1;
    prec_t s_2, s_3;

// ----- Main Loop

    int niter;
    prec_t tol = 1.0e-8, tolerance = 1.0e-6, error = 0; 

/* */

    plotStream("set contour");
    plotStream("unset surface");
    plotStream("set grid");
    plotStream("set view 0,0,1.5,1");
    plotStream("set size 0.73,1.0");
    plotStream("set cntrparam level discrete -0.1, -0.09, -0.07, -0.05, -0.03, -0.01, -0.001, -0.0001, 0.0");
    //plotStream("set cntrparam level 20");

    plotVort("set contour");
    plotVort("unset surface");
    plotVort("set grid");
    plotVort("set view 0,0,1.5,1");
    plotVort("set size 0.73,1.0");
    plotVort("set cntrparam level discrete -5, -4, -3, -2, -1, 0, 0.5, 1, 2, 3");
    //plotVort("set cntrparam level 20");

    for(int t = 1; t <= max_iter; ++t) 
    {
	niter = gmres (Gstream, lstream, b, preStream, N-1 , tol);    
	std::cout << "\n | " << t <<"\t Streamfuntion : GMRES Iterations : " 
		  << niter;
	
	error = evalSolution(stream, x, y, lstream, c, 
			     beta, sgrid, "solStream.dat");
	std::cout << "\n | Streafunction : error = " << error; 
	
	evalVelocity(uvel, vvel, x, y, xgrid, ygrid, stream, c);
	
	plotStream("splot [-0.1:1.1][-0.1:1.1][] \"solStream.dat\" w l");
	std::cout << "\n | Hit <1> and <enter> "; std::cin >> pausa; 

	vortEquation(Gvort, x, y, 
		     RBF::MQ_2DX<prec_t,2>(c), 
		     RBF::MQ_2DY<prec_t,2>(c),
		     RBF::MQ<prec_t,2>(c),
		     RBF::MQ_1DX<prec_t,2>(c),
		     RBF::MQ_1DY<prec_t,2>(c),
		     uvel, vvel, Re);

// ----- Vorticity boundary condtions
	b = 0;
// Edges 1: Dirichlet
	bi = NI + 1;
	ei = bi + side_y;
	for(int i = bi, ii = 0; i < ei; ++i, ++ii) {	  	  	    
	    s_2 = stream(ii * side_x + 1);
	    s_3 = stream(ii * side_x + 2);	    
//	    b(i) = (-8 * s_2 + s_3) / (2 * dx * dx);	    
	    b(i) = -2 * s_2 / (dx * dx);
	}

// Edge 2: Dirichlet
	bi = ei;      
	ei = bi + side_y;
	for(int i = bi, ii = 1; i < ei; ++i, ++ii) {	  	  	  
	    s_2 = stream(ii * side_x);
	    s_3 = stream(ii * side_x - 1);	    
//	    b(i) = (-8 * s_2 + s_3 ) / (2 * dx * dx);
	    b(i) = -2 * s_2 / (dx * dx);	    
	}

// Edges 3: Dirichlet
	bi = ei;
	ei = bi + side_x;
	for(int i = bi, ii = 1; i < ei; ++i, ++ii) {	  	  
	    s_2 = stream(ii);
	    s_3 = stream(ii + side_x);	    
//	    b(i) = (-8 * s_2 + s_3 ) / (2 * dy * dy);
	    b(i) = -2 * s_2 / (dy * dy);	    
	}

// Edges 4: Dirichlet (sliding wall)
	bi = ei;
	ei = bi + side_x;	
	for(int i = bi, ii = 1; i < ei; ++i, ++ii) {	  	  	  
	    s_2 = stream(ii + side_x * (side_y - 1));
	    s_3 = stream(ii + side_x * (side_y - 2));	    
//	    b(i) = (-8 * s_2 + s_3 ) / (2 * dy * dy) - 3 / dy;
	    b(i) = -2 * s_2 / (dy * dy) - 2 / dy;	    
	}
// Corners:
	b(N-3) = 0.5 * ( b(point[0]) + b(point[4]) );
	b(N-2) = 0.5 * ( b(point[1]) + b(point[6]) );
	b(N-1) = 0.5 * ( b(point[2]) + b(point[5]) );
	b(N  ) = 0.5 * ( b(point[3]) + b(point[7]) );
/*
	b(N-3) = 0.5 * ( vort(point[0]) + vort(point[4]) );
	b(N-2) = 0.5 * ( vort(point[1]) + vort(point[6]) );
	b(N-1) = 0.5 * ( vort(point[2]) + vort(point[5]) ) ;
	b(N  ) = 0.5 * ( vort(point[3]) + vort(point[7]) );
*/
/* */
	std::ofstream b_file ("b1.dat");    
	for(int i = NI + 1; i <= N; ++i)
	    b_file << x(i) << "\t" << y(i) << "\t" << b(i) << "\n";
	b_file.close();
/* */
	DiagonalPrec<Mat> preVort(Gvort);
	niter = gmres (Gvort, lvort, b, preVort, N-1 , tol);    
	std::cout << "\n | " << t <<"\t Vorticity : GMRES Iterations : " 
		  << niter;
	error = evalSolution(vort, x, y, lvort, c, 
			     beta, vgrid, "solVort.dat");
	std::cout << "\n | Vorticity : error = " << error;	

	plotVort("splot [-0.1:1.1][-0.1:1.1][] \"solVort.dat\" w l");
	std::cout << "\n | Hit <1> and <enter> "; std::cin >> pausa;  

	b = 0;
	for(int i = 1; i <= NI; ++i) 
	    b(i) = -vort(i);			    
    }
/* */

    std::cout << "\n | Hit <1> and <enter> "; std::cin >> pausa;  

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
template<typename RBF, typename RBF_1DY, typename RBF_1DX,
	 typename RBF_2DX, typename RBF_2DY>
void streamLaplace(Mat& G, const Vec& x, const Vec& y,
		   RBF_2DX rbf_2dx, RBF_2DY rbf_2dy, 
		   RBF rbf, RBF_1DY rbf_1dy, RBF_1DX rbf_1dx)
//		   const Vec& vort)
{
// ----- Wl matrix

    for(int j = 1; j <= N; ++j) 
	for(int i = 1; i <= NI; ++i) {
	    G(i,j) = rbf_2dx ( x(i), y(i), x(j), y(j) ) + 
	      rbf_2dy ( x(i), y(i), x(j), y(j) );// + vort(i);
	}

// ----- Wb matrix
    for(int j = 1; j <= N; ++j) 
	for(int i = NI+1; i <= N; ++i)
	    G(i,j) = rbf ( x(i), y(i), x(j), y(j) );   
}
//--------------------------------------------------------------------------
// Fill the matrix of the system
//
template<typename RBF, typename RBF_1DY, typename RBF_1DX,
	 typename RBF_2DX, typename RBF_2DY>
void vortEquation(Mat& G, const Vec& x, const Vec& y,
		  RBF_2DX rbf_2dx, RBF_2DY rbf_2dy, 
		  RBF rbf, RBF_1DY rbf_1dy, RBF_1DX rbf_1dx,
		  const Vec& u, const Vec& v, prec_t Re)
{
// ----- Wl matrix
    for(int j = 1; j <= N; ++j) 
	for(int i = 1; i <= NI; ++i) {
	    G(i,j) = rbf_2dx ( x(i), y(i), x(j), y(j) ) + 
		     rbf_2dy ( x(i), y(i), x(j), y(j) ) -
		Re * ( u(i) * rbf_1dx( x(i), y(i), x(j), y(j) ) +
		       v(i) * rbf_1dy( x(i), y(i), x(j), y(j) ) );
	}

// ----- Wb matrix
    for(int j = 1; j <= N; ++j) 
	for(int i = NI+1; i <= N; ++i)
	    G(i,j) = rbf ( x(i), y(i), x(j), y(j) );   
}


//--------------------------------------------------------------------------
// Evaluate the velocity
//
void evalVelocity(Vec& u, Vec& v, const Vec& x, const Vec& y, 
		  const Vec& xg, const Vec& yg,
		  const Vec& stream, prec_t c)
{
    std::ofstream u_file ("uvel.dat");    
    std::ofstream v_file ("vvel.dat");    
    std::ofstream s_file ("stemp.dat");    

    int p;
    /* */
    prec_t dx = hx / (Nx - 1);
    prec_t dy = hy / (Ny - 1);
    prec_t s_m1, s_p1;
    int id = Nx + 2, ii = 1;
    for(int j = 2; j <= Ny-1; ++j) {
	for(int i = 2; i <= Nx-1; ++i, ++id, ++ii) {
	    s_m1 = stream(id - Nx);
	    s_p1 = stream(id + Nx);
	    u(ii) = 0.5 * (s_p1 - s_m1) / dy;

//	    cout << " u = " << s_m1 << "\t" << s_p1 << "\t" << u(ii) << endl;

	    s_m1 = stream(id-1);
	    s_p1 = stream(id+1); 
	    v(ii) = -0.5 * (s_p1 - s_m1) / dx;	

//	    cout << " v = "<< s_m1 << "\t" << s_p1 << "\t" << v(ii) << endl;
	}  
	id += 2;
//	u_file << "\n";
//	v_file << "\n";
    }

//    cin >> p;

    for(int i = 1; i <= NI; ++i) {
	s_file << x(i) << "\t" << y(i) << "\t" << stream(i) << "\n";
	u_file << x(i) << "\t" << y(i) << "\t" << u(i) << "\n";
	v_file << x(i) << "\t" << y(i) << "\t" << v(i) << "\n";
    }

    /* *
    for(int i = 1; i <= NI; ++i) {
	u(i) = RBF::eval(x(i), y(i), x, y, lstream,RBF::MQ_1DY<prec_t,2>(c),0);
	v(i) =-RBF::eval(x(i), y(i), x, y, lstream,RBF::MQ_1DX<prec_t,2>(c),0);
	u_file << x(i) << "\t" << y(i) << "\t" << u(i) << "\n";
	v_file << x(i) << "\t" << y(i) << "\t" << v(i) << "\n";
    }
/* */
    u_file.close();
    v_file.close();
    s_file.close();
}

//--------------------------------------------------------------------------
// Evaluate the numerical solution on a mesh of Nx by Ny
//

prec_t evalSolution(Vec& f, const Vec& x, const Vec& y, const Vec& lambda,
		    prec_t c, 
		    prec_t beta, Vec& fgrid, string name)
{
    std::ofstream b_file ("solRand.dat");    
    prec_t fnew, error = 0;
    for(int i = 1; i <= NI; ++i) {	  
	fnew = RBF::eval(x(i), y(i), x, y, lambda, RBF::MQ<prec_t,2>(c),0);
	fnew = beta * fnew + (1 -  beta) * f(i);
	error += fabs( fnew - f(i) );
	f(i) = fnew;
	b_file << x(i) << "\t" << y(i) << "\t" << f(i) << "\n";
    }
    b_file.close();

    int nb = 1, id = 1;
    std::ofstream grid_file (name.c_str());    
    prec_t ddx = 1. / Ngrid;
    prec_t ddy = 1. / Ngrid;
    for(int i = nb; i <= Ngrid - nb; ++i) {
      for(int j = nb; j <= Ngrid - nb; ++j, ++id) { 
	    fnew = RBF::eval( i * ddx, j * ddy, x, y, lambda, 
			      RBF::MQ<prec_t,2>(c), 0 );
	    fnew = beta * fnew + (1 -  beta) * fgrid(id);
	    grid_file << i * ddx << "\t" << j * ddy << "\t" << fnew  << "\n";
	    fgrid(id) = fnew;
      }
      grid_file << "\n";
    }
    grid_file.close();

    return error;
} 

/* *
bool dumpFiles(const Mat& G, const Vec& u, const Vec& x, const Vec& y) 
{
    std::cout << "\n +-----+  ";
    
    std::ofstream WLI_file ("WLI.dat");
    for (int j = 1; j <= NI; ++j){
	for (int i = 1; i <= NI; ++i) 
	    WLI_file << i << " " << j << " " << G(i,j) << std::endl;
	WLI_file << std::endl;
    }

    std::cout << "\n | Hit <1> and <enter> "; std::cin >> pausa;  
    plotter("splot \"WLI.dat\" w l");

    std::ofstream WLG_file ("WLG.dat");
    for (int j = NI+1; j <= NI+NG; ++j){
	for (int i = 1; i <= NI; ++i) 
	    WLG_file << i << " " << j << " " << G(i,j) << std::endl;
	WLG_file << std::endl;
    }
    std::cout << "\n | Hit <1> and <enter> "; std::cin >> pausa;  
    plotter("splot \"WLG.dat\" w l");

    std::ofstream WLB_file ("WLB.dat");
    for (int j = NI+1+NG; j <= N; ++j){
	for (int i = 1; i <= NI; ++i) 
	    WLB_file << i << " " << j << " " << G(i,j) << std::endl;
	WLB_file << std::endl;
    }

    std::cout << "\n | Hit <1> and <enter> "; std::cin >> pausa;  
    plotter("splot \"WLB.dat\" w l");


    std::ofstream WBIG_file ("WBIG.dat");
    for (int j = 1; j <= N; ++j){
	for (int i = NI+1; i <= NI+NG; ++i) 
	    WBIG_file << i << " " << j << " " << G(i,j) << std::endl;
	WBIG_file << std::endl;
    }

    std::cout << "\n | Hit <1> and <enter> "; std::cin >> pausa;  
    plotter("splot \"WBIG.dat\" w l");

    std::ofstream WBB_file ("WBB.dat");
    for (int j = 1; j <= N; ++j){
	for (int i = NI+1+NG; i <= N; ++i) 
	    WBB_file << i << " " << j << " " << G(i,j) << std::endl;
	WBB_file << std::endl;
    }

    std::cout << "\n | Hit <1> and <enter> "; std::cin >> pausa;  
    plotter("splot \"WBB.dat\" w l");

    WLI_file.close();
    WLG_file.close();
    WLB_file.close();
    WBIG_file.close();
    WBB_file.close();

    std::ofstream UI_file ("UI.dat");
    for(int i = 1; i <= NI; ++i) {
	UI_file << x(i) << " " << y(i) << " " << u(i) << "\n";
    }

    std::ofstream UG_file ("UG.dat");
    for(int i = NI+1; i <= NI+NG; ++i) {
	UG_file << x(i) << " " << y(i) << " " << u(i) << "\n";
    }

    std::ofstream UB_file ("UB.dat");
    for(int i = NI+1+NG; i <= N; ++i) {
	UB_file << x(i) << " " << y(i) << " " << u(i) << "\n";
    }
    UI_file.close();
    UG_file.close();
    UB_file.close();


    std::ofstream lambda_file ("lambda.dat");
    for(int i = 1; i <= NI; ++i) {
	lambda_file << x(i) << "\t" << y(i) << "\t" << lambda(i) << std::endl;
	if ( !(i % (Nx-1)) ) lambda_file << "\n";
    }

return 0;
}
/* */

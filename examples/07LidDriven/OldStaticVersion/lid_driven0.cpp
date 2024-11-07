/*! 
 ***************************************************************************
 *  \file lid_driven.cpp
 *  Lid Driven Cavity in 2D.
 *  In this example the lid-driven cavity is solved in the vorticity -
 *  stream function formulation. The equation to solve are:
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
template<typename RBF, typename RBF_2DX, typename RBF_2DY>
void streamLaplace(Mat&, const Vec&, const Vec&, RBF, RBF_2DX, RBF_2DY);

template<typename RBF, 
	 typename RBF_1DX, typename RBF_1DY,
	 typename RBF_2DX, typename RBF_2DY>
void vortEquation(Mat&, const Vec&, const Vec&, RBF,
		  RBF_1DX, RBF_1DY, RBF_2DX, RBF_2DY, 
		  const Vec& u, const Vec& v, prec_t);

void evalVelocity(Vec&, Vec&, const Vec&, const Vec&, const Mat&, const Mat&,
		  const Vec&, prec_t);

prec_t evalSolution(Mat&, const Vec&, const Vec&, const Vec&, prec_t, 
		    prec_t, Mat&, string);

//==========================================================================
//                            GLOBAL DATA
//
int Nx, Ny, N, NI, NB, NG;
int Ngrid = 40; //<--- deprecated
int stop;
prec_t hx, hy, dx, dy;

GNUplot plotStream, plotVort;


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

    int max_iter, max_knots, rtype;
    prec_t ep, c, Re, beta;

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

// ----- Fill the matrices using MQ and its derivatives.  

    std::cout << "\n +-----+  ";
    std::cout << "\n | Filling the linear system ...  ";
    
    Mat Gstream(N, N), Gvort(N,N);
    Vec b(N), lstream(N), lvort(N);
    Mat sgrid(Nx, Ny), sout(Ngrid, Ngrid);
    Mat vgrid(Nx, Ny), vout(Ngrid, Ngrid);
    Vec uvel(N), vvel(N);

    time.tic();
    streamLaplace( Gstream, x, y, RBF::MQ<prec_t,2>(c),
		   RBF::MQ_2DX<prec_t,2>(c), RBF::MQ_2DY<prec_t,2>(c) );
    
    std::cout << "\n | Elapsed time = " << time.toc();

    time.tic();
// ----- Preconditioner construction
    std::cout << "\n | Construction of the preconditioner....";
    ACBFPrec<Mat, RectangleKnots<prec_t> > preStream(Gstream);
    preStream.construct(rect, max_knots, 8);
//    DiagonalPrec<Mat> preStream(Gstream);
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

    int niter, id;
    prec_t tol = 1.0e-8, tolerance = 1.0e-6, error = 0; 

#ifdef WITH_GNUPLOT
/* */
    plotStream("set contour");
    plotStream("unset surface");
    plotStream("set grid");
    plotStream("set view 0,0,1.5,1");
    plotStream("set size 0.73,1.0");
    plotStream("set cntrparam level discrete -0.1, -0.09, -0.07, -0.05, -0.03, -0.01, -0.001, -0.0001, -0.00001, 0.0");
//    plotStream("set cntrparam level 40");

    plotVort("set contour");
    plotVort("unset surface");
    plotVort("set grid");
    plotVort("set view 0,0,1.5,1");
    plotVort("set size 0.73,1.0");
    plotVort("set cntrparam level discrete -5, -4, -3, -2, -1, -0.5, 0, 0.5, 1, 2, 3");
//    plotVort("set cntrparam level 40");
#endif


/* *
    vortEquation(Gvort, x, y, 
		 RBF::MQ<prec_t,2>(c),
		 RBF::MQ_1DX<prec_t,2>(c),
		 RBF::MQ_1DY<prec_t,2>(c),
		 RBF::MQ_2DX<prec_t,2>(c), 
		 RBF::MQ_2DY<prec_t,2>(c),
		 uvel, vvel, Re);
    ACBFPrec<Mat, RectangleKnots<prec_t> > preVort(Gvort);
    preVort.construct(rect, max_knots, 0);
//    DiagonalPrec<Mat> preVort(Gvort);
/* */
    for(int t = 1; t <= max_iter; ++t) 
    {
	niter = gmres (Gstream, lstream, b, preStream, N-1 , tol);    
	std::cout << "\n | " << t <<"\t Streamfuntion : GMRES Iterations : " 
		  << niter;
	
	error = evalSolution(sgrid, x, y, lstream, c, beta, sout, 
			     "solStream.dat");
	std::cout << "\n | Streamfunction : error = " << error; 

	evalVelocity(uvel, vvel, x, y, sgrid, vgrid, 
		     lstream, c);

#ifdef WITH_GNUPLOT
	plotStream("splot [-0.1:1.1][-0.1:1.1][] \"solStream.dat\" w l");
//	std::cout << "\n | Hit <1> and <enter> "; std::cin >> stop; 
#endif

	vortEquation(Gvort, x, y, 
		     RBF::MQ<prec_t,2>(c),
		     RBF::MQ_1DX<prec_t,2>(c),
		     RBF::MQ_1DY<prec_t,2>(c),
		     RBF::MQ_2DX<prec_t,2>(c), 
		     RBF::MQ_2DY<prec_t,2>(c),
		     uvel, vvel, Re);

// ----- Vorticity boundary condtions
	b = 0;
// Edges 1: Dirichlet
	bi = NI + 1;
	ei = bi + side_y;
	for(int i = bi, k = 2; i < ei; ++i, ++k) {	  	  	    
	    b(i) = -0.2 * (4 * sgrid(2,k) - 0.5 * sgrid(3,k)) / (dx * dx);
//	    b(i) = -2 * sgrid(2, k) / (dx * dx);
//	    b(i) = - sgrid(3, k) / (2 * dx * dx);
/* *
	    b(i) = ( -( sgrid(2,k-1) + 
			sgrid(2,k) + 
			sgrid(2,k+1) )
		     / (3*dx*dx) 
		     -( 0.5 * vgrid(1,k-1) + 
			0.5 * vgrid(1,k+1) + 
			0.25 * vgrid(2,k-1) + 
			vgrid(2,k) + 
			0.25 * vgrid(2,k+1)) 
		     / 9 
		) * 9 / 2;
/* */
	}

// Edge 2: Dirichlet
	bi = ei;      
	ei = bi + side_y;
	for(int i = bi, k = 2; i < ei; ++i, ++k) {
	    b(i) = -0.2 * (4 * sgrid(2, k) - 0.5 * sgrid(3, k)) / (dx * dx);
//	    b(i) = -2 * sgrid(Nx-1, k) / (dx * dx);
//	    b(i) = - sgrid(Nx-2, k) / (2 * dx * dx);
/* *
	    b(i) = ( -( sgrid(Nx-1,k-1) + 
			sgrid(Nx-1,k) + 
			sgrid(Nx-1,k+1) ) 
		     / (3*dx*dx) 
		     -( 0.5 * vgrid(Nx,k-1) + 
			0.5 * vgrid(Nx,k+1) + 
		        0.25 * vgrid(Nx-1,k-1) + 
			vgrid(Nx-1,k) + 
		        0.25 * vgrid(Nx-1,k+1) ) 
		     / 9 
		) * 9 / 2;
/* */
	}

// Edges 3: Dirichlet
	bi = ei;
	ei = bi + side_x;
	for(int i = bi, k = 2; i < ei; ++i, ++k) {	  	  
	    b(i) = -0.2 * (4 * sgrid(k, 2) - 0.5 * sgrid(k, 3)) / (dy * dy);
//	    b(i) = -2 * sgrid(k, 2) / (dy * dy);	    
//	    b(i) = - sgrid(k, 3) / (2 * dy * dy);	    
/* *
	    b(i) = ( -( sgrid(k-1,2) + 
			sgrid(k,2) + 
		        sgrid(k+1,2) ) 
		     / (3*dy*dy) 
		     -(
			 0.5 * vgrid(k-1,1) + 
			 0.5 * vgrid(k+1,1) + 
			 0.25 * vgrid(k-1,2) + 
			 vgrid(k,2) + 
			 0.25 * vgrid(k+1,2)
			 )
		     / 9 
		) * 9 / 2;
/* */
	}

// Edges 4: Dirichlet (sliding wall)
	bi = ei;
	ei = bi + side_x;	
	for(int i = bi, k = 2; i < ei; ++i, ++k) {	  	  	  
	    b(i) = -0.2 * ( (4 * sgrid(k, 2) - 0.5 * sgrid(k, 3)) / (dy*dy) + 3/dy );
//	    b(i) = -2 * sgrid(k, Ny-1) / (dy * dy) - 2 / dy;
//	    b(i) = - sgrid(k, Ny-2) / (2 * dy * dy) - 1 / (2 * dy);
/* *
	    b(i) = ( -1.0 / dy 
		     -( sgrid(k-1,Ny-1) + 
			sgrid(k,Ny-1) + 
			sgrid(k+1,Ny-1) ) 
		     / (3*dy*dy)
		     -( 
			 0.5 * vgrid(k-1,Ny) + 
			 0.5 * vgrid(k+1,Ny) +
			 0.25 * vgrid(k-1,Ny-1) + 
			 vgrid(k,Ny-1) + 
			 0.25 * vgrid(k+1,Ny-1) 
			 ) 
		     / 9 
		) * 9 / 2;
/* */
	}
// Corners:
/* *
	b(N-3) = 0.5 * ( b(point[0]) + b(point[4]) );
	b(N-2) = 0.5 * ( b(point[1]) + b(point[6]) );
	b(N-1) = 0.5 * ( b(point[2]) + b(point[5]) );
	b(N  ) = 0.5 * ( b(point[3]) + b(point[7]) );
/* */
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
/* */
	bi = NI + 1; ei = bi + side_y;
	for(int i = bi, k = 2; i < ei; ++i, ++k)
	    vgrid(1, k) = b(i);
	bi = ei; ei = bi + side_y;
	for(int i = bi, k = 2; i < ei; ++i, ++k) 
	    vgrid(Nx, k) = b(i);
	bi = ei; ei = bi + side_x;
	for(int i = bi, k = 2; i < ei; ++i, ++k)
	    vgrid(k, 1) = b(i);
	bi = ei; ei = bi + side_x;
	for(int i = bi, k = 2; i < ei; ++i, ++k)
	    vgrid(k, Ny) = b(i);	
	vgrid(1 , 1 ) = b(N-3);
	vgrid(1 , Ny) = b(N-2);
	vgrid(Nx, 1 ) = b(N-1);
	vgrid(Nx, Ny) = b(N  );

/* */
	std::ofstream b_file ("b0.dat");    
	for(int i = 1; i <= N; ++i)
	    b_file << x(i) << "\t" << y(i) << "\t" << b(i) << "\n";
	b_file.close();	
/* */
//	ACBFPrec<Mat, RectangleKnots<prec_t> > preVort(Gvort);
//	preVort.construct(rect, max_knots, 8);
	DiagonalPrec<Mat> preVort(Gvort);
	niter = gmres (Gvort, lvort, b, preVort, N-1 , tol);    
	std::cout << "\n | " << t <<"\t Vorticity : GMRES Iterations : " 
		  << niter;

	error = evalSolution(vgrid, x, y, lvort, c, beta, vout, "solVort.dat");
	std::cout << "\n | Vorticity : error = " << error;	

#ifdef WITH_GNUPLOT
	plotVort("splot [-0.1:1.1][-0.1:1.1][] \"solVort.dat\" w l");
//	std::cout << "\n | Hit <1> and <enter> "; std::cin >> stop;  
#endif

	b = 0;
	id = 1;
	for(int j = 2; j < Ny; ++j)
	    for(int i = 2; i < Nx; ++i, ++id) 
		b(id) = -vgrid(i,j);

    }
/* *

    std::ofstream sgrid_file ("solStream.dat");    
    std::ofstream vgrid_file ("solVort.dat");    
    prec_t xd, yd;
    prec_t ddx = 1. / (Nx - 1);
    prec_t ddy = 1. / (Ny - 1);
    id = 1;
    for(int j = 1; j <= Ny; ++j) { 
	for(int i = 1; i <= Nx; ++i, ++id) {
	    xd = (i-1) * ddx;
	    yd = (j-1) * ddy;
	    sgrid_file << xd <<"\t"<< yd << "\t" << sgrid(i,j) <<"\n";
	    vgrid_file << xd <<"\t"<< yd << "\t" << vgrid(i,j) <<"\n";
	}
	sgrid_file << "\n";	    
	vgrid_file << "\n";	    
    }
    sgrid_file.close();
    vgrid_file.close();
/* */


    std::cout << "\n Total time ( Ax = b ) = " << total_time;
    std::cout << "\n +-----+  ";
    std::cout << "\n + Happy finish :-) ";
    std::cout << "\n +-----+ \n\n";    

    std::cin >> stop;
    return 0;
}

//==========================================================================
//                            FUNCTION DEFINITIONS
//--------------------------------------------------------------------------
// Fill the matrix of the system
//
template<typename RBF, typename RBF_2DX, typename RBF_2DY>
void streamLaplace(Mat& G, const Vec& x, const Vec& y, RBF rbf,
		   RBF_2DX rbf_2dx, RBF_2DY rbf_2dy)
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
template<typename RBF, typename RBF_1DX, typename RBF_1DY,
	 typename RBF_2DX, typename RBF_2DY>
void vortEquation(Mat& G, const Vec& x, const Vec& y,
		  RBF rbf, RBF_1DX rbf_1dx, RBF_1DY rbf_1dy, 
		  RBF_2DX rbf_2dx, RBF_2DY rbf_2dy, 
		  const Vec& u, const Vec& v, prec_t Re)
{
/*
    std::ofstream u_m ("umat.dat");    
    std::ofstream v_m ("vmat.dat");    
    for(int i = 1; i <= NI; ++i) {
	u_m << x(i) << "\t" << y(i) << "\t" << u(i) << "\n";
	v_m << x(i) << "\t" << y(i) << "\t" << v(i) << "\n";
    }
    u_m.close();
    v_m.close();
    cout << "\n | Re = " <<  Re;
*/
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
		  const Mat& sgrid, const Mat& vgrid,
		  const Vec& lstream, prec_t c)
{
/* *
    std::ofstream u_file ("uvel.dat");    
    std::ofstream v_file ("vvel.dat");    
/* */
    
/* *
    int id = 1;
    for(int j = 2; j <= Ny-1; ++j) {
	for(int i = 2; i <= Nx-1; ++i, ++id) {
	    u(id) = (sgrid(i,j+1) - sgrid(i,j-1)) / (3 * dy) +
		(sgrid(i+1,j+1) + sgrid(i-1,j+1) -
		 sgrid(i-1,j-1) - sgrid(i+1,j-1)) / (12 * dy) +
		dy * (vgrid(i,j+1) - vgrid(i,j-1)) / 12 ;

	    v(id) = (sgrid(i-1,j) - sgrid(i+1,j)) / (3 * dx) -
		(sgrid(i+1,j+1) - sgrid(i-1,j+1) -
		 sgrid(i-1,j-1) + sgrid(i+1,j-1)) / (12 * dy) +
		dy * (vgrid(i-1,j) - vgrid(i+1,j)) / 12 ;

//	    u(id) = 0.5 * ( sgrid(i,j+1) - sgrid(i,j-1) ) / dy;
//	    v(id) = -0.5 * ( sgrid(i+1,j) - sgrid(i-1,j) ) / dx;	

	    u_file << x(id) << "\t" << y(id) << "\t" << u(id) << "\n";
	    v_file << x(id) << "\t" << y(id) << "\t" << v(id) << "\n";
	}  
	u_file << "\n";
	v_file << "\n";
    }
    
/* */
    for(int i = 1; i <= NI; ++i) {
	u(i) = RBF::eval(x(i), y(i), x, y, lstream,RBF::MQ_1DY<prec_t,2>(c),0);
	v(i) =-RBF::eval(x(i), y(i), x, y, lstream,RBF::MQ_1DX<prec_t,2>(c),0);
//	u_file << x(i) << "\t" << y(i) << "\t" << u(i) << "\n";
//	v_file << x(i) << "\t" << y(i) << "\t" << v(i) << "\n";
    }
/* *
    u_file.close();
    v_file.close();
/* */
}

//--------------------------------------------------------------------------
// Evaluate the numerical solution on a mesh of Nx by Ny
//

prec_t evalSolution(Mat& fgrid, const Vec& x, const Vec& y, const Vec& lambda,
		    prec_t c, prec_t beta, Mat& fout, string name)
{

/* *
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

/* */

    std::ofstream grid_file (name.c_str());    

    prec_t fnew, error = 0, xp, yp;
    int id = 1;
    for(int j = 2; j < Ny; ++j) { 
	for(int i = 2; i < Nx; ++i, ++id) {
	    fnew = RBF::eval( x(id), y(id), x, y, lambda, 
			      RBF::MQ<prec_t,2>(c), 0 );
	    fnew = beta * fnew + (1 -  beta) * fgrid(i, j);	    
	    error += fabs( fnew - fgrid(i, j) );
	    fgrid(i, j) = fnew;
	}    
    }

/* */
    prec_t xd, yd;
    prec_t ddx = 1. / (Nx - 1);
    prec_t ddy = 1. / (Ny - 1);

/* */
    id = 1;
    for(int j = 1; j <= Ny; ++j) { 
	for(int i = 1; i <= Nx; ++i, ++id) {
	    xd = (i-1) * ddx;
	    yd = (j-1) * ddy;
	    grid_file << xd <<"\t"<< yd << "\t" << fgrid(i,j) <<"\n";
	}
	grid_file << "\n";	    
    }

/* *
    for(int j = 1; j <= Ngrid; ++j) {
	for(int i = 1; i <= Ngrid; ++i) {
	    xd = (i-1) * dh;
	    yd = (j-1) * dh;
	    fnew = RBF::eval( xd, yd, x, y, lambda, RBF::MQ<prec_t,2>(c), 0 );
	    fnew = beta * fnew + (1 -  beta) * fout(i, j);
	    fout(i, j) = fnew;		    
	    grid_file << xd <<"\t"<< yd << "\t" << fout(i,j) <<"\n";
	}
	grid_file << "\n";	    
    }
/* */
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

    std::cout << "\n | Hit <1> and <enter> "; std::cin >> stop;  
    plotter("splot \"WLI.dat\" w l");

    std::ofstream WLG_file ("WLG.dat");
    for (int j = NI+1; j <= NI+NG; ++j){
	for (int i = 1; i <= NI; ++i) 
	    WLG_file << i << " " << j << " " << G(i,j) << std::endl;
	WLG_file << std::endl;
    }
    std::cout << "\n | Hit <1> and <enter> "; std::cin >> stop;  
    plotter("splot \"WLG.dat\" w l");

    std::ofstream WLB_file ("WLB.dat");
    for (int j = NI+1+NG; j <= N; ++j){
	for (int i = 1; i <= NI; ++i) 
	    WLB_file << i << " " << j << " " << G(i,j) << std::endl;
	WLB_file << std::endl;
    }

    std::cout << "\n | Hit <1> and <enter> "; std::cin >> stop;  
    plotter("splot \"WLB.dat\" w l");


    std::ofstream WBIG_file ("WBIG.dat");
    for (int j = 1; j <= N; ++j){
	for (int i = NI+1; i <= NI+NG; ++i) 
	    WBIG_file << i << " " << j << " " << G(i,j) << std::endl;
	WBIG_file << std::endl;
    }

    std::cout << "\n | Hit <1> and <enter> "; std::cin >> stop;  
    plotter("splot \"WBIG.dat\" w l");

    std::ofstream WBB_file ("WBB.dat");
    for (int j = 1; j <= N; ++j){
	for (int i = NI+1+NG; i <= N; ++i) 
	    WBB_file << i << " " << j << " " << G(i,j) << std::endl;
	WBB_file << std::endl;
    }

    std::cout << "\n | Hit <1> and <enter> "; std::cin >> stop;  
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

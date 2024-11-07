/*! 
 ***************************************************************************
 *
 *  \file lid_driven_par.cpp
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

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>

#include "Traits.hpp"
#include "Knots/OVRectangleKnots.hpp"
#include "DD/CartComm.hpp"
#include "RBF/RBF.hpp"
#include "RBF/MQ.hpp"
#include "Solvers/Krylov/diagonalprec.hpp"
#include "Solvers/Krylov/acbfprec.hpp"
#include "Solvers/Krylov/gmres.hpp"

#include "GNUplot.hpp"

//==========================================================================
//                            FUNCTION PROTOTYPES
//

bool read_data_and_Bcast(CartComm&);
void print_info();

template<typename RBF, typename RBF_2DX, typename RBF_2DY>
void streamLaplace(Mat&, const Vec&, const Vec&, RBF,
		   RBF_2DX, RBF_2DY);

template<typename RBF, typename RBF_1DX, typename RBF_1DY,
	 typename RBF_2DX, typename RBF_2DY>
void vortEquation(Mat&, const Vec&, const Vec&, RBF,
		  RBF_1DX, RBF_1DY, RBF_2DX, RBF_2DY, 
		  const Vec& u, const Vec& v, prec_t);

prec_t evalSolution(Mat&, const Vec&, const Vec&, const Vec&, prec_t, 
		    prec_t, string);
void evalVelocity(Vec&, Vec&, const Vec&, const Vec&, const Vec&, prec_t);

prec_t exchangeInfo(Mat& , const Vec&, const Vec&, const Vec&, 
		    prec_t, CartComm&, Vec&, Vec&, Vec&, Vec&);


//==========================================================================
//                            GLOBAL DATA
//
// Data from "input" file
prec_t hx, hy, over_x, over_y, ep, Re, tol_gmres, tol_ddm, beta;
int Nx, Ny, rtype, max_iter, Ngx, Ngy, max_knots, max_ddm;

// Other useful global information
prec_t sx, sy, dx, dy, c;
int N, NI, NB;
int size, rank, I, J, NP_I, NP_J;
std::string filename, istr, jstr;
std::ostringstream inum, jnum;

#ifdef WITH_GNUPLOT
GNUplot plotStream, plotVort;
#endif
int pausa;

//==========================================================================
//                            MAIN FUNCTION
//
int main( int argc, char * argv[])
{
    MPI::Init(argc, argv);
    rank = MPI::COMM_WORLD.Get_rank();       
    size = MPI::COMM_WORLD.Get_size();       

// ----- Construction of a Cartesian topology
    CartComm cart(argc, argv, rank);    
    I    = cart.get_I();
    J    = cart.get_J();
    NP_I = cart.getNumProc_I();
    NP_J = cart.getNumProc_J();
    const int *neighbor = cart.getNeighbors();

// ----- istr and jstr are used to filenames construction
    inum.width(1); inum.fill('0'); inum << I; istr = inum.str();
    jnum.width(1); jnum.fill('0'); jnum << J; jstr = jnum.str(); 

// ----- Reading and broadcasting info from the input file
    read_data_and_Bcast(cart); 

// ----- Point generation 
    random_t RT = static_cast<random_t>(rtype);
    OVRectangleKnots<prec_t> rect(hx, Nx, hy, Ny, over_x, over_y, cart, RT);
    rect.setRandomness(ep);
    rect.constructKnots();
    rect.constructOverlapping(cart);    
    rect.writeToFile("xyz_", I, J);

    Vec x  = rect.getKnots(X);
    Vec y  = rect.getKnots(Y);
    c      = rect.getShapeParameter();
    N      = rect.getTotalKnots();
    NI     = rect.getInteriorKnots();
    hx     = rect.getLx();            // Size of the subdomain
    hy     = rect.getLy();
    Nx     = rect.getNx();            // Number of points in this subdomain.
    Ny     = rect.getNy();
    sx     = rect.getShift(X);        // Shift for this subdomain.
    sy     = rect.getShift(Y);
    dx     = rect.getDx();            // Size of the "grid"
    dy     = rect.getDy();
    over_x = rect.getOverlap(X);      // Overlap size for this subdomain.
    over_y = rect.getOverlap(Y);

    Vec left  = rect.getArrayL();     // Overlapping arrays
    Vec right = rect.getArrayR();
    Vec up    = rect.getArrayU();
    Vec down  = rect.getArrayD();

    print_info();


// ----- Fill the matrices using MQ and its derivatives.      
    Mat Gstream(N, N), Gvort(N, N);
    Vec b(N), lstream(N), lvort(N);
    Mat sgrid(Nx, Ny);//, sout(Ngx, Ngy);
    Mat vgrid(Nx, Ny);//, vout(Ngx, Ngy);
    Vec uvel(N), vvel(N);

    streamLaplace( Gstream, x, y, RBF::MQ<prec_t,2>(c),
		   RBF::MQ_2DX<prec_t,2>(c), RBF::MQ_2DY<prec_t,2>(c) );

// ----- Preconditioner construction
    DiagonalPrec<Mat> preStream(Gstream);

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
    prec_t tol = 1.0e-8, tolerance = 1.0e-6, error = 0, local_error; 

    for(int t = 1; t <= max_iter; ++t) 
    {
	niter = gmres (Gstream, lstream, b, preStream, N-1 , tol);    	


	evalSolution(sgrid, x, y, lstream, c, beta, "strm_");

	local_error = exchangeInfo(sgrid, x, y, lstream, c, cart,
				   left, right, down, up);
	cart.comm.Allreduce(&local_error, &error, 1, MPI::DOUBLE, MPI_MAX);

	std::cout << "\n | Rank = " << rank << " | iter = " << t 
		  << " | S : G It : " << niter
		  << " | EL = " << local_error
		  << " | EG = " << error;

	evalVelocity(uvel, vvel, x, y, lstream, c);

	vortEquation(Gvort, x, y, 
		     RBF::MQ<prec_t,2>(c),
		     RBF::MQ_1DX<prec_t,2>(c),
		     RBF::MQ_1DY<prec_t,2>(c),
		     RBF::MQ_2DX<prec_t,2>(c), 
		     RBF::MQ_2DY<prec_t,2>(c),
		     uvel, vvel, Re);

// ----- Vorticity boundary condtions
	b = 0;
	if (neighbor[LEFT] == -1) {
// Edges 1: Dirichlet
	    bi = point[0];
	    ei = bi + side_y;
	    for(int i = bi, k = 2; i < ei; ++i, ++k) {	  	  	    
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
	    }
	}

	
	if (neighbor[RIGHT] == -1) {
// Edge 2: Dirichlet
	    bi = point[2];      
	    ei = bi + side_y;
	    for(int i = bi, k = 2; i < ei; ++i, ++k) {
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
	    }
	}

	if (neighbor[DOWN] == -1) {
// Edges 3: Dirichlet
	    bi = point[4];
	    ei = bi + side_x;
	    for(int i = bi, k = 2; i < ei; ++i, ++k) {	  	  
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
	    }
	}

	if (neighbor[UP] == -1) {
// Edges 4: Dirichlet (sliding wall)
	    bi = point[6];
	    ei = bi + side_x;	
	    for(int i = bi, k = 2; i < ei; ++i, ++k) {	  	  	  
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
	    }
	}
// Corners:
/* * 
// Average
	b(N-3) = 0.5 * ( b(point[0]) + b(point[4]) );
	b(N-2) = 0.5 * ( b(point[1]) + b(point[6]) );
	b(N-1) = 0.5 * ( b(point[2]) + b(point[5]) );
	b(N  ) = 0.5 * ( b(point[3]) + b(point[7]) );
/* */

	if (neighbor[DOWN] == -1 && neighbor[LEFT] == -1) 
	    b(N-3) = -3 * sgrid(2,2) / (dx * dx) 
		- ( 0.5 * vgrid(2,1) + 0.5 * vgrid(1,2) + 
		    0.25 * vgrid(2,2) );

	if (neighbor[DOWN] == -1 && neighbor[RIGHT] == -1) 
	    b(N-1) = -3 * sgrid(Nx-1,2) / (dx * dx) 
		- ( 0.5 * vgrid(Nx-1,1) + 0.5 * vgrid(Nx,2) + 
		    0.25 * vgrid(Nx-1,2) );	

	if (neighbor[UP] == -1 && neighbor[LEFT] == -1)
	    b(N-2) = -4.5 / dx - 3 * sgrid(2,Ny-1) / (dx*dx)
		- ( 0.5 * vgrid(2,Ny) + 
		    0.5 * vgrid(1,Ny-1) + 
		    0.25 * vgrid(2,Ny-1) );
	    
	if (neighbor[UP] == -1 && neighbor[RIGHT] == -1)
	    b(N  ) = -4.5 / dx - 3 * sgrid(Nx-1,Ny-1) / (dx*dx) 
		- ( 0.5 * vgrid(Nx-1,Ny) + 
		    0.5 * vgrid(Nx,Ny-1) + 
		    0.25 * vgrid(Nx-1,Ny-1) );

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

	if (neighbor[DOWN] == -1 && neighbor[LEFT] != -1) b(N-3) = vgrid(2,1);
	if (neighbor[DOWN]==-1 && neighbor[RIGHT]!=-1) b(N-1) = vgrid(Nx-1,1);
	if (neighbor[UP] == -1 && neighbor[LEFT] != -1) b(N-2) = vgrid(2,Ny);
	if (neighbor[UP] == -1 && neighbor[RIGHT] != -1) b(N) = vgrid(Nx-1,Ny);
	
	vgrid(1 , 1 ) = b(N-3);
	vgrid(1 , Ny) = b(N-2);
	vgrid(Nx, 1 ) = b(N-1);
	vgrid(Nx, Ny) = b(N  );

	DiagonalPrec<Mat> preVort(Gvort);
	niter = gmres (Gvort, lvort, b, preVort, N-1 , tol);
	evalSolution(vgrid, x, y, lvort, c, beta, "vort_");

	local_error = exchangeInfo(vgrid, x, y, lvort, c, cart,
				   left, right, down, up);
	cart.comm.Allreduce(&local_error, &error, 1, MPI::DOUBLE, MPI_MAX);

	std::cout << "\n | Rank = " << rank << " | iter = " << t 
		  << " | V : G It : " << niter
		  << " | EL = " << local_error
		  << " | EG = " << error;

	b = 0;
	id = 1;
	for(int j = 2; j < Ny; ++j)
	    for(int i = 2; i < Nx; ++i, ++id) 
		b(id) = -vgrid(i,j);

// ----- Artificial boundary conditions for streamfunction
	id = point[0];
	for(int i = 2; i <= Ny - 1; ++i, ++id) 
	    b(id) = sgrid(1,i);

	id = point[2];
	for(int i = 2; i <= Ny - 1; ++i, ++id) 
		b(id) = sgrid(Nx,i);   

	id = point[4];
	for(int i = 2; i <= Nx - 1; ++i, ++id) 
	    b(id) = sgrid(i, 1);     	    

	id = point[6];
	for(int i = 2; i <= Nx - 1; ++i, ++id) 
	    b(id) = sgrid(i, Ny);	    

	b(N-3) = sgrid(1,1);
	b(N-2) = sgrid(1,Ny);
	b(N-1) = sgrid(Nx,1);
	b(N  ) = sgrid(Nx,Ny);
    }


    std::cout << "\n +-----+  ";
    std::cout << "\n + Happy finish :-) ";
    std::cout << "\n +-----+ \n\n";    

    MPI::Finalize();
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
		  const Vec& lstream, prec_t c)
{
    filename = "u_" + istr + "_" + jstr + ".dat";
    std::ofstream u_file (filename.c_str());
    filename = "v_" + istr + "_" + jstr + ".dat" ;
    std::ofstream v_file (filename.c_str());

    for(int i = 1; i <= NI; ++i) {
	u(i) = RBF::eval(x(i), y(i), x, y, lstream,RBF::MQ_1DY<prec_t,2>(c),0);
	v(i) =-RBF::eval(x(i), y(i), x, y, lstream,RBF::MQ_1DX<prec_t,2>(c),0);
	u_file << x(i) << "\t" << y(i) << "\t" << u(i) << "\n";
	v_file << x(i) << "\t" << y(i) << "\t" << v(i) << "\n";
    }
/* */
    u_file.close();
    v_file.close();
/* */
}

//--------------------------------------------------------------------------
// Evaluate the numerical solution on a mesh of Nx by Ny
//

prec_t evalSolution(Mat& fgrid, const Vec& x, const Vec& y, const Vec& lambda,
		    prec_t c, prec_t beta, string name)
{
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


    filename = name + istr + "_" + jstr + ".dat";
    std::ofstream fgrid_file (filename.c_str());  
    for(int j = 1; j <= Ny; ++j) { 
	for(int i = 1; i <= Nx; ++i) {
	    fgrid_file << (i-1) * dx + sx << "\t" 
		       << (j-1) * dy + sy << "\t" 
		       << fgrid(i,j) << "\n";
	}
	fgrid_file << "\n";
    }
    fgrid_file.close();

    return error;
} 


//--------------------------------------------------------------------------
// Get information in the overlapping from neigbors processors.
//
prec_t exchangeInfo(Mat &f, const Vec& x, const Vec& y, 
		    const Vec& lambda, prec_t c, 
		    CartComm& cart, 
		    Vec& left, Vec& right, Vec& down, Vec& up)

{    

    filename = "a_" + istr + "_" + jstr + ".dat" ;
    std::ofstream a_file (filename.c_str());

    prec_t sol;
    int ei, bi;

// ----- Get the solution in the overlapping. An interpolation is done
//       using the lambda's calculated for the subdomain, to the overlapping
//       points stored in arrayL, arrayR, arrayU and arrayD from the class
//       OVRectangleKnots. The left, right, up and down vectors defined here
//       will contain the interpolation values at last.
    ei = left.length();
    Vec s_left(ei);
    for(int i = 1; i < ei; ++i) {
	s_left(i) = RBF::eval(left(ei), left(i), x, y, 
			      lambda, RBF::MQ<prec_t, 2>(c));
    }

    ei = right.length();
    Vec s_right(ei);
    for(int i = 1; i < ei; ++i) {
	s_right(i) = RBF::eval(right(ei), right(i), x, y, 
			       lambda,RBF::MQ<prec_t, 2>(c));
    }

    ei = down.length();
    Vec s_down(ei);
    for(int i = 1; i < ei; ++i) {
	s_down(i) = RBF::eval(down(i), down(ei), x, y, 
			      lambda,RBF::MQ<prec_t, 2>(c));
    }

    ei = up.length();
    Vec s_up(ei);
    for(int i = 1; i < ei; ++i) {
	s_up(i) = RBF::eval(up(i), up(ei), x, y, 
			    lambda,RBF::MQ<prec_t, 2>(c));
    }

// Vectors used to receive the information from neigbor subdomains.
    Vec r_left(Ny), r_right(Ny), r_up(Nx), r_down(Ny);

    MPI::Request req[8];
    MPI::Status  stat[8];

    MPI::Datatype Edge_x = MPI::DOUBLE.Create_contiguous(Nx);
    MPI::Datatype Edge_y = MPI::DOUBLE.Create_contiguous(Ny);
    Edge_x.Commit();
    Edge_y.Commit();

    const int *neighbor = cart.getNeighbors();

    req[0] = cart.comm.Isend(&s_down(1) , 1, Edge_x, neighbor[DOWN],  0);
    req[1] = cart.comm.Isend(&s_up(1)   , 1, Edge_x, neighbor[UP],    0);
    req[2] = cart.comm.Isend(&s_right(1), 1, Edge_y, neighbor[RIGHT], 0);
    req[3] = cart.comm.Isend(&s_left(1) , 1, Edge_y, neighbor[LEFT],  0);

    req[4] = cart.comm.Irecv(&r_up(1)   , 1, Edge_x, neighbor[UP],    0);
    req[5] = cart.comm.Irecv(&r_down(1) , 1, Edge_x, neighbor[DOWN],  0); 
    req[6] = cart.comm.Irecv(&r_left(1) , 1, Edge_y, neighbor[LEFT],  0);
    req[7] = cart.comm.Irecv(&r_right(1), 1, Edge_y, neighbor[RIGHT], 0);

    MPI::Request::Waitall(8, req, stat);    


    prec_t error = 0.0, avrg;
    if (neighbor[LEFT] != -1) {
	for(int i = 2; i <= Ny - 1; ++i) {
	    r_left(i) = beta * r_left(i) + (1 - beta) * f(1,i);
	    error += fabs( f(1,i) - r_left(i));
	    f(1,i) = r_left(i);
	    a_file << left(Ny+1) << "\t" << left(i) << "\t" << f(1,i) << "\n";
	}
	
// ----- Corner 1   //
	if (neighbor[DOWN] != -1) {  
	    avrg = (r_left(1) + r_down(1) ) * 0.5;
	    avrg = beta * avrg + (1 - beta) * f(1,1);
	    error += fabs(f(1,1) - avrg);    
	    f(1,1) = avrg;                   
	    a_file << left(Ny+1) << "\t" << left(1) << "\t"
		   << f(1,1) << "\n";	    
	} else {
	    r_left(1) = beta * r_left(1) + (1 - beta) * f(1,1);
	    error += fabs(f(1,1) - r_left(1));  
	    f(1,1) = r_left(1);
	    a_file << left(Ny+1) << "\t" << left(1) << "\t" << f(1,1) << "\n";
	} 

// ----- Corner 2   //
	if (neighbor[UP] != -1) {  
	    avrg = (r_left(Ny) + r_up(1) ) * 0.5;
	    avrg = beta * avrg + (1 - beta) * f(1,Ny);
	    error += fabs(f(1,Ny) - avrg);    
	    f(1,Ny) = avrg;                   
	    a_file << left(Ny+1) << "\t" << left(Ny) << "\t" << f(1,Ny) << "\n";	    
	} else {
	    r_left(Ny) = beta * r_left(Ny) + (1 - beta) * f(1,Ny);
	    error += fabs(f(1,Ny) - r_left(Ny));  
	    f(1,Ny) = r_left(Ny);
	    a_file << left(Ny+1) << "\t" << left(Ny) << "\t" << f(1,Ny) << "\n";	    
	} 	
    }
	

    if (neighbor[RIGHT] != -1) {
	for(int i = 2; i <= Ny - 1; ++i) {
	    r_right(i) = beta * r_right(i) + (1 - beta) * f(Nx,i);
	    error += fabs(f(Nx,i) - r_right(i));
	    f(Nx,i) = r_right(i);
	    a_file << right(Ny+1) << "\t" << right(i) << "\t"<< f(Nx,i) << "\n";
	}
// ----- Corner 3   //
	if (neighbor[DOWN] != -1) {  
	    avrg = (r_right(1) + r_down(Nx) ) * 0.5;
	    avrg = beta * avrg + (1 - beta) * f(Nx,1);
	    error += fabs(f(Nx,1) - avrg);    
	    f(Nx,1) = avrg;                   
	    a_file << right(Ny+1) << "\t" << right(1) << "\t"<< f(Nx,1) << "\n";	    
	} else {
	    r_right(1) = beta * r_right(1) + (1 - beta) * f(Nx,1);
	    error += fabs(f(Nx,1) - r_right(1));  
	    f(Nx,1) = r_right(1);
	    a_file << right(Ny+1) << "\t" << right(1) << "\t"<< f(Nx,1) << "\n";
	} 

// ----- Corner 4   //
	if (neighbor[UP] != -1) {  
	    avrg = (r_right(Ny) + r_up(Nx) ) * 0.5;
	    avrg = beta * avrg + (1 - beta) * f(Nx,Ny);
	    error += fabs(f(Nx,Ny) - avrg);    
	    f(Nx,Ny) = avrg;                   
	    a_file << right(Ny+1) << "\t" << right(Ny) << "\t"<< f(Nx,Ny) << "\n";	    
	} else {
	    r_right(Ny) = beta * r_right(Ny) + (1 - beta) * f(Nx,Ny);
	    error += fabs(f(Nx,Ny) - r_right(Ny));  
	    f(Nx,Ny) = r_right(Ny);
	    a_file << right(Ny+1) << "\t" << right(Ny) << "\t"<< f(Nx,Ny) << "\n";	    
	} 	
    }

    if (neighbor[DOWN] != -1) {
	for(int i = 2; i <= Nx - 1; ++i) {
	    r_down(i) = beta * r_down(i) + (1 - beta) * f(i,1);
	    error += fabs(f(i,1) - r_down(i));
	    f(i,1) = r_down(i);
	    a_file << down(i) << "\t" << down(Nx+1) << "\t"<< f(i,1) << "\n";
	}

// Corners 1 and 3 just in case there is neither left nor right neigbors
	if (neighbor[LEFT] == -1 && neighbor[RIGHT] == -1) {
	    r_down(1) = beta * r_down(1) + (1 - beta) * f(1,1);
	    error += fabs(f(1,1) - r_down(1));  
	    f(1,1) = r_down(1);
	    a_file << down(1) << "\t" << down(Nx+1) << "\t"<< f(1,1) << "\n";

	    r_down(Nx) = beta * r_down(Nx) + (1 - beta) * f(Nx,1);
	    error += fabs(f(Nx,1) - r_down(Nx));  
	    f(Nx,1) = r_down(Nx);
	    a_file << down(Nx) << "\t" << down(Nx+1) << "\t"<< f(Nx,1) << "\n";
	}
    }

    if (neighbor[UP] != -1) {
	for(int i = 2; i <= Nx - 1; ++i) {
	    r_up(i) = beta * r_up(i) + (1 - beta) * f(i,Ny);
	    error += fabs(f(i,Ny) - r_up(i));
	    f(i,Ny) = r_up(i);
	    a_file << up(i) << "\t" << up(Nx+1) << "\t"<< f(i,Ny) << "\n";
	}
// Corners 2 and 4 just in case there is neither left nor right neigbors
	if (neighbor[LEFT] == -1 && neighbor[RIGHT] == -1) {
	    r_up(1) = beta * r_up(1) + (1 - beta) * f(1,Ny);
	    error += fabs(f(1,Ny) - r_up(1));  
	    f(1,Ny) = r_up(1);
	    a_file << up(1) << "\t" << up(Nx+1) << "\t"<< f(1,Ny) << "\n";

	    r_up(Nx) = beta * r_up(Nx) + (1 - beta) * f(Nx,Ny);
	    error += fabs(f(Nx,Ny) - r_up(Nx));  
	    f(Nx,Ny) = r_up(Nx);
	    a_file << up(Nx) << "\t" << up(Nx+1) << "\t"<< f(Nx,Ny) << "\n";
	}

    }

    a_file.close();

    return error;
}

//--------------------------------------------------------------------------
//
bool read_data_and_Bcast(CartComm& cart) 
{
    int    dataBCast_Int[8];
    prec_t dataBCast_Flt[9];

    if (rank == 0) {
	cout << "\n\n"
	     << "                  MASTER PROCESS rank = " << rank << "\n"
	     << " +----------------------------------------------------+\n"
	     << " |       RADIAL BASIS FUNCTION FOR PDE SOLVING        |\n"
	     << " +----------------------------------------------------+\n"
	     << " | Author: Luis M. de la Cruz                         |\n"
	     << " | Date  : Tue Mar 18 12:12:41 GMT 2008               |\n"
	     << " +----------------------------------------------------+\n";
	
// ----- Reading data from "input" file	
	ifstream input_file ("input");
	input_file >> dataBCast_Flt[0]; // length in x-axis
	input_file >> dataBCast_Flt[1]; // length in y-axis
	input_file >> dataBCast_Int[0]; // Num. of pts in x-axis
	input_file >> dataBCast_Int[1]; // Num. of pts in y-axis
	input_file >> dataBCast_Flt[2]; // Size of overlapping in x-axis
	input_file >> dataBCast_Flt[3]; // Size of overlapping in y-axis
	input_file >> dataBCast_Int[2]; // 0 - Unif, 1 - Rand, 2 - Rand Unif
	input_file >> dataBCast_Flt[4]; // Randomness (0-1.0)
	input_file >> dataBCast_Int[3]; // Maximum num of iterations
	input_file >> dataBCast_Flt[5]; // Reynolds number
	input_file >> dataBCast_Int[4]; // Grid for evaluation: x-axis
	input_file >> dataBCast_Int[5]; // Grid for evaluation: y-axis
	input_file >> dataBCast_Int[6]; // Max knots for ACBF precond
	input_file >> dataBCast_Flt[6]; // Tolerance for GMRES
	input_file >> dataBCast_Flt[7]; // Tolerance for DDM
	input_file >> dataBCast_Int[7]; // Maximum num of DDM steps
	input_file >> dataBCast_Flt[8]; // Beta, under-relaxation factor
	input_file.close();
    }
// ----- Broadcast the info to all processors
    cart.comm.Bcast(dataBCast_Int, 8, MPI::INT, 0);    
    cart.comm.Bcast(dataBCast_Flt, 9, MPI::DOUBLE, 0); 

// ----- Using global variables.
    hx        = dataBCast_Flt[0];
    hy        = dataBCast_Flt[1];
    over_x    = dataBCast_Flt[2];
    over_y    = dataBCast_Flt[3];
    ep        = dataBCast_Flt[4];
    Re        = dataBCast_Flt[5];
    tol_gmres = dataBCast_Flt[6];
    tol_ddm   = dataBCast_Flt[7];
    beta      = dataBCast_Flt[8];
    Nx            = dataBCast_Int[0];
    Ny            = dataBCast_Int[1];
    rtype         = dataBCast_Int[2];	
    max_iter      = dataBCast_Int[3];
    Ngx           = dataBCast_Int[4];
    Ngy           = dataBCast_Int[5];
    max_knots     = dataBCast_Int[6];
    max_ddm       = dataBCast_Int[7];    
    return 0;
}

void print_info() {
    std::cout << "\n rank = " << rank 
	      << "  (I,J) = (" << I << "," << J << ")"
	      << "  hx = " << hx << " hy = " << hy 
	      << "  Nx = " << Nx << " Ny = " << Ny
	      << "  dx = " << dx << " dy = " << dy
	      << "\n";
}

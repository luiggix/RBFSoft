/*! 
 ***************************************************************************
 *
 *  \file vtkformat.cpp
 *  Reading and writting VTK format.
 ***************************************************************************
 *  \author  Luis M. de la Cruz Mon Mar 10 10:55:42 GMT 2008
 ***************************************************************************
 */

#include <ctime>    
#include <cstdlib>  
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

#include "Traits.hpp"
#include "Knots/BoxKnots.hpp"

//==========================================================================
//                            FUNCTION PROTOTYPES
//
void initialVelocity(Vec&, Vec&, const Vec&, const Vec&, 
		     prec_t, prec_t);

void writeToFile_VTK(Vec&, Vec&, std::string);

//==========================================================================
//                            GLOBAL DATA
//
int Nx, Ny, Nz, N, NI, NB;
prec_t hx, hy, hz;

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

    int max_iter, max_knots;
    prec_t ep, c, dt, A;

    std::ifstream input_file ("input");	
    input_file >> hx;   // length in x-axis
    input_file >> hy;   // length in y-axis
    input_file >> hz;   // length in z-axis
    input_file >> Nx; // Boundary points in x-axis
    input_file >> Ny; // Boundary points in y-axis
    input_file >> Nz; // Boundary points in y-axis
    input_file >> ep;   // Randomness of the interior points (0-0.5)
    input_file >> A; // amplitude;
    input_file.close();

 // ----- Point generation 

    std::cout << "\n +-----+  "
	      << "\n | Calculating knots ... ";  
   
    BoxKnots<prec_t> box(hx, Nx, hy, Ny, hz, Nz);
    N  = box.getTotalKnots();
    NI = box.getInteriorKnots();
    NB = box.getBoundaryKnots();
    box.print();
    time.tic();
    box.constructKnots();
    box.writeToFile("xy_knots.dat");

    Vec x = box.getKnots(X);
    Vec y = box.getKnots(Y);
    Vec z = box.getKnots(Z);

    std::cout << "\n | Knots generation elapsed time = " << time.toc()
	      << "\n +-----+  ";

    c = 1.0 / sqrt (1.0 * N);

// ----- Initial condition
    Vec uvel(N), vvel(N), wvel(N); 
    prec_t rolls = 1.0;
    initialVelocity(uvel, vvel, x, y, A, rolls);
    writeToFile_VTK(uvel, vvel, "velocity.vtk");

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
    std::ofstream uvel_file ("uvel.dat");
    std::ofstream vvel_file ("vvel.dat");
    int cnt = 1;
    for(int k = 2; k < Nz; ++k) {
	for(int j = 2; j < Ny; ++j) {
	    for(int i = 2; i < Nx; ++i, ++cnt) {
		u(cnt) = -A * cos(PI * y(cnt)) * sin(PI * rolls * x(cnt) /hx);
		v(cnt) = A2 * sin(PI * y(cnt)) * cos(PI * rolls * x(cnt) /hx);
		uvel_file << x(cnt) << "\t" << y(cnt) << "\t" << u(cnt) <<"\n";
		vvel_file << x(cnt) << "\t" << y(cnt) << "\t" << v(cnt) <<"\n";
	    }	    
	    uvel_file << "\n";
	    vvel_file << "\n";	
	}
	uvel_file << "\n";
	vvel_file << "\n";
    }
    uvel_file.close();
    vvel_file.close();
}

void writeToFile_VTK(Vec& u, Vec& v, std::string name)
{
    std::ofstream file (name.c_str());
    
    prec_t dx = hx / (Nx - 1);
    prec_t dy = hy / (Ny - 1);
    prec_t dz = hz / (Nz - 1);
        
    file << "# vtk DataFile Version 2.0 \n"
	 << "Velocity in internal nodes\n" 
	 << "ASCII \n"
	 << "DATASET RECTILINEAR_GRID \n" 
	 << "DIMENSIONS " << Nx-2 << " " << Ny-2 << " " << Nz-2 << "\n";

    file << "X_COORDINATES " << Nx-2 << " float \n";
    for(int i = 1; i < Nx - 1; ++i) {
	file << i * dx << " ";
	if ( !(i % 6) ) file << "\n";
    }

    file << "Y_COORDINATES " << Ny-2 << " float \n"; 
    for(int i = 1; i < Ny - 1; ++i) {
	file << i * dy << " ";
	if ( !(i % 6) ) file << "\n";
    }

    file << "Z_COORDINATES " << Nz-2 << " float \n"; 
    for(int i = 1; i < Nz - 1; ++i) {
	file << i * dz << " ";
	if ( !(i % 6) ) file << "\n";
    }
	
    file << "POINT_DATA " << NI << "\n"
	 << "VECTORS velocity float \n";
//	 << "LOOKUP_TABLE default \n";
    for(int i = 1; i <= NI; ++i) {
	file << u(i) << " " << v(i) << " " << 0 << "\n";
    }
        
    file.close();

}


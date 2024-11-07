
/*! 
 ***************************************************************************
 *
 *  \file GNUplotToVTK.cpp
 *  Reading GNUplot format and writting VTK format in 2D.
 ***************************************************************************
 *  \author  Luis M. de la Cruz [ Mon May 12 13:48:03 BST 2008 ]
 ***************************************************************************
 */

#include <iostream>
#include <fstream>
#include <string>

int main( int argc, char * argv[])
{
    int Nvtk_x, Nvtk_y, Nvtk_z, Ntotal;
    double Lvtk_x, Lvtk_y, Lvtk_z;
    double dx, dy, dz;
    double x_dummy, y_dummy, data;
    std::string in_name, out_name, variable;

    std::cout << "\n +---- GNUPLOT =-> VTK -----+";
    std::cout << "\n | Nx = "; std::cin >> Nvtk_x;
    std::cout << "\n | Ny = "; std::cin >> Nvtk_y;
    std::cout << "\n | Nz = "; std::cin >> Nvtk_z;
    std::cout << "\n | Lx = "; std::cin >> Lvtk_x;
    std::cout << "\n | Ly = "; std::cin >> Lvtk_y;
    std::cout << "\n | Lz = "; std::cin >> Lvtk_z;
    std::cout << "\n | Input file name = "; std::cin >> in_name;    
    std::cout << "\n | Ouput file name = "; std::cin >> out_name;
    std::cout << "\n | Variable name = "; std::cin >> variable;
    std::cout << "\n +-----+ ";

    dx = Lvtk_x / (Nvtk_x - 1);
    dy = Lvtk_y / (Nvtk_y - 1);
//    dz = Lvtk_z / (Nvtk_z - 1);
    dz = 1.0;
    Ntotal = Nvtk_x * Nvtk_y * Nvtk_z;

    std::ifstream in_file (in_name.c_str());        
    std::ofstream out_file (out_name.c_str());      
  
    std::cout << "\n | Transforming file ... ";

    out_file << "# vtk DataFile Version 2.0 \n"
	 << "Lid_driven streamfuntion-vorticity \n"
	 << "ASCII \n"
	 << "DATASET RECTILINEAR_GRID \n" 
	 << "DIMENSIONS " 
	 << Nvtk_x << " " << Nvtk_y << " " << Nvtk_z;

    out_file << "\nX_COORDINATES " << Nvtk_x << " float \n";
    for(int i = 0; i < Nvtk_x; ++i) {
	out_file << i * dx << " ";
	if ( !(i % 6) ) out_file << "\n";
    }

    out_file << "\nY_COORDINATES " << Nvtk_y << " float \n"; 
    for(int i = 0; i < Nvtk_y; ++i) {
	out_file << i * dy << " ";
	if ( !(i % 6) ) out_file << "\n";
    }

    out_file << "\nZ_COORDINATES " << Nvtk_z << " float \n"; 
    for(int i = 0; i < Nvtk_z; ++i) {
	out_file << i * dz << " ";
	if ( !(i % 6) ) out_file << "\n";
    }
	
    out_file << "\nPOINT_DATA " << Ntotal
	 << "\nSCALARS " << variable << " float"
	 << "\nLOOKUP_TABLE default \n";

    double mat[41][41];
    for(int k = 1; k <= Nvtk_z ; ++k) {
	for(int j = 1; j <= Nvtk_y ; ++j) {
	    for(int i = 1; i <= Nvtk_x; ++i) {
//		in_file >> x_dummy >> y_dummy >> data;
		in_file >> x_dummy >> y_dummy >> mat[i-1][j-1];
//		out_file << data << " ";
//		if ( !(i % 6) ) out_file << "\n";
	    }
	}
    }

    for(int i = 0; i < Nvtk_x; ++i) {
	for(int j = 0; j < Nvtk_y ; ++j) {
	    out_file << mat[i][j] << " ";
	    if ( !(i % 6) ) out_file << "\n";
        }
    }
    out_file.close();
    in_file.close();

    std::cout << "\n | File " << out_name << " generated !!!";
    std::cout << "\n +-----+ ";

    return 0;
}

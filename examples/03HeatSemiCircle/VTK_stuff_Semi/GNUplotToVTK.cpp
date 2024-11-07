
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
    double x_dummy, y_dummy;
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
    int ncells = (Nvtk_x - 1) * (Nvtk_y - 1);


    std::cout << "\n NX = " << Nvtk_x
	      << " NY = " << Nvtk_y
	      << " NZ = " << Nvtk_z
	      << " Points = " << Ntotal
	      << " Cells = " << ncells << "\n"
;
    float *data = new float[Ntotal];

    std::ifstream in_file (in_name.c_str());        
    std::ofstream out_file (out_name.c_str());      
  
    std::cout << "\n | Transforming file ... ";

    out_file << "# vtk DataFile Version 2.0 \n"
	     << "Unstructured example \n"
	     << "ASCII \n"
	     << "DATASET UNSTRUCTURED_GRID \n" 
	     << "POINTS " << Ntotal << " float \n";

    int id = 0;
    for(int k = 1; k <= Nvtk_z ; ++k) {
	for(int j = 1; j <= Nvtk_y ; ++j) {
	    for(int i = 1; i <= Nvtk_x; ++i, ++id) {
		in_file >> x_dummy >> y_dummy >> data[id];		
		out_file << x_dummy << " " << y_dummy << " " << 0.0 << "\n";
	    }
	}
    }

    int NI = Nvtk_x;

    out_file << "\nCELLS " << ncells << " " << ncells * 5 <<" \n";
    for(int j = 0; j < Nvtk_y -1 ; ++j) 
	for(int i = 0; i < Nvtk_x -1; ++i) 
	    out_file << 4 << " " 
		     << j * NI + i << " "
		     << j * NI + i + 1 << " "
		     << (j+1) * NI + i + 1 << " "
		     << (j+1) * NI + i << " "
		     << "\n";

    out_file << "\nCELL_TYPES " << ncells << "\n";
    
    for(int i=1; i<=ncells; ++i)
	out_file << 9 << " ";

    out_file << "\nPOINT_DATA " << Ntotal
	     << "\nSCALARS " << variable << " float 1"
	     << "\nLOOKUP_TABLE default \n";
    for(int i = 0; i < Ntotal; ++i) {
	out_file << data[i] << " ";
	if ( !(i % 6) ) out_file << "\n";
    }
        
    out_file.close();
    in_file.close();

    std::cout << "\n | File " << out_name << " generated !!!";
    std::cout << "\n +-----+ ";

    return 0;
}

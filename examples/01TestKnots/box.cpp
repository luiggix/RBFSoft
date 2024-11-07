/*! 
 ***************************************************************************
 *  \file    box.cpp
 *  Testing the BoxKnots class.
 * In this example some features of the class BoxKnots are tested.
 * Particularly the function Knots::findNeighbors() is used to find the
 * neighborhood of a target point.
 * \par Input   
 * The \c inputBox file contains the input data required for this program:
 * \c hx length in x-axis;
 * \c hy length in y-axis;
 * \c hz length in z-axis;
 * \c Nx number of points in x-axis;
 * \c Ny number of points in y-axis;
 * \c Nz number of points in z-axis;
 * \par Output
 * \c xyzBox.dat coordinates of random points;
 * \c tarBox.dat the target point;
 * \c neiBox.dat list of neighbors.
 * \par Post-procesing
 * You can plot the results using the next command in gnuplot:
 *  \verbatim
    % gnuplot> splot "neiBox.dat" w lp, "xyzBox.dat" w p, "tarBox.dat" w p \endverbatim  

 * <p>
 * \image html  boxtest.png "Knots, target and neighborhood." width=5cm 
 * \image latex boxtest.eps "Knots, target and neighborhood." width=5cm
 ***************************************************************************
 *  \author  Luis M. de la Cruz [ Mon Apr 28 10:31:36 BST 2008 ]
 ***************************************************************************
 */

#include <ctime>    
#include <iostream>
#include <fstream>
#include <vector>
#include "Traits.hpp"
#include "Knots/BoxKnots.hpp"
#include "GNUplot.hpp"

int main( int argc, char * argv[])
{
    timer time;

    std::cout << "\n\n"
	      << " +----------------------------------------+\n"
	      << " |       TUNA::RBF FOR PDE SOLVING        |\n"
	      << " +----------------------------------------+\n"
	      << "\n";

    int Nx, Ny, Nz;
    prec_t hx, hy, hz;

    std::ifstream input_file("inputBox");
    input_file >> hx
	       >> hy
	       >> hz
	       >> Nx   
	       >> Ny
	       >> Nz;
    input_file.close();

    BoxKnots<prec_t> box(hx, Nx, hy, Ny, hz, Nz);
    int N  = box.getTotalKnots();
    int NI = box.getInteriorKnots();
    int NB = box.getBoundaryKnots();
    box.print();

    Vec x(N), y(N), z(N);      
    std::cout << "\n +-----+  ";
    std::cout << "\n | Calculating knots ... ";
    time.tic();
    box.constructKnots();
    std::cout << "\n | Knots generation elapsed time = " << time.toc() ;
    x = box.getKnots(X);
    y = box.getKnots(Y);
    z = box.getKnots(Z);
    box.writeToFile("xyzBox.dat");

    std::cout << "\n +-----+  "; 
    std::cout << "\n | Constructing KDTree ... ";
    time.tic();
    box.calcKDTree();
    std::cout << "\n | KDTree elapsed time = " << time.toc() ;
    std::cout << "\n +-----+  "; 

    int target = N/2;
    prec_t range = hx * 0.2;
    std::vector<knot3D> ne;
    std::cout << "\n +-----+  "; 
    std::cout << "\n | Finding the neighbors points to ( "
	      << x(target) << ",\t" << y(target) << ",\t" << z(target) << " )";
    time.tic();
    ne = box.findNeighbors(x(target), y(target), z(target), range);
    std::cout << "\n | KDTree elapsed time = " << time.toc() ;
    std::cout << "\n +-----+  "; 

    std::ofstream file_target ("tarBox.dat");
    file_target << x(target) << "\t" << y(target) << "\t" << z(target) << "\n";
    file_target.close();

    std::ofstream file_out ("neiBox.dat");
    std::vector<knot3D>::const_iterator ne_i;
    ne_i = ne.begin();
    for (; ne_i != ne.end(); ++ne_i) 
	file_out << *ne_i;    
    file_out.close();

    std::cout << "\n\n >---> Files Generated :  xyzBox.dat \n"
              << "                          neiBox.dat \n" 
	      << "                          tarBox.dat \n";
    
#ifdef WITH_GNUPLOT
    int pausa;
    GNUplot plotter;
    plotter("splot \"neiBox.dat\" w lp,\"xyzBox.dat\" w p,\"tarBox.dat\" w p");
    std::cout << "\n\n >---> Press any key and then <enter> to finish " ;
    std::cin >> pausa;
#endif

    return 0;

}


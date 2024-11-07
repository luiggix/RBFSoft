/*! 
 ***************************************************************************
 *  \file    line.cpp
 *  Testing the LineKnots class.
 * In this example some features of the class LineKnots are tested.
 * Particularly the function Knots::findNeighbors() is used to find the
 * neighborhood of a target point. 
 * \par Input   
 *  The \c inputLine contains the input data required for this program:
 * \c hx length in x-axis;
 * \c Nx number of points in x-axis;
 * \par Output
 * \c xyzLine.dat coordinates of random points;
 * \c tarLine.dat the target point;
 * \c neiLine.dat list of neighbors of the target.
 * \par Post-procesing
 * You can plot the results using the next command in gnuplot:
 *  \verbatim
    % gnuplot> p "neiLine.dat" w lp, "xyzLine.dat" w p, "tarLine.dat" w p \endverbatim  
 * <p>
 ***************************************************************************
 *  \author  Luis M. de la Cruz [ Thu Sep  6 14:35:41 BST 2007 ]
 ***************************************************************************
 */

#include <ctime>    
#include <iostream>
#include <fstream>
#include <vector>
#include "Traits.hpp"
#include "Knots/LineKnots.hpp"
#include "GNUplot.hpp"

int main( int argc, char * argv[])
{
    timer time;

    std::cout << "\n\n"
	      << " +----------------------------------------+\n"
	      << " |       TUNA::RBF FOR PDE SOLVING        |\n"
	      << " +----------------------------------------+\n"
	      << "\n";
    int Nx;
    prec_t hx;

    std::ifstream input_file("inputLine");
    input_file >> hx >> Nx;   
    input_file.close();

    LineKnots<prec_t> line(hx, Nx);
    int N  = line.getTotalKnots();
    int NI = line.getInteriorKnots();
    int NB = line.getBoundaryKnots();
    line.print();

    Vec x(N);      
    std::cout << "\n +-----+  ";
    std::cout << "\n | Calculating knots ... ";
    time.tic();
//    line.readFromFile("xyzRec.dat"); // You can also read points from a file
    line.constructKnots();
    std::cout << "\n | Knots generation elapsed time = " << time.toc() ;
    x = line.getKnots(X);
    line.writeToFile("xyzLine.dat");

    std::cout << "\n +-----+  "; 
    std::cout << "\n | Constructing KDTree ... ";
    time.tic();
    line.calcKDTree();
    std::cout << "\n | KDTree elapsed time = " << time.toc() ;
    std::cout << "\n +-----+  "; 

    int target = N/4;
    prec_t range = hx * 0.1;
    std::vector<knot1D> ne;
    std::cout << "\n +-----+  "; 
    std::cout << "\n | Finding the neighbors points to ( "
	      << x(target) << " ) ";
    time.tic();
    ne = line.findNeighbors(x(target), range);
    std::cout << "\n | KDTree elapsed time = " << time.toc() ;
    std::cout << "\n +-----+  "; 

    std::ofstream file_target ("tarLine.dat");
    file_target << x(target) << "\t" << 0.0 << "\n";
    file_target.close();

    std::ofstream file_out ("neiLine.dat");
    std::vector<knot1D>::const_iterator ne_i;
    ne_i = ne.begin();
    double temp;
    for (; ne_i != ne.end(); ++ne_i)
      file_out << ne_i->coord[0] << "\t" << 0.0 << "\n";
    file_out.close();

    std::cout << "\n\n >---> Files Generated :  xyzLine.dat \n"
              << "                          neiLine.dat \n" 
	      << "                          tarLine.dat \n";

#ifdef WITH_GNUPLOT
    int pausa;
    GNUplot plotter;
    //    plotter("set xlabel \"Point Id\"");
    plotter("set xlabel \"x-coordinate\"");
    plotter("p \"neiLine.dat\" w lp,\"xyzLine.dat\" w p,\"tarLine.dat\" w p");
    std::cout << "\n\n >---> Press any key and then <enter> to finish " ;
    std::cin >> pausa;
#endif


    return 0;

}


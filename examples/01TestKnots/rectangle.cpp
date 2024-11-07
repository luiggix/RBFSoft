/*! 
 ***************************************************************************
 *  \file    rectangle.cpp
 *  Testing the RectangleKnots class.
 * In this example some features of the class RectangleKnots are tested.
 * Particularly the function Knots::findNeighbors() is used to find the
 * neighborhood of a target point.
 * \par Input   
 * The \c inputRectangle file contains the input data required for this program:
 * \c hx length in x-axis;
 * \c hy length in y-axis;
 * \c Nx number of points in x-axis;
 * \c Ny number of points in y-axis;
 * \c rtype knots distribution, (0 unif, 1 rand, 2 rand unif);
 * \c ep degree of randomness.
 * \c layer 0 no ghost points, 1 one layer of ghost points
 * \par Output
 * \c xyzRect.dat coordinates of random points;
 * \c tarRect.dat the target point;
 * \c neiRect.dat list of neighbors of the target.
 * \par Post-procesing
 * You can plot the results using the next command in gnuplot:
 *  \verbatim
    % gnuplot> p "neiRect.dat" w lp, "xyzRect.dat" w p, "tarRect.dat" w p \endverbatim     
 * <p>
 * \image html  rectest.png "Knots, target and neighborhood." width=5cm 
 * \image latex rectest.eps "Knots, target and neighborhood." width=5cm
 ***************************************************************************
 *  \author  Luis M. de la Cruz [ Thu Sep  6 14:35:41 BST 2007 ]
 ***************************************************************************
 */

#include <ctime>    
#include <iostream>
#include <fstream>
#include <vector>
#include "Traits.hpp"
#include "Knots/RectangleKnots.hpp"
#include "GNUplot.hpp"

int main( int argc, char * argv[])
{
    timer time;

    std::cout << "\n\n"
	      << " +----------------------------------------+\n"
	      << " |       TUNA::RBF FOR PDE SOLVING        |\n"
	      << " +----------------------------------------+\n"
	      << "\n";

    int Nx, Ny, rtype, layer;
    prec_t hx, hy, ep;

    std::ifstream input_file("inputRectangle");
    input_file >> hx
	       >> hy   
	       >> Nx   
	       >> Ny   
	       >> rtype
	       >> ep
	       >> layer;
    input_file.close();
    random_t RT = static_cast<random_t>(rtype);

    RectangleKnots<prec_t> rect(hx, Nx, hy, Ny, RT, layer);

//    RectangleKnots<prec_t> rect(hx, Nx, hy, Ny, 2, 2, 2, 2, -1, -1, RT, 1);
    int N  = rect.getTotalKnots();
    int NI = rect.getInteriorKnots();
    int NB = rect.getBoundaryKnots();
    rect.print();

    Vec x(N), y(N);      
    std::cout << "\n +-----+  ";
    std::cout << "\n | Calculating knots ... ";
    time.tic();
//    rect.readFromFile("xyzRec.dat"); // You can also read points from a file
    rect.setRandomness(ep);
    rect.constructKnots();
    std::cout << "\n | Knots generation elapsed time = " << time.toc() ;
    x = rect.getKnots(X);
    y = rect.getKnots(Y);
    rect.writeToFile("xyzRect.dat");

    std::cout << "\n +-----+  "; 
    std::cout << "\n | Constructing KDTree ... ";
    time.tic();
    rect.calcKDTree();
    std::cout << "\n | KDTree elapsed time = " << time.toc() ;
    std::cout << "\n +-----+  "; 

    int target = N/2;
    prec_t range = hx * 0.1;
    std::vector<knot2D> ne;
    std::cout << "\n +-----+  "; 
    std::cout << "\n | Finding the neighbors points to ( "
	      << x(target) << ",\t" << y(target) << " ) ";
    time.tic();
    ne = rect.findNeighbors(x(target), y(target), range);
    std::cout << "\n | KDTree elapsed time = " << time.toc() ;
    std::cout << "\n +-----+  "; 

    std::ofstream file_target ("tarRect.dat");
    file_target << x(target) << "\t" << y(target) << "\n";
    file_target.close();

    std::ofstream file_out ("neiRect.dat");
    std::vector<knot2D>::const_iterator ne_i;
    ne_i = ne.begin();
    for (; ne_i != ne.end(); ++ne_i) 
	file_out << *ne_i;    
    file_out.close();

    std::cout << "\n\n >---> Files Generated :  xyzRect.dat \n"
              << "                          neiRect.dat \n" 
	      << "                          tarRect.dat \n";

#ifdef WITH_GNUPLOT
    int pausa;
    GNUplot plotter;
    plotter("p \"neiRect.dat\" w lp,\"xyzRect.dat\" w p,\"tarRect.dat\" w p");
    std::cout << "\n\n >---> Press any key and then <enter> to finish " ;
    std::cin >> pausa;
#endif

    return 0;

}


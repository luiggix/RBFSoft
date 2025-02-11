/*! 
 ***************************************************************************
 *  \file    semicircle.cpp
 *  Testing the SemiCircleKnots class.
 * In this example some features of the class SemiCircleKnots are tested.
 * Particularly the function Knots::findNeighbors() is used to find the
 * neighborhood of a target point.
 * \par Input
 * The \c inputSemiCircle file contains the input data required for this program:
 * \c r1 and \c r2 Radii for defining the semicircle;
 * \c t1 and \c t2 Angles for defining the semicircle;
 * \c Nr number of points in R;
 * \c Nt number of points in theta;
 * \c rtype knots distribution (0 unif, 1 rand, 2 rand unif);
 * \c ep degree of randomness.
 * \par Output
 * \c xyzSemi.dat coordinates of random points;
 * \c tarSemi.dat the target point;
 * \c neiSemi.dat list of neighbors of the target.
 * \par Post-procesing
 * You can plot the results using the next command in gnuplot:
 *  \verbatim
    % gnuplot> p "neiSemi.dat" w lp, "xyzSemi.dat" w p, "tarSemi.dat" w p \endverbatim 
 * <p>
 * \image html  semtest.png "Knots, target and neighborhood." width=5cm 
 * \image latex semtest.eps "Knots, target and neighborhood." width=5cm
 ***************************************************************************
 *  \author  Luis M. de la Cruz [ Thu Sep  6 14:35:41 BST 2007 ]
 ***************************************************************************
 */

#include <iostream>
#include <fstream>
#include <vector>
#include "Traits.hpp"
#include "Knots/SemiCircleKnots.hpp"
#include "GNUplot.hpp"

int main( int argc, char * argv[])
{
    timer time;

    std::cout << "\n\n"
	      << " +----------------------------------------+\n"
	      << " |       TUNA::RBF FOR PDE SOLVING        |\n"
	      << " +----------------------------------------+\n"
	      << "\n";

    int Nt, Nr, rtype;
    prec_t r1, r2, t1, t2;  
	
    std::ifstream input_file("inputSemiCircle");
    input_file >> r1
	       >> r2
	       >> t1   	
	       >> t2    	
	       >> Nr
	       >> Nt
	       >> rtype;
    input_file.close();

    random_t RT = static_cast<random_t>(rtype);
   
    SemiCircleKnots<prec_t> semi(r1, r2, Nr, t1, t2, Nt, RT);
    int N  = semi.getTotalKnots();
    int NI = semi.getInteriorKnots();
    int NB = semi.getBoundaryKnots();
    semi.print();

    Vec x(N), y(N);      
    std::cout << "\n +-----+  ";
    std::cout << "\n | Calculating knots ... ";
    time.tic();
//    semi.readFromFile("xyzSemi.dat"); // You can also read points from a file
    semi.constructKnots();
    std::cout << "\n | Knots generation elapsed time = " << time.toc() ;
    x = semi.getKnots(X);
    y = semi.getKnots(Y);
    semi.writeToFile("xyzSemi.dat");

    std::cout << "\n +-----+  "; 
    std::cout << "\n | Constructing KDTree ... ";
    time.tic();
    semi.calcKDTree();
    std::cout << "\n | KDTree elapsed time = " << time.toc() ;
    std::cout << "\n +-----+  "; 

    int target = N/2 + Nr/2;
    prec_t range = (r2 - r1) * 0.1;
    std::vector<knot2D> ne;
    std::cout << "\n +-----+  "; 
    std::cout << "\n | Finding the neighbors points to ( "
	      << x(target) << ",\t" << y(target) << " ) ";
    time.tic();
    ne = semi.findNeighbors(x(target), y(target), range);
    std::cout << "\n | KDTree elapsed time = " << time.toc() ;
    std::cout << "\n +-----+  "; 

    std::ofstream file_target ("tarSemi.dat");
    file_target << x(target) << "\t" << y(target) << "\n";
    file_target.close();

    std::ofstream file_out ("neiSemi.dat");
    std::vector<knot2D>::const_iterator ne_i;
    ne_i = ne.begin();
    for (; ne_i != ne.end(); ++ne_i) 
	file_out << *ne_i;
    file_out.close();

    std::cout << "\n\n >---> Files Generated :  xyzSemi.dat \n"
              << "                          neiSemi.dat \n" 
	      << "                          tarSemi.dat \n";

#ifdef WITH_GNUPLOT
    int pausa;
    GNUplot plotter;
    plotter("p \"neiSemi.dat\" w lp,\"xyzSemi.dat\" w p,\"tarSemi.dat\" w p");
    std::cout << "\n\n >---> Press any key and then <enter> to finish " ;
    std::cin >> pausa;
#endif

    return 0;

}


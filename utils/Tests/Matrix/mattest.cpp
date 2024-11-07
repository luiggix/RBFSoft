#include <ctime>    
#include <cstdlib>  
#include <cmath>

#include <iostream>
#include <fstream>
#include <string>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
using namespace boost::numeric::ublas;


typedef double prec_t; 
typedef matrix<prec_t,row_major> MatrixRow;
typedef matrix<prec_t,column_major> MatrixCol;


//==========================================================================
//                            MAIN FUNCTION
//
int main( int argc, char * argv[])
{
  if (argc < 2){ 
      std::cout << "\n Use: " << argv[0] << " N \n\n"
		<< " N = Matrix dimension \n\n";
      exit(1);
    }

// Matrices and vectors declaration.
  int N = atoi(argv[1]);  

  clock_t start,finish;
  double time1, time2, time3, time4;

  MatrixRow AR(N, N); 
  MatrixCol AC(N, N);

  start = clock();  
  for(int j = 0; j < N; ++j) {
      for(int i = 0; i < N; ++i) {
	  AR(i,j) = i * j + sin(i*j);
      }
  }
  finish = clock();
  time1 = (double(finish)-double(start))/CLOCKS_PER_SEC;

  start = clock();  
  for(int i = 0; i < N; ++i) {
      for(int j = 0; j < N; ++j) {
	  AR(i,j) = i * j + sin(i*j);
      }
  }
  finish = clock();
  time2 = (double(finish)-double(start))/CLOCKS_PER_SEC;


  start = clock();  
  for(int j = 0; j < N; ++j) {
      for(int i = 0; i < N; ++i) {
	  AC(i,j) = i * j + sin(i*j);
      }
  }
  finish = clock();
  time3 = (double(finish)-double(start))/CLOCKS_PER_SEC;

  start = clock();  
  for(int i = 0; i < N; ++i) {
      for(int j = 0; j < N; ++j) {
	  AC(i,j) = i * j + sin(i*j);
      }
  }
  finish = clock();
  time4 = (double(finish)-double(start))/CLOCKS_PER_SEC;

  
  std::cout << "\n ----- TIMES REPORT ----- \n\n"
	    << "\n ROW MAJOR MATRIX: "
	    << " \t (J (I) ) : " << time1 
	    << " \t (I (J) ) : " << time2 
	    << "\n COLUMN MAJOR MATRIX: "
	    << " \t (J (I) ) : " << time3 
	    << " \t (I (J) ) : " << time4 
	    << "\n";
  
  return 0;
}


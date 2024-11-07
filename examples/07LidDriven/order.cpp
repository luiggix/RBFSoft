

#include <iostream>
#include <cmath>

int main()
{
/* */
// Eurtuk : V max
//    double h1 = 1.0 / 41, u1 = 0.18649;
    double h1 = 1.0 / 31, u1 = 0.180868;
    double h2 = 1.0 / 21, u2 = 0.172595;
    double h3 = 1.0 / 11, u3 = 0.143695;
/* *
// Jensen : V max
//    double h1 = 1.0 / 41, u1 = 0.191721;
    double h1 = 1.0 / 31, u1 = 0.18837;
    double h2 = 1.0 / 21, u2 = 0.184521;
    double h3 = 1.0 / 11, u3 = 0.151909;

/* *
// Eurtuk : V min
//    double h1 = 1.0 / 41, u1 = -0.263961;
    double h1 = 1.0 / 31, u1 = -0.258305;
    double h2 = 1.0 / 21, u2 = -0.250718;
    double h3 = 1.0 / 11, u3 = -0.198528;
/* *
// Eurtuk : U min
//    double h1 = 1.0 / 41, u1 = -0.227871;
    double h2 = 1.0 / 31, u2 = -0.228078;
    double h3 = 1.0 / 21, u3 = -0.227327;
    double h3 = 1.0 / 11, u3 = -0.202053;

/* *
    double h1 = 1.0 / 64, u1 = 0.17896;
    double h2 = 1.0 / 44, u2 = 0.17829;
    double h3 = 1.0 / 32, u3 = 0.17709;
/* */
    double error = 1.0, alpha, a_old = 0.01;
    double tolerance = 1.0e-6;
    int iter = 1, max_iter = 1000;

    double A = (u3 - u1) / (u2 - u1); 
    double h2_h1 = h2 / h1;
    double div = log ( h3 / h1 );
    double dif, sign;

    while ( iter < max_iter) {
	alpha = log ( A * pow(h2_h1, a_old) - A + 1) / div;       
	dif = alpha - a_old;	

	if ( iter == 1 ) sign = dif;
	if ( dif < 0 && sign > 0) break;
	if ( dif > 0 && sign < 0) break;

	error = fabs( dif );

	std::cout << "\n | Iter = " << iter 
		  << " | Alpha old = " << a_old
		  << " | Alpha = " << alpha
		  << " | Error = " << error << "\n";
	a_old = a_old + 0.01;
	iter++;       	
    }

    double C = (u2 - u1) / ( pow(h2, alpha) - pow(h1, alpha) );

    double Uex = u1 - C * pow(h1, alpha);

    std::cout << "\n Uex = " << Uex; 

    return 0;
}

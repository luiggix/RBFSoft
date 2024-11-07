
#include <fstream>
#include <iostream>
#include <string>
using namespace std;

int main() {

    ifstream ifile ("solStream.dat");



    int NI = 60;
    int NJ = 15;

    double x[NI][NJ], y[NI][NJ], s[NI][NJ], u;

    for(int j = 0; j < NJ; ++j) 
	for(int i = 0; i < NI; ++i)
	    ifile >> x[i][j] >> y[i][j] >> s[i][j];
       
    int mon;
    string name;
    cout << "\n >---> Name = "; cin >> name;
    cout << "\n >---> Mon = "; cin >> mon;
    ofstream ofile (name.c_str());

    double shift = mon * 30 / 59;
    ofile << shift << " " << 0 << "\n";
    for(int j = 1; j < NJ - 1; ++j) {
	u = (s[mon][j+1] - s[mon][j-1]) / (y[mon][j+1] - y[mon][j-1]);  
	ofile << u + shift << "\t" << y[mon][j] << "\n";
    }
    ofile << shift << " " << 1 << "\n";

    return 0;
}

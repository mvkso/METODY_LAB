#include <iostream>
#include <cmath>

using namespace std;
const double e = 1e-8;
const int i_max = 100;

double max(double a, double b, double c){
	return (a < b) ? ((b < c) ? c : b) : ((a < c) ? c : a);
}

double f1(double x, double y, double z){
	return x * x + y * y + z * z - 2.0;
}

double f2(double x, double y){
	return x * x + y * y - 1.0;
}

double f3(double x, double y){
	return x * x - y;
}

int main()
{
    cout.setf(ios::scientific);
	cout.precision(16);
    return 0;
}

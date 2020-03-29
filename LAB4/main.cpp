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
double maximum(double x,double y, double z){
    return (x<y)?((y<z)?z:y):((x<z)?z:x);
}

void Newton(double x, double y, double z)
{
    double est, rez, det_a,w_odwrotny;
    double *values = new double[3];

    //jakobian przeksztalcenia oznaczony J^-1
    double **macierz_jakob = new double*[3];

    for (int i = 0; i < 3; i++){
		macierz_jakob[i] = new double[3];
	}
	//poprawka jaka odejmiemy od ->xn zeby uzyskaæ ->xn+1
	double *delta = new double[3];

	values[0] = f1(x, y, z);
	values[1] = f2(x, y);
	values[2] = f3(x, y);

	cout << "i" << " " << "x" << "\t\t\t   " << "y" << "\t\t\t   " << "z" << "\t\t\t   " << "est"
		 << "\t\t\t  " << "rez" << endl;

	det_a=(-4 * x * z) - (8 * x * y * z);
	w_odwrotny = 1 / det_a;

	for (int i = 0; i < i_max; i++){
		//zerujemy macierz delta[]
		fill(delta, delta + 3, 0);
		if (!x || !z){
			cout << "Dzielenie przez zero" << endl;
			break;
		}

		macierz_jakob[0][0] = 0;
		macierz_jakob[0][1] = (-2 * z) * w_odwrotny;
		macierz_jakob[0][2] = (-4 * y * z) * w_odwrotny;
		macierz_jakob[1][0] = 0;
		macierz_jakob[1][1] = (-4 * x * z) * w_odwrotny;
		macierz_jakob[1][2] = (4 * x * z) * w_odwrotny;
		macierz_jakob[2][0] = (-2 * x - 4 * x * y) * w_odwrotny;
		macierz_jakob[2][1] = (2 * x + 4 * x * y) * w_odwrotny;
		macierz_jakob[2][2] = 0;

		/*
		 * Wyznaczamy deltê poprzez mno¿enie macierzy jakobianu i wartoœci, mno¿ymy macierze
		 */
		for (int j = 0; j < 3; j++){
			for (int k = 0; k < 3; k++){
				delta[j] = delta[j] + (macierz_jakob[j][k]*values[k]);
			}
		}
		//odejmujemy roznice
		x-=delta[0];
		y-=delta[1];
		z-=delta[2];
		//przypisujemy nowe wartosci
		values[0]=f1(x, y, z);
		values[1]=f2(x, y);
		values[2]=f3(x, y);

		//obliczamy estymator i reziduum dla normy inf
		est = maximum(fabs(delta[0]), fabs(delta[1]), fabs(delta[2]));
		rez = maximum(fabs(values[0]), fabs(values[1]), fabs(values[2]));

		cout << i << " " << x << " " << y << " " << z << " " << est << " " << rez << endl;

		//badamy warunki zakonczenia iteracji
		if (fabs(rez) < e) break;
		if (fabs(est) < e) break;
	}

}

int main()
{
    cout.setf(ios::scientific);
	cout.precision(16);
	Newton(1.0, 2.0, 3.0);
    return 0;
}

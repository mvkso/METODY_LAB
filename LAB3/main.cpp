#include <iostream>
#include <cmath>

using namespace std;


const char *f1 = "Funkcja SINUSA: sin(x/4)^2-x";
const char *f2 = "Funkcja TANGESA: tg(2x)-x-1";

//zdefiniowanie wlasnego typu aby przekazywac dowolna funkcje zwracajaca double o argumencie double
typedef double (*function)(double);

//wartosc funkcji numer1
double f1(double x){return sin(x / 4.0) * sin(x / 4.0) - x;} //sin((x)/4)^2 -x
//wartosc funckji numer2
double f2(double x){return tan(2.0 * x) - x - 1.0;} //tg(2x)-x-1

//pochodne
double fp1(double x){return (1.0 / 4.0) * (sin(x / 2.0)) - 1;} // (1/4)*(sin(x/2)-1)
double fp2(double x){return -1.0 + 2.0 / (cos(2.0 * x) * cos(2.0 * x));}// -1 + 2/(cos(2x)^2)


//FI(X) - do metody picarda wymagane jest przekszta³cenie f(x)=0 do fi(x)=x gdzie fi(x)=f(x)+x
double f1_wpic(double x){return sin(x / 4.0) * sin(x / 4.0);} //picard dla sin(x/4)^2
double f2_wpic(double x){return tan(2.0 * x) - 1.0;} //picard dla tg(2x)-1

// FI '(X) - dodatkowe pochodne dla metody picarda aby sprawdzic zbieznosci
double f1p_wpic(double x){return 1.0 / 4.0 * sin(x / 2.0);}
double f2p_wpic(double x){return 2.0 / (cos(2.0 * x) * cos(2.0 * x));}



double Picard(function f, function f_iter, function fp, double x, int n_max, double TOLX, double TOLF){
	if (fabs(fp(x)) >= 1){
		cout << "Pochodna FI'(x) >= 1 - \rozbie¿noœæ, Metoda Picarda nie przyblizy pierwiastka";
		return 0;

	}
	double est = 0,rez = 0;
	double new_approx = x;

	cout << "Wartosc x=" << x << endl;
	//iteracja
	for (int i = 0; i < n_max; i++){
		cout << "iteracja i=" << i;

		//przypisanie nowej wartosci x-nastepne = wartosc funkcji, bo y=x
		new_approx = f_iter(new_approx);
		cout << ", wartosc przyblizona=" << new_approx;

		//zmiana estymatora i reziduum
		est = fabs(new_approx - x);
		x = new_approx;
		rez = fabs(f(x));

		cout << ", residuum= " << rez << ", estymator= " << est << endl;

//		warunki zakoczenia iteracji
		if ((rez <= TOLF) or (est <= TOLX))
			break;
	}
	return new_approx;
}


int main()
{
    cout << "Hello world!" << endl;
    return 0;
}

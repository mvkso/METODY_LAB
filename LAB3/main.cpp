#include <iostream>
#include <cmath>

using namespace std;



//zdefiniowanie wlasnego typu aby przekazywac dowolna funkcje zwracajaca double o argumencie double
typedef double (*function)(double);

//wartosc funkcji numer1
double f1(double x){return sin(x / 4.0) * sin(x / 4.0) - x;} //sin((x)/4)^2 -x
//wartosc funckji numer2
double f2(double x){return tan(2.0 * x) - x - 1.0;} //tg(2x)-x-1

//pochodne
double fp1(double x){return (1.0 / 4.0) * (sin(x / 2.0)) - 1;} // (1/4)*(sin(x/2)-1)
double fp2(double x){return -1.0 + 2.0 / (cos(2.0 * x) * cos(2.0 * x));}// -1 + 2/(cos(2x)^2)


//FI(X) - do metody picarda przekszta³cenie f(x)=0 na fi(x)=x - fi(x)=f(x)+x
double fi_1(double x){return sin(x / 4.0) * sin(x / 4.0);} //picard dla sin(x/4)^2
double fi_2(double x){return tan(2.0 * x) - 1.0;} //picard dla tg(2x)-1

// FI '(X)
double fi_p1(double x){return 1.0 / 4.0 * sin(x / 2.0);}
double fi_p2(double x){return 2.0 / (cos(2.0 * x) * cos(2.0 * x));}



double Picard(function f, function f_iter, function fp, double x, int nmax, double TOLX, double TOLF){
	if (fabs(fp(x)) >= 1)
        {
		cout << "FI'(x) >= 1 -> rozbieznosc";
		return 0;
		}else cout<<"FI'(X)<1 -> ZBIEZNOSC"<<endl;
	double est = 0,rez = 0;
	double new_approx = x;

	cout << "x=" << x << endl;
	//iteracja
	for (int i = 0; i < nmax; i++){
		cout << "i=" << i;

		//przypisanie nowej wartosci x-nastepne = wartosc funkcji, bo y=x
		new_approx = f_iter(new_approx);
		cout <<  "\t"<<", przyblizenie=" << new_approx;

		//zmiana estymatora i reziduum
		est = fabs(new_approx - x);
		x = new_approx;
		rez = fabs(f(x));

		cout <<  "\t"<<", rez= " << rez <<  "\t"<<", est= " << est << endl;

//		warunki zakoczenia iteracji
		if ((rez <= TOLF) or (est <= TOLX))
			break;
	}
	return new_approx;
}

double Bisekcja(function fun, double a, double b, int nmax, double TOLX, double TOLF){
	double xn = 0; //srodek przedzialu
    double est=0,rez=0;
	for (int i = 0; i < nmax; i++){
		cout << "ITERACJA i=" << i<<endl;
//		sprawdzenie czy konce przedzialow sa roznych znakow
		if ((fun(a) > 0 and fun(b) > 0) or (fun(a) < 0 and fun(b) <0)){
			cout << "Przedzial nie spelnia wymagan" << endl;
			return EXIT_FAILURE;
		}
		else{
			cout << "a= " << a <<", b= "<<b;
//			srodek przedzialu
			xn = (a + b) / 2.0;
//			estymatora
			est = fabs((b - a) / 2.0);
			rez=fabs(fun(xn));

			cout << ", xn= " << xn <<  "\t"<<",est= " <<est << "\t"<< ", rez= " <<rez << endl;
			if ((fun(a) < 0 and fun(xn) > 0) or (fun(xn) < 0 and fun(a) > 0)){
				//istnieje zero w przedziale [a;xn]
				b = xn;
			}
			else{
				a = xn;
			}

		}
//		warunki zakoczenia iteracji
		if ((fabs(fun(xn)) <= TOLF) or (fabs((b - a) / 2) <= TOLX)) break;
	}
	if (!fun(xn)) cout << "x0= " << xn << "(" << fun(xn) << ")" << endl;
}

double Newton(function f, function p, double x, int nmax, double TOLX, double TOLF)//metoda stycznych
{
	double x0 = x,x1;
	double est = 0,rez=0;
	for (int i = 0; i < nmax; i++){
//		ze wzoru
		x1 = x0 - (f(x0) / p(x0));
		est = fabs(x0 - x1);
		x0 = x1;
		rez=(fabs(f(x0)));
            cout << "i=" << i<< ", x1= " << x1 << ", est= " << est << ", rez=" << f(x0) << endl;;
		if (rez<= TOLF or est <= TOLX) break;
	}
}

double Sieczne(function f, double x0, double x1, int nmax, double TOLX, double TOLF){
	double x2 = 0;
	double est = 0,rez=0;
	cout << "x0=" << x0 << ", x1=" << x1 << endl;
	for (int i = 0; i < nmax; i++){
//		ze wzoru
		x2 = x1 - f(x1) / ((f(x1) - f(x0)) / (x1 - x0));
		est = fabs(x2 - x1);
		rez=fabs(f(x2));
            cout << "i=" << i<< "\t,x2=" << x2 << "\t,est= " << fabs(est) << "\t,rez=" << rez << endl;
		x0 = x1;
		x1 = x2;
            if (rez <= TOLF or est <= TOLX) break;
	}
	return x2;
}



int main()
{
    cout<<"KAMIL MAKSYMOWICZ GL05"<<endl;
    cout.setf(ios::scientific);
	cout.precision(16);
	double e=1e-8;
	const int n_max=100;
	cout<<"METODA PICARDA: "<<endl<<endl;
	cout<<"Funkcja: sin(x/4)^2-x"<<endl;
	Picard(f1, fi_1, fi_p1, 0.5, n_max, e, e);
	cout<<endl<<"Funkcja: tg(2x)-x-1"<<endl;
	Picard(f2, fi_2, fi_p2, 0.5, n_max, e, e);
	cout << endl << endl << "Metoda Bisekcji:" <<endl<< endl;
	cout<<"Funkcja: sin(x/4)^2-x"<<endl;
	Bisekcja(f1, -0.5, 0.9, n_max, e, e);
	cout<<endl<<"Funkcja: tg(2x)-x-1"<<endl;
	Bisekcja(f2, 0.4, 0.5, n_max, e, e);
	cout << endl << endl << "Metoda Newtona"<<endl << endl;
	cout<<"Funkcja: sin(x/4)^2-x"<<endl;
	Newton(f1, fp1, -0.3, n_max, e, e);
	cout<<endl<<"Funkcja: tg(2x)-x-1"<<endl;
	Newton(f2, fp2, 0.6, n_max, e, e);
	cout << endl << endl << "Metoda Siecznych"<<endl  << endl;
	cout<<"Funkcja: sin(x/4)^2-x"<<endl;
	Sieczne(f1, -0.5, 0.6, n_max, e, e);
	cout<<endl<<"Funkcja: tg(2x)-x-1"<<endl;
	Sieczne(f2, 0.4, 0.5, n_max, e, e);




    return 0;
}

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

using namespace std;

double funkcja(double x){
	 return ( exp(2.0-2.0*x) - 4.0*exp(4.0-2.0*x) + 4.0*exp(2.0*x) -exp(2.0+2.0*x) - x + x*exp(4.0) ) / ( 4.0 - 4.0*exp(4.0) );
}

const double p = 1.0, q = 0.0, r = -4.0;
const double Alfa = 0.0, Beta = 1.0, Gamma = -1.0;
const double Fi = 0.0, Psi = 1.0, Theta = -0.5;
const double xPoczatkowe = 0.0, xKoncowe = 1.0;

double maksymalny_blad(double *blad, const int rozmiar);
double trzypunktowa_dyskretyzacja_konwencjonalna(double h, const int max);
double dyskretyzacja_numerowa_jednorodna(double h, const int max);
void rozwiaz(double *d, double *u, double *l, double *b, double*x, int max);
void thomas(double *d, double *u, double *l,int max);



int main(){
	int max;
	double h;
	double bladK, bladN;
	fstream BLEDY;

	BLEDY.open("Wynik.txt", ios::out);
	if(!BLEDY.good()) {
		perror("Nie udalo sie otworzyc pliku");
	}
    BLEDY.precision(8);

	for(max = 10; max < 30000; max += 40) {

		h = (xKoncowe - xPoczatkowe) / (max - 1);
		bladK = log10(trzypunktowa_dyskretyzacja_konwencjonalna(h, max));
		bladN = log10(dyskretyzacja_numerowa_jednorodna(h, max));
		BLEDY << log10(h) << " " << bladK << " " << bladN << endl;
	}
	BLEDY.close();
	return 0;


}



void rozwiaz(double *d, double *u, double *l, double *b, double*x, int max){
	for(int i = 1; i < max; i++) {
		b[i] = b[i] - ((l[i - 1] / d[i - 1]) * b[i - 1]);
	}

	x[max - 1] = b[max - 1] / d[max - 1];
	for(int i = max - 2; i >= 0; i--) {
		x[i] = (b[i] - u[i] * x[i + 1]) / d[i];
	}
}
void thomas(double *d, double *u, double *l,int max){
	for(int i = 1; i < max; i++) {
		d[i] = d[i] - ((l[i - 1] / d[i - 1]) * u[i - 1]);
	}
}


double maksymalny_blad(double *blad, int rozmiar){
	double max = abs(blad[0]);
	for(int i = 0; i < rozmiar; i++) {
		if(fabs(blad[i] > max))
			max = abs(blad[i]);
	}
	return max;
}




double trzypunktowa_dyskretyzacja_konwencjonalna(double h, int max){
    double *d, *u, *l, *b, *x, *blad, s, *tmp;
    double x1=xPoczatkowe, x2=xPoczatkowe,zwroc;
    fstream PLIK;
    int i;
    PLIK.precision(8);
    tmp=new double[max];
    for (i=0 ; i<max ; i++)
        tmp[i]=0;

    d = new double[max];
	u = new double[max];
	l = new double[max];
	b = new double[max];
	x = new double[max];
	blad = new double[max];

    s = xPoczatkowe;

	d[0] = Beta-Alfa/h;
	u[0] = Alfa/h;
	//l[0] = 0;
	b[0] = -Gamma;

	for (i = 1 ; i < max-1 ; i++)
	{
		l[i-1] = p / (h*h) - q/(2.0*h);
		d[i] = -2.0*p / (h*h) + r;
		u[i] = p / (h*h) + q / (2.0*h);
		b[i] = (s+h*i);
	}

	l[max-2] = -Fi / h;
	d[max-1] = Fi/h + Psi;
	b[max-1] = -Theta;

	//algorytm Thomasa
	//thomas(l,d,u,b,x,N);
	//thomasMacierz(d,u,l,N,tmp);
	//thomasWektor(tmp,b,N);
	//rozwiaz(d,u,l,b,x,N);
	thomas(d,u,l,max);
	rozwiaz(d,u,l,b,x,max);

	for ( int i = 0; i < max; i++ )
	{
		blad[i] = fabs( x[i] - funkcja( x1 ) );
		x1 += h;
	}


	if ( max == 50 )
	{
		PLIK.open( "konwencjonalna.txt", fstream::out );

		cout << endl << "Trzypunktowa dyskretyzacja konwencjonalna" << endl;
		cout.width(4);
		cout<<"i";
		cout.width(15);
		cout<<"punkt";
		cout.width(15);
		cout<<"x[n]";
		cout.width(15);
		cout<<"U(x)"<<endl;

		for ( int j = 0; j < max; j++ )
		{
			blad[j] = fabs( x[j] - funkcja( x2 ) );
			PLIK << x2 << " " << x[j] << " " << funkcja( x2 ) << " " << endl;
			cout.width(4);
			cout << j;
			cout.width(15);
			cout << x2;
			cout.width(15);
			cout << x[j];
			cout.width(15);
			cout << funkcja(x2)<< endl;

			x2 += h;
		}
		PLIK.close( );
	}

	delete[] d;
	delete[] u;
	delete[] l;
	delete[] b;
	delete[] x;
	delete[] tmp;

    zwroc=maksymalny_blad( blad, max );
    delete[] blad;
	return zwroc;
}




double dyskretyzacja_numerowa_jednorodna(double h, const int max){
    double *d, *u, *l, *b, *x, *blad, s, *tmp;
    double x1=xPoczatkowe, x2=xPoczatkowe;
    double zwroc;
    fstream plik;
    int i;
    plik.precision(8);
    tmp=new double[max];
    for (i=0 ; i<max ; i++)
        tmp[i]=0;

    d = new double[max];
	u = new double[max];
	l = new double[max];
	b = new double[max];
	x = new double[max];
	blad = new double[max];

    s = xPoczatkowe;

	d[0] = Beta-Alfa/h;
	u[0] = Alfa/h;
	//l[0] = 0;
	b[0] = -Gamma;

	for (i = 1 ; i < max-1 ; i++)
	{
		l[i-1] = p / (h*h) + r/(12.0);
		d[i] = (-2.0 * p) / (h * h) + r * (10.0 / 12.0);
		u[i] = p / (h * h) + r / 12.0;
		b[i] = (s+i*h-h) / 12.0 + (10.0/12.0) *  (s+i*h) + (s+i*h+h) / 12.0;
	}

	l[max - 2] = -Fi / h;
	d[max - 1] = -Fi / h + Psi;
	b[max - 1] = -Theta;

	//algorytm Thomasa
	//thomas(l,d,u,b,x,N);
	//thomasMacierz(d,u,l,N,tmp);
	//thomasWektor(tmp,b,N);
	//rozwiaz(d,u,l,b,x,N);
	thomas(d,u,l,max);
	rozwiaz(d,u,l,b,x,max);

	for ( int i = 0; i < max; i++ )
	{
		blad[i] = fabs( x[i] - funkcja( x1 ) );
		x1 += h;
	}


	if ( max == 50 )
	{
		plik.open( "Numerow.txt", fstream::out );

		cout << endl << "Dyskretyzacja Numerowa" << endl;
		cout.width(4);
		cout<<"i";
		cout.width(15);
		cout<<"punkt";
		cout.width(15);
		cout<<"x[n]";
		cout.width(15);
		cout<<"U(x)"<<endl;

		for ( int j = 0; j < max; j++ )
		{
			blad[j] = fabs( x[j] - funkcja( x2 ) );
			plik << x2 << " " << x[j] << " " << funkcja( x2 ) << " " << endl;
			cout.width(4);
			cout << j;
			cout.width(15);
			cout << x2;
			cout.width(15);
			cout << x[j];
			cout.width(15);
			cout << funkcja(x2)<< endl;

			x2 += h;
		}
		plik.close( );
	}

	delete[] d;
	delete[] u;
	delete[] l;
	delete[] b;
	delete[] x;
	delete[] tmp;

    zwroc=maksymalny_blad( blad, max );
    delete[] blad;
	return zwroc;
}

// projekt_MO.cpp : Defines the entry point for the console application.
//

//#include "stdafx.h"

#include <iostream>
#include <cmath>
#include <fstream>
#include "calerf.h"
#include "macierz.h"
#include "wektor.h"

void rozw_analityczne (double **analit,int ile_x, int ile_t, double h, double dt);
void lassonen(double **lasson,int ile_x, int ile_t, double h, double dt);
void blad (double **blad, double **analit, double **lasson,int ile_x, int ile_t, double h, double dt);
void thomas( double *a, double *d, double *c, double *b, double *u, int n );
using namespace std;

double D = 1.;   //wspolczynnik ciepla
double t_max = 2.;  // max t (na osi y)
double x_min = -6.*sqrt(D*t_max);
double x_max = 6.*sqrt(D*t_max);  // max na osi x
//double lambda = 1; //lambda dla metody bezposredniej

class metoda
{

    public:
      // double **laasonen, **analit, **blad, *blad_max;

	   macierz laasonen;
	   macierz analit;
	   macierz blad;
	   wektor blad_max;
       metoda ();
       int ile_x;
	   double lambda_A;
	   double lambda_L;
       double h;
       double dt;
       int ile_t;
       void rozw_analityczne ();
       void rozw_laasonen ();
     //  void thomas ( double *up, double *diag, double *low, double *b, double *wyn);
       void f_blad ();
       void zapis_do_pliku_m (char* nazwa,  macierz);
      // void zapis_do_pliku_w (char* nazwa,  wektor);
	 friend wektor algThomasa(wektor &b,wektor &a,wektor &d,wektor &c);

};

metoda::metoda () //konstruktor wylicza przedzialy i alokuje pamiec na macierze rozwiazan
{
	lambda_A = 0.4;
	lambda_L = 1.0;
    ile_x = 40;
    h = (x_max -x_min)/ ile_x;

	dt =  h*h*lambda_L / D;
    ile_t = static_cast<int> ( (t_max/dt) + 1 );

	laasonen= macierz(ile_t,ile_x);

  //  laasonen = new double*[ile_t];
  //  for (int i=0;i<ile_t;i++)
  //  laasonen[i] = new double[ile_x];

	dt =  h*h*lambda_L / D;
    ile_t = static_cast<int> ( (t_max/dt) + 1 );
  //  analit = new double*[ile_t];
	analit=macierz(ile_t,ile_x);
   // for (int i=0;i<ile_t;i++)
   // analit[i] = new double[ile_x];

	blad=macierz(ile_t,ile_x);
   // blad = new double*[ile_t];
   // for (int i=0;i<ile_t;i++)
   // blad[i] = new double[ile_x];

	blad_max=wektor(ile_t);
 //   blad_max = new double[ile_t];
}

void metoda::rozw_analityczne ()
{

//	ile_x = 40;
	//h = x_max / ile_x;
	dt =  h*h*lambda_L / D;
    ile_t = static_cast<int> ( (t_max/dt) + 1 );

	double t = 0,x = x_min;
    for (int j = 0; j < ile_t; j++)
    {
		for (int i = 0; i < ile_x ; i++)
		{
			analit(j,i) = erfc (x / (2.0*sqrt(D*t)))*0.5; // rozwiazanie analityczne
			x = x+h;
		}
    x =  x_min;
    t = t+dt;
    }
    zapis_do_pliku_m ("rozw_analityczne_l.txt", analit);
}

void metoda::rozw_laasonen ()
{

	//ile_x = 40;
	//h = x_max / ile_x;

	dt =  h*h*lambda_L / D;
    ile_t = static_cast<int> ( (t_max/dt) + 1 );

	double x = x_min,t=0;
  //   double *up = new double[ile_x];
	 wektor up(ile_x);
  //   double *b = new double[ile_x];
	 wektor b(ile_x);
   //  double *low = new double[ile_x];
	 wektor low(ile_x);
    // double *diag = new double[ile_x];
	  wektor diag(ile_x);
    // double *wyn = new double[ile_x];
	  wektor wyn(ile_x);

	  /*

	  warunek poczatkowe: U(x,0) = 1 if(x<0) else if(x>=0) 0.0

	  warunek brzegowy :  U(x_min,t)=1  U(x_max,t)=0;

	  */
	  double xx = x_min;
	  for( int i = 0; i < ile_x; i++ ){ //warunek poczatkowy
          if(xx<0)
			 laasonen(0,i) = 1.0;
		  else
			 laasonen(0,i) = 0.0;
		  xx=xx+h;

	  }


	  for( int i = 0; i < ile_t; i++ ) // warunki brzegowe
         laasonen(i,0) = 1.0;

	  for( int i = 0; i < ile_t; i++ )
          laasonen(i,ile_x-1) = 0.0;

          /*Wypelnianie macierzy */
	  for( int k = 1; k < ile_t; k++ )
     {
          up[0] = 0;
          diag[0] = 1;
          b[0] = laasonen(k-1,0);

          for( int i = 1; i < ile_x-1; i++ )
          {
               up[i] = lambda_L;
               diag[i] = -( 1 + (2*lambda_L) );
               low[i] = lambda_L;
               b[i] = -laasonen(k-1,i);
          }

          up[ile_x-1] = 0;
          diag[ile_x-1] = 1;
          b[ile_x-1] = laasonen(k-1,ile_x-1);
         // thomas( up, diag, low, b, wyn);
		  wyn = algThomasa(b,up,diag,low);
          for( int i = 1; i < ile_x-1; i++ )
          laasonen(k,i) = wyn[i];  // wyliczenie wektora rozwiazan Thomasem i zapisanie do kolejnego wiersza
     }
     zapis_do_pliku_m ("rozw_metoda_laasonen.txt", laasonen);
}
 /*
void metoda::thomas( double *up, double *diag, double *low, double *b, double *wyn)
{
	int i;
	low[0] = low[0] / diag[0];
	b[0] = b[0] / diag[0];
	for( i = 1; i < ile_x; i++ )
    {
		double id = (diag[i] - low[i-1] * up[i]);
		low[i] = low[i] / id;
		b[i] = ( b[i] - b[i-1] * up[i] ) / id;
	}
	// wyliczanie wyniku
	wyn[ile_x - 1] = b[ile_x - 1];
	for( i = ile_x - 2; i >= 0; i-- )
		wyn[i] = b[i] - low[i] * wyn[i+1];
}
*/
void metoda::f_blad ()
{
	dt =  h*h*lambda_L / D;
    ile_t = static_cast<int> ( (t_max/dt) + 1 );
     double   x =  x_min,t = 0;
     for (int j = 0; j < ile_t; j++)
     {
         for (int i = 0; i < ile_x ; i++)
         {
             blad(j,i)=fabs(analit(j,i) - laasonen(j,i));
             if (blad(j,i)>blad_max[j]) blad_max[j] = blad(j,i); //blad bezwzgledny
             x=x+h;
         }
     t=t+dt;
        x =  x_min;
     }
    // zapis_do_pliku_m ("blad_bezwzgledny_dla_laasonen.txt", blad);
     //zapis_do_pliku_w ("blad_max_dla_laasonen.txt", blad_max);
}

void metoda::zapis_do_pliku_m (char* nazwa,  macierz m)
{

     ofstream file( nazwa );
     double x = x_min, t = 0;
     file<<nazwa<<", h = "<<h<<", dt = "<<dt<<endl<<endl;
     for (int j = 0; j < ile_t; j++)
     {

		 x = x_min;
         file<<"Dla poziomu czasowego = "<<t<<endl<<endl;
         for (int i = 0; i < ile_x ; i++)
         {
             file<<"wartosc przestrzenna: "<<x<<",\t U(x,t) = "<<m(j,i)<<endl;
             x=x+h;
         }


     file<<"\n"<<endl;
     t = t+dt;

     }
     file.close();
}
 /*
void metoda::zapis_do_pliku_w (char* nazwa, double* wektor)
{
     ofstream file( nazwa );
     double t = 0;
     file<<nazwa<<", dt =  "<<dt<<endl<<endl;
     for (int j = 0; j < ile_t; j++)
     {
         file<<"Dla poziomu czasowego = "<<t<<"\t blad_max = "<<wektor[j]<<endl;
         t = t+dt;
     }
     file.close();
}
 */
int main ()
{
    metoda nowa;
    nowa.rozw_analityczne ();
    nowa.rozw_laasonen ();
    nowa.f_blad ();
  /*  cout<<"Program rozwiazuje rownanie rozniczkowe czastkowe opisujace trasport ciepla"<<endl;
    cout<<"METODA POSREDNIA LASSONEN ORAZ ALGORYTMEM THOMASA"<<endl<<endl;
    cout<<"Wzor rownania: dU(x,t) / dt = D d^2U(x,t) / dx^2"<<endl<<endl;
    cout<<"Wspolrzedna przestrzenna okreslona na zbiorze: [0; +nieskonczonosci)"<<endl;
    cout<<"Wspolrzedna czasowa okreslona na zbiorze: [0; t_max]"<<endl;
    cout<<"Warunek poczatkowy: U(x,0) = 0"<<endl;
    cout<<"Warunki brzegowe U(0,t) = 1, U(+nieskonczonosci,t) = 0"<<endl<<endl;;
    cout<<"t_max = 2,  D = 1, lambda = 1 \nprzedzial nieskonczonosc zastapiony przedzialem >= od 6*sqrt(D*t_max)"<<endl;
    cout<<"Przedzial ten wynosi: "<<x_max<<endl<<endl;
    cout<<"Ilosc podzialow na osi x opisujacej przestrzej wybrano: "<<nowa.ile_x<<endl;
    cout<<"Wyliczony krok na osi przestrzennej h = "<<nowa.h<<endl;
    cout<<"Krok na osi czasowej wyliczono ze wzoru: lambda = D*dt/h^2 i wyniosl dt = "<<nowa.dt<<endl;
    cout<<"Ilosc podzialow na osi t opisujacej czas: "<<nowa.ile_t<<endl;
    cout<<"\n\n\nPoszczegolne wyniki dla kolejnych poziomow czasowych zapisano w plikach oraz blad bezwzgledy";
    cout<<" w porownaniu do rozwiazania analitycznego"<<endl;
	*/
    system ("pause");
}

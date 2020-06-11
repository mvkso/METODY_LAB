
#include<iostream>
#include<cmath>
#include<fstream>
#include "calerf.h"

using namespace std;
//====================================== Sta³e ========================================
    const double D=1.0;
    const double tmax=2.0;
    const double a=6*sqrt(D*tmax);
    const double lambda=1;
    const double h=0.1; //krok przestrzenny
    const double dt=(lambda*(h*h)/D);
    const int N=static_cast<const int>((tmax/dt)+2); //ilosc poziomow czasowych / wierszy
    const int M = static_cast<const int>((a*2/h)+1); //ilosc poziomow przestrz. / kolumn


//============================ Program korzysta z funkcji =============================
    double *nowy_wektor( int n );
    double **nowa_macierz( int n, int m );
    void usun_macierz( double **matrix, int n );
    void usun_wektor( double *x );
    void wypisz_wektor( double *x, int k );
    void wypisz_macierz( double **x, int r, int c );
    void fill_dt(double *x, double delta, int k);
    void fill_h(double *x, double delta, int k);
    void zapisz_plik( char* filename, double **matrix, int r, int c);
    double analityczne(double **ANAL, int N, int M, double *dt_levels, double *h_levels, double dt, double h);
    void Thomas( double *l, double *d, double *u, double *b, double *x, int n );
    void dekompozycja_LU(double **U, double *b, double *wyn);
    void CrankNicolson( double **U, int N, int M, double *dt_levels, double *h_levels, double dt, double h, double lambda_l);
    void LaasonenThomas( double **U, int N, int M, double *dt_levels, double *h_levels, double dt, double h, double lambda_l);
    void LaasonenLU( double **U, int N, int M, double *dt_levels, double *h_levels, double dt, double h, double lambda_l);
    void blad(double **BLAD, double **ANAL, double **U, int N, int M);



//=====================================================================================
int main()
{
    cout.setf( ios::fixed, ios::floatfield );
    cout.precision(1);

    double *dt_levels = nowy_wektor( N ); // wektor poziomow czasowych
    double *h_levels = nowy_wektor( M ); // wektor poziomow przestrzennych
    fill_dt(dt_levels, dt, N);
    fill_h(h_levels, h, M);

    char *filename = "hdt_CN.txt"; // zapis poziomow czasowych i przestrzennych do pliku:

    ofstream file( filename );
    if( !file.is_open() )
    {
        cout << "Nie moge zapisac: " << filename <<" !!!"<< endl;
        exit(1);
    }
    file.setf(ios::fixed, ios::floatfield);
    file.precision(4);
    file << endl << "      POZIOMY dt i h:             " << endl << endl;
    file << "dt\\h: \t";
    for( int i=0; i<M; i++)
         file  << h_levels[i]<< "\t";
    file << endl;
    for(int i=0; i<N; i++)
         file <<  dt_levels[i]<<endl;
    file.close();
    cout << "Zapisano: " << filename  << endl;

    //--------------------------------------------------------------------------------------

    double **ANAL = nowa_macierz( N, M );
    analityczne( ANAL, N, M, dt_levels, h_levels, dt, h );


    double **U = nowa_macierz( N, M );
    LaasonenThomas(U,N,M,dt_levels, h_levels,dt,h,lambda);

    double **A = nowa_macierz(N,M);
    LaasonenLU(A,N,M,dt_levels,h_levels,dt,h,lambda);

    double **BLAD = nowa_macierz(N,M);
    blad(BLAD, ANAL, U, N, M);

    usun_macierz( ANAL, N );
    usun_macierz( U, N );
    usun_macierz(A,N)
    usun_wektor(dt_levels);
    usun_wektor(h_levels);
    system("PAUSE");
   return 0;
}

double *nowy_wektor( int n )
{
       double *x = new double[n];
       return x;
}

double **nowa_macierz( int n, int m ) // n - liczba wierszy, m - liczba kolumn
{
       double **matrix = new double *[n];

       if( !matrix )
       {
           cout << "!Blad alokacji!";
           exit(1);
       }
       for( int i = 0; i < n; i++ )
       {
           matrix[i] = new double [m];
           if( !matrix[i] )
           {
               cout << "!Blad alokacji!";
               exit(1);
           }
       }
       return matrix;
}

void usun_macierz( double **matrix, int n ) //n - liczba wierszy
{
     for( int i = n-1; i > 0; i-- )
          delete []matrix[i];
     delete []matrix;
}

void usun_wektor( double *x )
{
     delete []x;
}

void wypisz_wektor( double *x, int k ) // x - wektor, k - kolumn
{
     for( int i = 0; i < k; i++ )
          cout << x[i] << " | ";
     cout << endl;
}

void wypisz_macierz( double **x, int r, int c ) // c - wiersze, r - kolumny
{
     for( int i = 0; i < r; i++ )
     {
          for( int j = 0; j < c; j++ )
               cout << x[i][j] << "    ";
          cout << endl;
     }
}


void fill_dt(double *x, double delta, int k)
{
     x[0] = 0; // wynika z zalozen
     for(int i=1; i<k; i++)
          x[i]=x[i-1]+delta;
}

void fill_h(double *x, double delta, int k)
{
     x[0] = -a; // wynika z warunkow poczatkowych
     for(int i=1; i<k; i++)
          x[i]=x[i-1]+delta;
}


void zapisz_plik( char* filename, double **matrix, int r, int c)
{
     ofstream file( filename );
     if( !file.is_open() )
     {
         cout << "Nie moge zapisac: " << filename <<" !!!"<< endl;
         exit(1);
     }
     file.setf(ios::scientific, ios::floatfield);
     file.precision(4);
     file << endl << "                " << filename << ":                " << endl << endl;
     for( int i = 0; i < r; i++ )
     {
          file << endl;
          for( int j = 0; j < c; j++ )
              file << matrix[i][j]<<"\t";
     }
     file.close();
     cout << "Zapisano: " << filename  << endl;
}

double analityczne(double **ANAL, int N, int M, double *dt_levels, double *h_levels, double dt, double h)
{
     for(int i=0; i<N; i++)    //warunek brzegowy U(0,t)=0
     ANAL[i][0] = 1.;
     for(int i=0; i<N; i++)    //warunek brzegowy U(INF,t)=0
     ANAL[i][M-1]=0.;
     for(int i=0; i<M; i++)  //warunek poczatkowy  U(x,0)= (1 dla x<0) (1 dla x>=0)
        if(h_levels[i]<0.) ANAL[0][i]=1.;
                    else ANAL[0][i]=0.;

       for(int i=1; i<N; i++)
       {
            for( int j=1; j<M; j++)
                 ANAL[i][j]=0.5*erfc((h_levels[j])/(2*sqrt(D*dt_levels[i])));
       }
       zapisz_plik("an_CN.txt", ANAL, N, M);
}

void Thomas( double *l, double *d, double *u, double *b, double *x, int n )  // l,d,u,b - wektory danych. x - wektor wynikowy , n - dlugosc wektorow
{
       int i;

       // modyfikacja wspolczynnikow

       u[0] = u[0] / d[0];                           // ryzyko dzielenia przez 0 /
       b[0] = u[0] / d[0];                           // dzielenei przez 0 oznaczaloby pojedyncza macierz (mac. prosta?)
       for( i = 1; i < n; i++ )
   {
               double mianownik = (d[i] - u[i-1] * l[i]);     // ryzyko dzielnia przez 0
               u[i] = u[i] / mianownik;                             //ostatnia wartosc jest liczona nadmiarowo
               b[i] = ( b[i] - b[i-1] * l[i] ) / mianownik;
       }
       // zwroc wynik
       x[n - 1] = b[n - 1];
       for( i = n - 2; i >= 0; i-- )
               x[i] = b[i] - u[i] * x[i+1];
}
void dekompozycja_LU(double **A, double *b, double *wyn)
{
	//Dekompozycja LU macierzy- metoda eliminacji Gaussa
	double x;
	for(int k=0; k<M-1; k++)
	{
		for(int i=k+1; i<M; i++)
		{
			x = A[i][k]/A[k][k];
			A[i][k] = x;
			for(int j=k+1; j<M; j++)
			{
				A[i][j] = A[i][j] - (x*A[k][j]);
			}
		}
	}

	//Rozwi¹zywanie uk³adu równañ
	double suma;
	double *z = new double[M];

	//podstawianie w przód
	for(int i=0; i<M; i++)
	{
		suma = 0;
		for(int j=0; j<=i-1; j++)
		{
			suma += A[i][j]*z[j];
		}

		z[i] = b[i]-suma;
	}

	 //podstawianie w ty³
	for(int i=M-1; i>=0; i--)
	{
		suma = 0;
		for(int j=i+1; j<M; j++)
		{
			suma +=A[i][j]*wyn[j];
		}

		wyn[i] = (z[i]-suma)/A[i][i];
	}
}
void blad(double **BLAD, double **ANAL, double **U, int N, int M)
{
    for(int i=0; i<N; i++)
    {
            for(int j=0; j<M; j++)
            {
                BLAD[i][j]=fabs(U[i][j]-ANAL[i][j]);
            }
            }
    zapisz_plik("bledy_CN.txt", BLAD, N, M);

}
/*
void CrankNicolson( double **U, int N, int M, double *dt_levels, double *h_levels, double dt, double h, double )
{
     double *a=nowy_wektor(M);
     double *d=nowy_wektor(M);
     double *c=nowy_wektor(M);
     double *u=nowy_wektor(M);
     double *b=nowy_wektor(M);

     for(int i=0; i<N; i++)    //warunek brzegowy Xm
     U[i][0] = 1.;
     for(int i=0; i<N; i++)    //warunek brzegowy X0
     U[i][M-1]=0.;
     for(int i=0; i<M; i++)  //warunek poczatkowy  U(x,0)= (1 dla x<0) (1 dla x>=0)
        if(h_levels[i]<0.) U[0][i]=1.;
                    else U[0][i]=0.;


     for(int k=1; k<N; k++)
     {
          a[0]=0.0;
          d[0]=1.0;
          c[0]=0.0;
          b[0]=U[k-1][0];

          for(int i=1; i<M-1; i++)
          {
               a[i]=lambda/2.0;
               d[i]=-(1+lambda);
               c[i]=lambda/2.0;
               b[i]=-(lambda/2.0*U[k-1][i-1]+(1.0-lambda)*U[k-1][i]+(lambda/2.0)*U[k-1][i+1]);
          }

          a[M-1]=0.0;
          d[M-1]=1.0;
          c[M-2]=lambda/2.0;
          b[M-1]=U[k-1][M-1];

          Thomas(a, d, c, b, u, M);
          for(int i=1; i<M-1; i++)
               U[k][i]=u[i];
     }
     zapisz_plik("num_CN.txt", U, N, M );
}

*/

void LaasonenThomas(double **U, int N, int M, double *dt_levels, double *h_levels, double dt, double h,double lambda_l){
     double *a=nowy_wektor(M);
     double *d=nowy_wektor(M);
     double *c=nowy_wektor(M);
     double *u=nowy_wektor(M);
     double *b=nowy_wektor(M);

     for(int i=0; i<N; i++)    //warunek brzegowy Xm
     U[i][0] = 1.;
     for(int i=0; i<N; i++)    //warunek brzegowy X0
     U[i][M-1]=0.;
     for(int i=0; i<M; i++)  //warunek poczatkowy  U(x,0)= (1 dla x<0) (1 dla x>=0)
        if(h_levels[i]<0.) U[0][i]=1.;
                    else U[0][i]=0.;

    /*
    double **A = nowa_macierz(N,M);
    for(int i=0;i<N;i++) for(int j=0;j<M;j++) A[i][j]=0;
    */
    for( int k = 1; k < N; k++ )
     {
		   a[0]=0.0;
          d[0]=1.0;
          c[0]=0.0;
          b[0]=U[k-1][0];

          for( int i = 1; i < M-1; i++ )
		  {
               a[i]=lambda/2.0;
               d[i]=-(1+lambda);
               c[i]=lambda/2.0;
               b[i]=-(lambda/2.0*U[k-1][i-1]+(1.0-lambda)*U[k-1][i]+(lambda/2.0)*U[k-1][i+1]);
		  }

          a[M-1]=0.0;
          d[M-1]=1.0;
          c[M-2]=lambda/2.0;
          b[M-1]=U[k-1][M-1];

          Thomas(a, d, c, b, u, M);
          for(int i=1; i<M-1; i++)
               U[k][i]=u[i];
     }
     zapisz_plik("num_MPL_thomas.txt", U, N, M );






}
void LaasonenLU(double **U, int N, int M, double *dt_levels, double *h_levels, double dt, double h,double lambda_l)
{
    double *b=nowy_wektor(M);
    double *wyn=nowy_wektor(M);
    double **A = nowa_macierz(N,M);
    for(int i=0;i<N;i++) for(int j=0;j<M;j++) A[i][j]=0;

    for(int i=0; i<N; i++)    //warunek brzegowy Xm
     U[i][0] = 1.;
     for(int i=0; i<N; i++)    //warunek brzegowy X0
     U[i][M-1]=0.;
     for(int i=0; i<M; i++)  //warunek poczatkowy  U(x,0)= (1 dla x<0) (1 dla x>=0)
        if(h_levels[i]<0.) U[0][i]=1.;
                    else U[0][i]=0.;

    for( int k = 1; k < N; k++ )
     {

          A[0][0]=1.0;
          b[0]=U[k-1][0];

          for( int i = 1; i < M-1; i++ )
		  {
               A[i][i-1]=lambda/2.0;
               A[i][i]=-(1+lambda);
               A[i][i+1]=lambda/2.0;
               b[i]=-(lambda/2.0*U[k-1][i-1]+(1.0-lambda)*U[k-1][i]+(lambda/2.0)*U[k-1][i+1]);
		  }

          A[N-2][M-1]=0.0;
          A[N-1][M-1]=1.0;
          A[N-1][M-2]=lambda/2.0;
          b[M-1]=U[k-1][M-1];


        dekompozycja_LU(A, b, wyn);
        for(int i=1; i<M-1; i++)
               U[k][i]=b[i];


        zapisz_plik("num_MPL_LU.txt", U, N, M );



}


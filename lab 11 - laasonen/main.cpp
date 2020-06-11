#include <fstream>
#include<iostream>
#include <iomanip>
#include<cmath>
#include <math.h>
#include<fstream>
#include "calerf.h"

using namespace std;
//====================================== Sta³e ========================================
    const double D=1.0;
    const double tmax=2.0;
    const double a=6*sqrt(D*tmax);
    const double lambda=1;
    const double r=1.0;
    const double h=0.1; //krok przestrzenny
    const double dt=(lambda*(h*h)/D);
    const int N=static_cast<const int>((tmax/dt)+2); //ilosc poziomow czasowych /
    const int M = static_cast<const int>((a*2/h)+1); //ilosc poziomow przestrz. /
    double blad_wykresy = 0.0;


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
    void blad(double **BLAD, double **ANAL, double **U, int N, int M,char* filename);
    double maxblad(double *A, int n, int m,double h,int t,double *h_levels,double *dt_levels);
    double maxblad1(double *BLAD, int n, int m);

//=====================================================================================
int main()
{
    //double blad_wykresy=0.0;
    double last_blad;
    double last_h;
    double h2;
    cout.setf( ios::fixed, ios::floatfield );
    cout.precision(1);

    double *dt_levels = nowy_wektor( N ); // wektor poziomow czasowych
    double *h_levels = nowy_wektor( M ); // wektor poziomow przestrzennych
    fill_dt(dt_levels, dt, N);
    fill_h(h_levels, h, M);
    char *filename2 = "Laasonen.txt";


    ofstream fs2( filename2 );
    if( !fs2.is_open() )
    {
        cout << "Nie moge zapisac: " << filename2 <<" !!!"<< endl;
        exit(1);
    }

   //file.close();
   // cout << "Zapisano: " << filename  << endl;

    //--------------------------------------------------------------------------------------

    double **ANAL = nowa_macierz( N, M );
    analityczne( ANAL, N, M, dt_levels, h_levels, dt, h );


    double **U = nowa_macierz( N, M );
    LaasonenThomas(U,N,M,dt_levels, h_levels,dt,h,lambda);

    double **A = nowa_macierz(N,M);
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<M;j++)
        {
            A[i][j]=0.0;
            //cout<<A[i][j]<<" ";
        }
        //cout<<endl;
    }

    LaasonenLU(A,N,M,dt_levels,h_levels,dt,h,lambda);


    double **BLAD = nowa_macierz(N,M);
    blad(BLAD, ANAL, U, N, M,"blad_laasonen_Thomas.txt");

    double **BLAD2 = nowa_macierz(N,M);
    blad(BLAD2, ANAL, A, N, M,"blad_laasonen_LU.txt");

    /** Laasonen Thomas*/


    //zapis wynikow Laasonen - Thomas do pliku
    fs2 << setprecision(15);
    fs2 << setw(25) << left << "x" << '\t';
    fs2 << setw(25) << left << "t" << '\t';
	fs2 << setw(25) << left << "Ux,t" << '\t';
	//fs2 << setw(25) << left << "Ax,t" << '\t';
	fs2 << setw(25) << left << "Analityczne" << '\t';
	fs2 << setw(25) << left << "BladThomas" << '\n';
	//fs2 << setw(25) << left << "BladLU" << endl;

	for(int i=0;i<N;i++){
        for(int j=0;j<M;j++){
            fs2<<setw(25)<<left<<h_levels[j]<<setw(25)<<left<<dt_levels[i]<<setw(25)<<left<<U[i][j]<<'\t'/*<<setw(25)<<left<<A[i][j]<<'\t'*/
            <<setw(25)<<left<<ANAL[i][j]<<'\t'<< scientific<<setw(25)<<left<<BLAD[i][j]/*<<'\t'<<setw(25)<<left<<BLAD2[i][j]<<'\t'*/<<endl<< fixed;
        }
	}

    cout << "Zapisano: " << filename2  << endl;
    fs2.close();

    //Obliczenie maksymalnego bezwzglêdnego b³êdu dla ka¿dego kroku czasowego t z zapisem do pliku//
    fs2.open("bladt_lassonen_thomas.txt", fstream::out);
    for (int i = 1; i < N; i++) fs2 << setw(25) << left << dt_levels[i]<< '\t' << setw(25) << left << maxblad1(BLAD[i], N, M) << endl;
	fs2.close();
    cout<<"Zapisano bladt_lassonen_thomas.txt"<<endl;




    //Zapis do pliku wartoœci funkcji uzyskanej i wartoœci analityczne dla 3 ró¿nych chwili czasu t
	fs2.open("met_lassonenwyn_thomas.txt", fstream::out);
	fs2 << "#t1 = " << dt_levels[5] << endl;
	fs2 << "#t2 = " << dt_levels[50] << endl;
	fs2 << "#t3 = " << dt_levels[100] << endl;
	for (int i = 0; i < M; i++) {
		fs2 << setw(25) << left << h_levels[i]<< '\t' << setw(25) << left << U[5][i] << '\t' << setw(25) << left << ANAL[5][i] << '\t';
		fs2 <<  '\t' << setw(25) << left << U[50][i] << '\t' << setw(25) << left << ANAL[50][i] << '\t';
		fs2 <<  '\t' << setw(25) << left << U[100][i] << setw(25) << left << ANAL[100][i] << endl;
	}
	cout<<"Zapisano: met_lassonenwyn_thomas.txt"<<endl;
	fs2.close();
    usun_macierz( U, N );
	//Obliczanie maksymalnego bezwzglêdnego b³êdu dla czasu tmax w zale¿noœci od kroku h dla metody Lassonen. Uruchamiam algorytm
	//zmienajaæ iloœæ wez³ów na osi x
    /*
    h2=0.1;
    blad_wykresy=0.0;
    double dt2;

        //Obliczanie maksymalnego bezwzglêdnego b³êdu dla czasu tmax w zale¿noœci od kroku h dla metody Lassonen. Uruchamiam algorytm
	//zmienajaæ iloœæ wez³ów na osi x

	fs2 << setprecision(15);
	fs2.open("blad_lassonen_thomas.txt", fstream::out);
	for (int i = 100; i > 0; i = i - 5) {
		int N1 = 400;
		int M1 = i;
		last_blad = blad_wykresy;
		last_h = h2;
		h2 = fabs((-a) / ((double)M1 - 1.0 ));
		double *h_vector = nowy_wektor( M1 );
		fill_h(h_vector,h2,M1);
		dt2=(lambda*(h2*h2)/D);
		double *dt_vector = nowy_wektor( N1 );
		fill_dt(dt_vector,dt2,M1);
		double **U=nowa_macierz(N1,M1);
		for(int i=0; i<N1; i++)    //warunek brzegowy Xm
        U[i][0] = 1.;
     for(int i=0; i<N1; i++)    //warunek brzegowy X0
        U[i][M1-1]=0.;
     for(int i=0; i<M1; i++)  //warunek poczatkowy  U(x,0)= (1 dla x<0) (1 dla x>=0)
        if(h_vector[i]<0.) U[0][i]=1.;
                    else U[0][i]=0.;

		LaasonenThomas(U,N1,M1,dt_levels,h_vector,dt,h2,lambda);

		double blad_wykresy= maxblad(U[N1-1], N1, M1,h2,(N1-1),h_vector,dt_vector);
		cout<<"h2= "<<h2<<" blad_wykresy = "<<blad_wykresy<<endl;
		fs2 << setw(15) << left << log10(h2) <<'\t'<< setw(15) << left << log10(blad_wykresy) << endl;
        cout<<"h2 logarytm= "<<log10(h2)<<" blad_wykresy logarytm= "<<log10(blad_wykresy)<<endl;
		usun_macierz( U, N1 );
		usun_wektor(h_vector);
		usun_wektor(dt_vector);

	}
	cout<<"Zapisano: blad_lassonen_thomas.txt"<<endl;
	fs2.close();
    */

	/**Laasonen - LU */


	 //zapis wynikow Laasonen - Thomas do pliku
    fs2 << setprecision(15);
    fs2 << setw(25) << left << "x" << '\t';
    fs2 << setw(25) << left << "t" << '\t';
	fs2 << setw(25) << left << "Ax,t" << '\t';
	//fs2 << setw(25) << left << "Ax,t" << '\t';
	fs2 << setw(25) << left << "Analityczne" << '\t';
	fs2 << setw(25) << left << "BladLU" << '\t';
	//fs2 << setw(25) << left << "BladLU" << endl;

	for(int i=0;i<N;i++){
        for(int j=0;j<M;j++){
            fs2<<setw(25)<<left<<h_levels[j]<<setw(25)<<left<<dt_levels[i]<<setw(25)<<left<<A[i][j]<<'\t'/*<<setw(25)<<left<<A[i][j]<<'\t'*/
            <<setw(25)<<left<<ANAL[i][j]<<'\t'<< scientific<<setw(25)<<left<<BLAD2[i][j]/*<<'\t'<<setw(25)<<left<<BLAD2[i][j]<<'\t'*/<<endl<< fixed;
        }
	}

    cout << "Zapisano: " << filename2  << endl;
    fs2.close();

    //Obliczenie maksymalnego bezwzglêdnego b³êdu dla ka¿dego kroku czasowego t z zapisem do pliku//
    fs2.open("bladt_lassonen_LU.txt", fstream::out);
    for (int i = 1; i < N; i++) fs2 << setw(25) << left << dt_levels[i]<< '\t' << setw(25) << left << maxblad1(BLAD2[i], N, M) << endl;
	fs2.close();
    cout<<"Zapisano bladt_lassonen_LU.txt"<<endl;




    //Zapis do pliku wartoœci funkcji uzyskanej i wartoœci analityczne dla 3 ró¿nych chwili czasu t
	fs2.open("met_lassonenwyn_LU.txt", fstream::out);
	fs2 << "#t1 = " << dt_levels[5] << endl;
	fs2 << "#t2 = " << dt_levels[50] << endl;
	fs2 << "#t3 = " << dt_levels[100] << endl;
	for (int i = 0; i < M; i++) {
		fs2 << setw(25) << left << h_levels[i]<< '\t' << setw(25) << left << A[5][i] << '\t' << setw(25) << left << ANAL[5][i] << '\t';
		fs2 <<  '\t' << setw(25) << left << A[50][i] << '\t' << setw(25) << left << ANAL[50][i] << '\t';
		fs2 <<  '\t' << setw(25) << left << A[100][i] << setw(25) << left << ANAL[100][i] << endl;
	}
	cout<<"Zapisano: met_lassonenwyn_LU.txt"<<endl;
	fs2.close();
    //usun_macierz( A, N );

	//Obliczanie maksymalnego bezwzglêdnego b³êdu dla czasu tmax w zale¿noœci od kroku h dla metody Lassonen. Uruchamiam algorytm
	//zmienajaæ iloœæ wez³ów na osi x

   // h2=0.1;
    //blad_wykresy=0.0;
        //Obliczanie maksymalnego bezwzglêdnego b³êdu dla czasu tmax w zale¿noœci od kroku h dla metody Lassonen. Uruchamiam algorytm
	//zmienajaæ iloœæ wez³ów na osi x
	/*
	fs2 << setprecision(15);
	fs2.open("blad_lassonen_LU.txt", fstream::out);
	for (int i = 100; i > 0; i = i - 5) {
		int N1 = 400;
		int M1 = i;
		last_blad = blad_wykresy;
		last_h = h2;
		h2 = (- a) / ((double)M1 - 1.0 );
		double *h_vector = nowy_wektor( N1 );
		fill_h(h_vector,h2,M1);
		double **A=nowa_macierz(N1,M1);
		for(int i=0; i<N1; i++)    //warunek brzegowy Xm
     U[i][0] = 1.;
     for(int i=0; i<N1; i++)    //warunek brzegowy X0
     U[i][M1-1]=0.;
     for(int i=0; i<M1; i++)  //warunek poczatkowy  U(x,0)= (1 dla x<0) (1 dla x>=0)
        if(h_vector[i]<0.) A[0][i]=1.;
                    else A[0][i]=0.;

		LaasonenLU(A,N1,M1,dt_levels,h_vector,dt,h,lambda);

		blad_wykresy = maxblad1(A[N1-1], N1, M1);
		fs2 << setw(15) << left << log10(h2) <<'\t'<< setw(15) << left << log10(blad_wykresy) << endl;

		usun_macierz( A, N1 );
		usun_wektor(h_vector);

	}

	cout<<"Zapisano: blad_lassonen_LU.txt"<<endl;
	fs2.close();
    */















    usun_macierz( ANAL, N );

    usun_macierz(A,N);
    usun_macierz(U,N);
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
     //file << endl << "                " << filename << ":                " << endl << endl;
     for( int i = 0; i < r; i++ )
     {
            file<<endl;
          for( int j = 0; j < c; j++ )
              file << matrix[i][j]<<" ";
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
       zapisz_plik("an_laasonen.txt", ANAL, N, M);
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
void blad(double **BLAD, double **ANAL, double **U, int N, int M, char* filename)
{
    for(int i=0; i<N; i++)
    {
            for(int j=0; j<M; j++)
            {
                BLAD[i][j]=fabs(U[i][j]-ANAL[i][j]);
            }
            }
    //zapisz_plik(filename, BLAD, N, M);

}


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
               a[i]=lambda;
               d[i]=-(1+(2*lambda));
               c[i]=lambda;
               b[i]=-(lambda*U[k-1][i-1]+(1.0-(2*lambda))*U[k-1][i]+(lambda)*U[k-1][i+1]);
		  }

          a[M-1]=0.0;
          d[M-1]=1.0;
          c[M-2]=lambda;
          b[M-1]=U[k-1][M-1];

          Thomas(a, d, c, b, u, M);
          for(int i=1; i<M-1; i++)
               U[k][i]=u[i];
     }
     //zapisz_plik("num_MPL_thomas.txt", U, N, M );






}
void LaasonenLU(double **A, int N, int M, double *dt_levels, double *h_levels, double dt, double h,double lambda_l)
{

    double *b=nowy_wektor(M);
    double *wyn=nowy_wektor(M);
    double **U = nowa_macierz(N,M);
    for(int i=0;i<N;i++) for(int j=0;j<M;j++) A[i][j]=0;

    for(int i=0; i<N; i++)    //warunek brzegowy Xm
     A[i][0] = 1.;
     for(int i=0; i<N; i++)    //warunek brzegowy X0
    A[i][M-1]=0.;
     for(int i=0; i<M; i++)  //warunek poczatkowy  U(x,0)= (1 dla x<0) (1 dla x>=0)
        if(h_levels[i]<0.) A[0][i]=1.;
                    else A[0][i]=0.;

    for( int k = 1; k < N; k++ )
     {

          U[0][0]=1.0;
          b[0]=A[k-1][0];

          for( int i = 1; i < M-1; i++ )
		  {
               U[i][i-1]=lambda;
               U[i][i]=-(1+(2*lambda));
               U[i][i+1]=lambda;
               b[i]=-(lambda*A[k-1][i-1]+(1.0-2*lambda)*A[k-1][i]+(lambda)*A[k-1][i+1]);
		  }

          U[N-1][M-2]=0.0;
          U[N-1][M-1]=1.0;
          U[N-2][M-1]=lambda;
          b[M-1]=A[k-1][M-1];


        dekompozycja_LU(U, b, wyn);
        for(int i=1; i<M-1; i++)
               //U[k][i]=b[i];
               A[k][i]=wyn[i];

     }
        //zapisz_plik("num_MPL_LU.txt", U, N, M );



}

double maxblad1(double *BLAD, int n, int m)
{
	double maxblad = 0.0;

	for (int j = 0; j < m; j++)
    {

    if(BLAD[j] > maxblad)
            maxblad = BLAD[j];
	}
	return maxblad;
}

double maxblad(double *A, int n, int m,double h,int t,double *h_levels,double *dt_levels)
{
	double maxblad = 0.0;
    double blad = 0.0;
	for (int j = 0; j < m; j++)
    {
		blad = fabsl(0.5*erfc((h_levels[j])/(2*sqrt(D*dt_levels[t]))) - A[j]);
		if (blad > maxblad)
            maxblad = blad;
	}
	return maxblad;
}

#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <cfloat>
#include <limits>

using namespace std;
double funkcja_exp(double x)
{
    if(x==0) return -1;
	return (1.0 - exp(-x)) / x;
}
double szeregtaylora(double x){
	/*
		 *  1) 1.0
		 *  2) x^1/silnia(2)
		 *  3) x^2/silnia(3)
		 *  4) x^3/silnia(4)
	 */
	double wynik = 1.0;
	double wyraz = 1.0;

	for (int i =2.0; i < 100.0; i++){
		wyraz *= -1*(x /(double)i);
		wynik += wyraz;

            cout<<"x: "<<x<<"   Wyraz: "<<wyraz<<"     "<<"Wynik: "<<wynik<<endl;
	}

	return wynik;
}


int main(){

	double x, wynik, wartosc_dok,blad_bezw,blad_wz, temp, lg,taylor, taylor_wzg, taylor_bezwz;
    string odstep = "\t\t";
	ifstream plik;
	ofstream zapis;
	fstream WYKRES;
    fstream WYKRES_TAYLOR;

	plik.open("dane.txt");
	zapis.open("wynik.txt");
	WYKRES.open("rysuj.dat", ios::out);
	WYKRES_TAYLOR.open("taylor.dat", ios::out);

	if (!plik.good() || !zapis.good() || !WYKRES.good() || !WYKRES_TAYLOR.good()) return EXIT_FAILURE;

	zapis << setprecision(21);
	WYKRES << setprecision(21);
    WYKRES_TAYLOR << setprecision(21);
	cout << setprecision(21);
	cout.setf(ios::scientific);
	zapis.setf(ios::scientific);
/*
    cout<<"x:\t\t\t\t\twartosc_dokladna:\t\t\t1.0-exp(-x))/x:\t\t\tblad bezwzgledny: \t\t\tblad wzgledny:"<<endl;
    zapis << odstep <<"log10" << odstep << odstep << "x" << odstep <<odstep  << "w_dokladny"
     << odstep << "f(x)" << odstep << odstep << "blad bezwgledny" << odstep
       << "blad wzgledny" <<odstep << "szereg" << odstep<<odstep << "roznica_szereg" << odstep << "blad wzgledny" << endl;
*/
     zapis << odstep <<"log10"
     << odstep << odstep << "x"
     << odstep <<odstep  << "w_dokladny"
      << odstep << "f(x)" << odstep
      << odstep << "blad bezwgledny" << odstep  << "blad wzgledny" <<odstep
      << "szereg taylora" << odstep<<odstep << "taylor b. bezwzgl" << odstep << "taylod b. wzgledny" << endl;

	while (!plik.eof()){



//		pobranie wartosci z pliku
		plik >> lg; //ln(x)
		plik >> x; //pobranie x
		plik >> wartosc_dok; //wartosc dokladna
		wynik=funkcja_exp(x);
        blad_bezw=fabs(wartosc_dok-wynik);
        blad_wz=fabs(blad_bezw/wartosc_dok);




        zapis <<setw(30)<< lg;
        zapis <<setw(30)<< x;
        zapis <<setw(30)<< wartosc_dok;
        zapis <<setw(30)<< wynik;
        zapis <<setw(30)<< blad_bezw;

         if(wartosc_dok!=0)
        {
           blad_wz=fabs(blad_bezw/wartosc_dok);
            //cout << setw(szerokosc) << wzgledna;
            zapis <<setw(30)<< blad_wz;

        WYKRES << lg << " ";
        WYKRES << log10(blad_wz) <<" ";
        }

        if(blad_bezw > numeric_limits<double>::epsilon())
        {
            taylor=szeregtaylora(x);
            //cout << setw(szerokosc) << taylor;
            zapis << setw(30) << taylor;

            taylor_bezwz=fabs(taylor-wartosc_dok);

            //cout << setw(szerokosc) << bezwgledna_sz;
            zapis << setw(30) << taylor_bezwz;

          if(wartosc_dok!=0)
        {
            taylor_wzg=fabs(taylor_bezwz/wartosc_dok);
            //cout << setw(szerokosc) << wzgledna;
            zapis <<setw(30)<< taylor_wzg;
            WYKRES_TAYLOR << lg << " ";
        WYKRES_TAYLOR << log10(fabs(taylor_wzg))<<" " ;

        }

        }




        zapis << "\n";
        WYKRES << "\n";
        WYKRES_TAYLOR << "\n";









		}
		//cout<<DBL_EPSILON;

		plik.close();
		zapis.close();

        WYKRES.close();
        WYKRES_TAYLOR.close();




	return 0;
}


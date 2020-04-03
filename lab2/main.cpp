#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <cfloat>
using namespace std;

double f1(double x)
{
    if(x==0)
        return -1;
    double zmienna=(1 - exp(-x))/x;
    return zmienna;
}
double taylorSeries(double x)
{
    double wynik = 1.0, wyraz = 1.0;

	for (int i = 1 ; i < 30 ; i++)
    {
		wyraz *= -(x / ((double)i + 1.0));
		wynik += wyraz;
	}
	return wynik;
}

int main()
{
    fstream PLIK;
    fstream WYNIK;
    fstream wykres;
    fstream wykres_taylor;
    string odstep="\t\t";
    const int szerokosc = 30;
    double zmienna,x,dokladny, wynik, bezwgledna, funkcja, lg, wzgledna,taylor,bezwgledna_sz, wzgledna_sz;


    PLIK.open("dane.txt", ios::in);
	WYNIK.open("wynik.txt", ios::out);
	wykres.open("rysuj.dat", ios::out);
	wykres_taylor.open("taylor.dat", ios::out);



	if (!PLIK.good() || !WYNIK.good() || !wykres.good() || !wykres_taylor.good())
    {
        cout<<"Blad otwarcia pliku"<<endl;
        return -1;
    }

    //Formatowanie tekstu
    WYNIK << setprecision(21);
    wykres << setprecision(21);
    wykres_taylor << setprecision(21);
	cout << setprecision(21);
	cout.setf(ios::scientific);
	WYNIK.setf(ios::scientific);

	WYNIK << odstep <<"log10" << odstep << odstep << "x" << odstep <<odstep  << "exact value" << odstep << "f1(x)" << odstep << odstep << "absolute error" << odstep  << "relative error" <<odstep << "szereg" << odstep<<odstep<< "roznica_szereg" << odstep << "blad wzgledny" << endl;

    while(!PLIK.eof())
    {
        PLIK >> lg;
		WYNIK << setw(szerokosc) << lg;
        PLIK >> x;
		WYNIK << setw(szerokosc) << x;
        PLIK >> dokladny;
		WYNIK << setw(szerokosc) << dokladny;

		funkcja=f1(x);
		WYNIK << setw(szerokosc) << f1;


		bezwgledna=fabs(funkcja-dokladny);
		WYNIK << setw(szerokosc) << bezwgledna;


        if(dokladny!=0)
        {
            wzgledna=fabs(bezwgledna/dokladny);
            WYNIK << setw(szerokosc) << wzgledna;


        }


        wykres << lg << " ";
        wykres << log10(wzgledna) <<" ";


        if(fabs(bezwgledna) >  DBL_EPSILON)
        {
            taylor=taylorSeries(x);
            WYNIK << setw(szerokosc) << taylor;
            bezwgledna_sz=fabs(taylor-dokladny);
            WYNIK << setw(szerokosc) << bezwgledna_sz;


            if(dokladny!=0)
            {
                wzgledna=fabs(bezwgledna_sz/dokladny);
                WYNIK << setw(szerokosc) << wzgledna_sz;

            }
        }
        wykres_taylor << lg << " ";
        wykres_taylor << log10(wzgledna)<<" \n" ;
        WYNIK << "\n";
        wykres << "\n";
        }

    cout<<DBL_EPSILON;
    PLIK.close();
	WYNIK.close();
	wykres.close();
	wykres_taylor.close();

    return 0;
}


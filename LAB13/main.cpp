#include <iostream>
#include <math.h>
#include <cmath>
#include <iomanip>
#include <fstream>
#define _USE_MATH_DEFINES
#define PI 3.14159265359
#define liczba_przedzialow 100000
#include "calerf.h"
#include "calerf.cpp"
using namespace std;

long double fun(double y);
void print();
long double prostokatyA(double x, double step);
long double prostokatyB(double x, double step);
long double prostokatyC(double x, double step);
long double kwadratura_trapezow(double x, double step);
long double kwadratura_parabol(double x, double step);
long double erf_podane(double x);

int main()
{
    int i;
    const int szerokosc = 30;
    long double x=1.0;
    cout<<scientific;
    fstream wykres1;
    //fstream wykres2;
    wykres1.open("wyniki.txt", ios::out);

    for(i = 0 ; i<3 ; i++)
    {
        cout<<"x="<<x<<endl<<endl;
        print();
        for (long double j = 0.5 ; j>10e-6 ; j/=2)
        {
            cout.width(szerokosc);
            cout<<left<<j;
            cout.width(szerokosc);
            cout<<left<<erf_podane(prostokatyA(x,j));
            cout.width(szerokosc);
            cout<<left<<erf_podane(prostokatyB(x,j));
            cout.width(szerokosc);
            cout<<left<<erf_podane(prostokatyC(x,j));
            cout.width(szerokosc);
            cout<<left<<erf_podane(kwadratura_trapezow(x,j));
            cout.width(szerokosc);
            cout<<left<<erf_podane(kwadratura_parabol(x,j))<<endl;
        }
        x=x+1.0;
        cout<<endl;

    }


    long double h;
    //long double a;
    long double wartosc_dokladna = calerf::ERFL(3.0);
    long double bladA, bladB, bladC, bladParabole,bladTrapezy;
    for(h = 0.5 ; h>10e-6 ; h/=2)
    {
        bladA = fabs(wartosc_dokladna - erf_podane(prostokatyA(3.0,h))) / wartosc_dokladna;
        bladB = fabs(wartosc_dokladna - erf_podane(prostokatyB(3.0,h))) / wartosc_dokladna;
        bladC = fabs(wartosc_dokladna - erf_podane(prostokatyC(3.0,h))) / wartosc_dokladna;
        bladTrapezy = fabs(wartosc_dokladna - erf_podane(kwadratura_trapezow(3.0,h))) / wartosc_dokladna;
        bladParabole = fabs(wartosc_dokladna - erf_podane(kwadratura_parabol(3.0,h))) / wartosc_dokladna;

        wykres1<<log10(h)<<" "<<log10(bladA)<<" "<<log10(bladB)<<" "<<log10(bladC)<<" "<<log10(bladTrapezy)<<" "<<log10(bladParabole)<<endl;
        //wykres2<<
    }
    return 0;
}



long double fun(double y)
{
    return exp(-y*y);
}

long double prostokatyA(double x, double step)
{
    long double sum=0;
    long double i;

    for(i=0 ; i<x ; i+=step)
    {
        sum += fun(i) * step;
    }
    return sum;
}

long double prostokatyB(double x, double step)
{
    long double sum=0;
    long double i;

    for(i=0 ; i<x ; i+=step)
    {
        sum += fun(i+step)*step;
    }
    return sum;
}

long double prostokatyC(double x, double step)
{
    long double sum=0;
    long double i;

    for(i=0 ; i<x ; i+=step)
    {
        sum += fun(i+step/2.0) * step;
    }
    return sum;
}
long double erf_podane(double x)
{
    return 2 / sqrt(PI) * x;
}


long double kwadratura_trapezow(double x, double step)
{
    long double sum=0;
    long double i;

    for(i=0 ; i<x ; i+=step)
    {
        sum += (fun(i) + fun(i+step)) /2.0*step;
    }
    return sum;
}

long double kwadratura_parabol(double x, double step)
{
    long double sum=0;
    long double i;
    //int a;
    for(i=0 ; i<x ; i+=step)
    {
        //
        sum += (1.0/6.0*fun(i) + 4.0/6.0*fun(i+step/2.0) + 1.0/6.0*fun(i+step))  * step;
    }
    return sum;
}
void print()
{
    const int szerokosc = 30;
    cout.width(szerokosc);
    cout<<left<<"h";
    cout.width(szerokosc);
    cout<<left<<"Prostokaty A";
    cout.width(szerokosc);
    cout<<left<<"Prostokaty B";
    cout.width(szerokosc);
    cout<<left<<"Prostokaty C";
    cout.width(szerokosc);
    cout<<left<<"Trapezy";
    cout.width(szerokosc);
    cout<<left<<"Parabole"<<endl;

}

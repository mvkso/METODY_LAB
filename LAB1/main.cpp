#include <iostream>
#include <conio.h>
#include<cfloat>

/*
Napisz program w jêzyku „C/C++”, umo¿liwiaj¹cy „doœwiadczalne” wyznaczenie liczby bitów
mantysy oraz tzw. epsylona maszynowego, dla zmiennych typu float i double, tj. najmniejszej liczby
e takiej, ze fl(e + 1) > 1. Jaki jest zwiazek e z precyzja arytmetyki?
*/

using namespace std;

int main()
{
    int m_f=0,m_d=0; // liczba bitów mantysy dla float i double
    float e_f=1.f; //float epsilon
    double e_d=1.L; //double epsilon
    float a=.5f*e_f + 1.f;
    double b=.5D*e_d + 1.D;
    // FLOAT
    cout<<"LICZE EPSILON TYPU FLOAT: "<<endl<<endl;
	while ((a) > 1.f) {
		e_f*=.5f;
		cout<<"e_f: "<<e_f<<endl;
		m_f++;
		a=.5f*e_f + 1.f;
	}

    cout<<endl<<"Liczba bitow mantysy dla float: "<<m_f<<endl;
	cout<<"Epsilon typu float: "<<e_f<<endl;

	// DOUBLE:
	cout<<"LICZE EPSILON TYPU DOUBLE: "<<endl<<endl;
	e_d = 1.D;
	while ((b) > 1) {
		e_d*=.5D;
		cout<<"e_d: "<<e_d<<endl;
		m_d++;
		b=.5D*e_d + 1.D;


	}

    cout<<endl<<"Liczba bitow mantysy dla double: "<<m_d<<endl;
	cout<<"Epsilon typu double: "<<e_d<<endl;







    return 0;
}

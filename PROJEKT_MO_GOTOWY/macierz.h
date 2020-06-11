#include <iostream>
#include <math.h>
#include "wektor.h"
class macierz{
public:
    void allokuj();
    macierz();
    ~macierz();
    macierz(int nn,int mm);
    macierz(int nn);
    macierz(const macierz &m1);
    macierz( long double **t , int nn,int mm);
    long double operator()(unsigned row, unsigned col);
    void wypisz();
    void wypisz() const;
    void wstaw(wektor w, int row);
    wektor operator*(wektor w1) const;
    macierz operator*(macierz m1)const;
    macierz operator*(long double x)const;
    macierz operator=(macierz ob);
    macierz operator+(macierz ob) const;
    macierz odwroc();
    void rozkladLU(macierz &l,macierz &u) const;
    wektor algThomasa(wektor b);
    void zapisz(string nazwa);
    wektor metodaGS(wektor b)const;
};



/*
wektor wektor::operator*(macierz m1) const {
		if(m1.n != n)
		{
			cout<<"rozne wymiary macierzy i wektora,zwracam pusty wektor!"<<endl;
			return wektor(n);
		}
		wektor temp(n);
		for(int j=0; j<n; j++)
				temp[j] = 0.0;

		for(int i=0 ; i< n ; i++)
		{
			for(int j=0 ; j<m1.m ; j++)
				temp[i] += tab[i]*m1(i,j);

		}

		return temp;
	}*/



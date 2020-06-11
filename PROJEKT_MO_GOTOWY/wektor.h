#pragma once

#include <iostream>
#include <iomanip>
#include <cmath>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <algorithm>
using namespace std;
class macierz;

class wektor{
	public: 
	 int n;
	long double *tab;
	wektor();
	//~wektor();
	wektor(int nn);
	wektor( long double *w,int nn);
	wektor(const wektor &w);
	wektor operator*(macierz m1)const;	
	wektor operator*(long double x) const;
	wektor & operator=(const wektor & ob1);
	wektor operator+(wektor ob) const;
	void wypisz() const;
	long double & operator[](int x);
	long double  operator[](int x) const;
	friend ostream & operator << (ostream &strm, const wektor &ob);
	
};
wektor algThomasa(wektor &b,wektor low,wektor diag,wektor up);

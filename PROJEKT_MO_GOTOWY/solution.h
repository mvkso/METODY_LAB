#pragma once
#include "macierz.h"
#include "calerf.h"
class solution
{
public:

	solution(int ile_xx);
	~solution(void){};

	macierz rozwiazanieT;
	macierz rozwiazanieGS;
	macierz rozwiazanieA;

	int ile_t;
	int ile_x;

	double alfa;
	double beta;

	double lambda;
	double h;
	double dt;

	double D;
	double x_min;
	double x_max;

	double t_max;

	macierz rozwiaz_laasonen_thomasa();
	macierz rozwiaz_laasonen_GS();
	macierz rozwiaz_analityczne();
	void warunekA();
	void warunekT();
	void warunekGS();

	macierz blad_bezwgl_T();
	macierz blad_bezwgl_GS();

	double  blad_max_T;
	double  blad_max_GS;

	macierz blad_T;
	macierz blad_GS;

	
	void saveT(string nazwa);
	void saveGS(string nazwa);
	void saveA(string nazwa);

	void saveT_gnuplot(string nazwa);
	void saveGS_gnuplot(string nazwa);
	void saveA_gnuplot(string nazwa);
};


#include <iostream>
#include <cmath>
#include <iomanip>
#include<fstream>

using namespace std;

double analityczne(double t){
	return 1.0 - exp(-10.0 * (t + atan(t)));
}

double metoda_be(double h, double t_Max) {
	double y = 0.0;
	for (double i = 0.0; i < t_Max; i += h) {
		y += h * (-((10.0 * i * i + 20.0) / (i * i + 1.0)) * (y - 1.0));
	}
	return y;
}

double metoda_pe(double h, double t_Max) {
	double y = 0.0, ulamek = 0.0;
	for (double i = 0.0; i < t_Max; i += h) {
		ulamek = (10.0 * (i + h) * (i + h) + 20.0) / ((i + h) * (i + h) + 1.0);
		y = (y + h * ulamek) / (1.0 + h * ulamek);
	}
	return y;
}

double metoda_t(double h, double t_Max) {
	double y = 0.0, ulamek = 0.0, ulamek_plus = 0.0;
	for (double i = 0.0; i < t_Max; i += h) {
		ulamek = (10.0 * i * i + 20.0) / (i * i + 1.0);
		ulamek_plus = (10.0 * (i + h) * (i + h) + 20.0) / ((i + h) * (i + h) + 1.0);
		y = ((-h / 2.0) * (ulamek * (y - 1.0) - ulamek_plus) + y) / (1.0 + (h / 2.0) * ulamek_plus);
	}
	return y;
}

double blad_mbe(double h, int N)  {
	double blad = 0.0, t = h, y = 0.0, wartosc_dokladna;
	blad = fabs(analityczne(t) - metoda_be(h, t));
	for (int i = 0; i < N; i++, t += h) {
		wartosc_dokladna = analityczne(t);
		y = metoda_be(h, t);
		wartosc_dokladna = fabs(wartosc_dokladna - y);

		if (wartosc_dokladna > blad) {
			blad = wartosc_dokladna;
		}
	}
	return blad;
}

double blad_mpe(double h, int N)  {
	double blad = 0.0, t = h, y = 0.0, wartosc_dokladna;
	blad = fabs(analityczne(t) - metoda_pe(h, t));

	for (int i = 0; i < N; i++, t += h) {
		wartosc_dokladna = analityczne(t);
		y = metoda_pe(h, t);
		wartosc_dokladna = fabs(wartosc_dokladna - y);

		if (wartosc_dokladna > blad) {
			blad = wartosc_dokladna;
		}
	}
	return blad;
}

double blad_mt(double h, int N) {
	double blad = 0.0, t = h, y = 0.0, wartosc_dokladna;
	blad = fabs(analityczne(t) - metoda_t(h, t));

	for (int i = 0; i < N; i++, t += h) {
		wartosc_dokladna = analityczne(t);
		y = metoda_t(h, t);
		wartosc_dokladna = fabs(wartosc_dokladna - y);

		if (wartosc_dokladna > blad) {
			blad = wartosc_dokladna;
		}
	}

	return blad;
}



int main(){
	double mbe_stabilny = 0.0, mbe_niestabilny = 0.0, mpe = 0.0, mt = 0.0;
	const int N = 1000;
	fstream bledy, stabilne, niestabilne;

	bledy.open("bledy.txt", ios::out);


	if (bledy.good()) {
		for (double krok= 20; krok > 1e-20; krok = krok / 1.2) {
			mbe_stabilny = log10(blad_mbe(krok, N));
			mpe = log10(blad_mpe(krok, N));
			mt = log10(blad_mt(krok, N));

			bledy << log10(krok) << "\t" << mbe_stabilny << "\t" << mpe << "\t" << mt << endl;
		}
		bledy.close();
	}
	stabilne.open("dane_stabilne.txt", ios::out);


	if (stabilne.good()) {
		for (double t = 0.0, krok = 0.01; t < 3.0; t += krok) {
			mbe_stabilny = metoda_be(krok, t);
			mpe = metoda_pe(krok, t);
			mt = metoda_t(krok, t);

			stabilne << t << "\t" << analityczne(t) << "\t" << mbe_stabilny << "\t" << mpe << "\t" << mt << endl;
		}
		stabilne.close();
	}
	niestabilne.open("dane_niestabilne.txt", ios::out);


	if (niestabilne.good()) {
		//stabilne dla kroku mniejszego od 0.1,
		//niestabilne dla kroku wieksze od 0.2,
		//posrednie dla kroku [0.1;0.2] (na poczatku niestabilne, ale zbiega sie w czasie)
		for (double t = 0.0, krok = 0.25; t < 5.0; t += krok) {
			mbe_niestabilny = metoda_be(krok, t);
			niestabilne << t << "\t" << analityczne(t) << "\t" << mbe_niestabilny << endl;
		}
		niestabilne.close();
	}
	return 0;
}


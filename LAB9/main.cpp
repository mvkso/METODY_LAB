#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

double funkcja(double x){
	return tan(x);
}
double rozwiazanie_analityczne(double x);

double maksymalny_blad(double *blad, const int rozmiar);


double trzypunktowa_dyskretyzacja_konwencjonalna(double h, const int max);

double dyskretyzacja_numerowa_jednorodna(double h, const int max);

//algorytm Thomasa
void eliminacja(double *u, double *d, double *l, const int rozmiar);

void operacja_rx(double *u, double *d, double *l, double *x, double *b, const int rozmiar);

const double p = 1.0, q = 0.0, r = 4.0;
const double Alfa = 0, Beta = 1.0, Gamma = 0;
// alfa*y'(a) + beta*y(a) + gama = 0 warunek brzegowy mieszany
const double Fi = 0.0, Psi = 1.0, Theta = -0.5;
// fi*y'(a) + psi*y(a) + teta = 0 warunek brzegowy mieszany
const double xPoczatkowe = 0, xKoncowe = M_PI_4;


void load_matrix_from_file(double matrix[][w], const int w, const int k, const char filename[]){
	ifstream file(filename, ios::in);
	if(file.is_open()) {
		for(int i = 0; i < w; i++) {
			for(int j = 0; j < k; j++) {
				file >> matrix[i][j];
			}
		}
		file.close();
	}
}


void load_vector_from_file(double vec[], const int w, const char filename[]) {
	ifstream file(filename, ios::in);
	if(file.is_open()) {
		for(int i = 0; i < w; i++) {
			file >> vec[i];
		}
		file.close();
	}
}

void string_to_file(string line, const char filename[]) {
	ofstream file(filename, ios::out);
	if(file.is_open()) {
		file << line;
		file.close();
	}
}

void append_string_to_file(string line, const char filename[]) {
	ofstream file(filename, ios::out | ios::app);
	if(file.is_open()) {
		file << line;
		file.close();
	}
}

void clear_file(const char filename[]) {
	ofstream file(filename, ios::out);
	if(file.is_open()) {
		file << "";
		file.close();
	}
}





double vector_allocate_number(const int w, int number) {
	double vector[w];
	for(int i = 0; i < w; i++) {
		vector[i] = number;
	}
	return vector;
}




void matrix_print(double matrix[][w], const int w, const int k) {
	for(int i = 0; i < w; i++) {
		for(int j = 0; j < k; j++) {
			cout << matrix[i][j] << " ";
		}
		cout << endl;
	}
}


string vector_to_string(double vec[], const int w) {
	stringstream stream;
	stream << "[ ";
	for(int i = 0; i < w; i++) {
		stream << vec[i] << " ";
	}
	stream << "]";
	return stream.str();
}

void vector_print(double vec[], const int w) {
	for(int i = 0; i < w; i++) {
		cout << vec[i] << " ";
	}
}


double vector_subtract(double vector_1[],double vector_2[], const int w) {
	double vec_buf[w];
	for(int i = 0; i < w; i++) {
		vec_buf[i] = vector_1[i] - vector_2[i];
	}
	return vec_buf;
}

double norm(double vec[])
{
    double next;
    double temp=fabs(vec[0]);
    for (int i=0 ; i<N-1 ; i++)
    {
        next=fabs(vec[i+1]);
        if(fabs(vec[i])<next) temp=next;
    }
    return temp;
}

template< typename T >
void vector_clone(T* vector_1, T* vector_2, const int w) {
	for(int i = 0; i < w; i++) {
		vector_2[i] = vector_1[i];
	}
}

template< typename T >
void matrix_clone(T** matrix_1, T** matrix_2, const int w, const int k) {
	for(int i = 0; i < w; i++) {
		for(int j = 0; j < k; j++) {
			matrix_2[i][j] = matrix_1[i][j];
		}
	}
}


double vector_residuum(double matrix[][w], double x[],double b[], const int w) {
	double sum;
	double vector_res[w];

	for(int i = 0; i < w; i++, sum = 0) {
		for(int j = 0; j < w; j++) {
			sum += matrix[i][j] * x[j];
		}
		vector_res[i] = sum - b[i];
	}
	return vector_res;
}


double estymator(double vector_1[],double vector_2[], const int w) {
	double vec_buf[w];
	vec_buf = vector_subtract(vector_1, vector_2, w);
	double error = norm(vec_buf, w);
	return error;
}


double residuum(double matrix[][w], double x[], double b[], const int w) {
	double vec_residuum[w] = vector_residuum(matrix, x, b, w);
	return vector_norm(vec_residuum, w);
}

int main(){
	int max;
	double h;
	double bladK, bladN;

	fstream bledy;
	bledy.open("Bledy.txt", ios::out);
	for(max = 10; max < 100000; max += 40) {
		//rowny przedzial
		cout << "\r" << max;
		h = (xKoncowe - xPoczatkowe) / (max - 1);
		bladK = log10(trzypunktowa_dyskretyzacja_konwencjonalna(h, max));
		bladN = log10(dyskretyzacja_numerowa_jednorodna(h, max));
		bledy << log10(h) << " " << bladK << " " << bladN << endl;
	}
	bledy.close();
	return 0;
}

double rozwiazanie_analityczne(double x){
	return (1.0 / 4.0) * ((2.0 * x * cos(2.0 * x)) + (2.0 * sin(2.0 * x)) -
	                      (log(2.0) * sin(2.0 * x)) - (2.0 * log(cos(x))) * sin(2.0 * x));
}

double maksymalny_blad(double *blad, const int rozmiar){
	double max = abs(blad[0]);
	for(int i = 0; i < rozmiar; i++) {
		if(fabs(blad[i] > max))
			max = abs(blad[i]);
	}
	return max;
}

double trzypunktowa_dyskretyzacja_konwencjonalna(double h, const int max){
	double u[max];
	double d[max];
	double l[max];
	double b[max];
	double xPrzyblizone[max];
	double bledy[max];
	double punkt = xPoczatkowe;
	double punktPomiarowy = xPoczatkowe;
	double max_blad;

	fstream plik;
	//przypisnie startowych
	d[0] = Beta - Alfa / h;
	u[0] = Alfa / h;
	b[0] = -Gamma;

	for(int i = 1; i < max - 1; i++) {
		//wzory do wypelnienia macierzy thomasa
		l[i - 1] = p / (h * h) - q / (2.0 * h);
		d[i] = (-2.0 * p) / (h * h) + r;
		u[i] = p / (h * h) - q / (2.0 * h);
		b[i] = -funkcja(punkt + h * i);
	}
	l[max - 2] = -Fi / h;
	d[max - 1] = -Fi / h + Psi;
	b[max - 1] = -Theta;

	//algorytm Thomasa
	eliminacja(u, d, l, max);
	operacja_rx(u, d, l, xPrzyblizone, b, max);

	for(int i = 0; i < max; i++) {
		bledy[i] = abs(xPrzyblizone[i] - rozwiazanie_analityczne(punkt));
		punkt += h;
	}

	if(max == 50) {
		plik.open("trzypunktowa_dyskretyzacja_konwencjonalna.txt", ios::out);
		for(int i = 0; i < max; i++) {
			bledy[i] = abs(xPrzyblizone[i] - rozwiazanie_analityczne(punktPomiarowy));
			plik << punktPomiarowy << "\t" << xPrzyblizone[i] << "\t" << rozwiazanie_analityczne(punktPomiarowy) << endl;
			punktPomiarowy += h;
		}
		plik.close();
	}
	max_blad = maksymalny_blad(bledy, max);


	return max_blad;
}

double dyskretyzacja_numerowa_jednorodna(double h, const int max){
	double u[max];
	double d[max];
	double l[max];
	double b[max];
	double xPrzyblizone[max];
	double bledy[max];

	double punkt = xPoczatkowe;
	double punktPomiarowy = xPoczatkowe;
	double max_blad;

	fstream plik;
	//algorytm Thomasa
	//przypisnie startowych
	d[0] = Beta - Alfa / h;
	u[0] = Alfa / h;
	b[0] = -Gamma;

	for(int i = 1; i < max - 1; i++) {
		//wzory do wypelnienia macierzy thomasa
		l[i - 1] = p / (h * h) + r / 12.0;
		d[i] = (-2.0 * p) / (h * h) + r * (10.0 / 12.0);
		u[i] = p / (h * h) + r / 12.0;
		b[i] = -funkcja(punkt + i * h - h) / 12.0 - (10.0 / 12.0) *  funkcja(punkt + i * h) - funkcja(punkt + i * h + h) / 12.0;
	}
	l[max - 2] = -Fi / h;
	d[max - 1] = -Fi / h + Psi;
	b[max - 1] = -Theta;

	//algorytm Thomasa
	eliminacja(u, d, l, max);
	operacja_rx(u, d, l, xPrzyblizone, b, max);

	for(int i = 0; i < max; i++) {
		bledy[i] = abs(xPrzyblizone[i] - rozwiazanie_analityczne(punkt));
		punkt += h;
	}

	if(max == 50) {
		plik.open("dyskretyzacja_numerowa_jednorodna.txt", ios::out);
		for(int i = 0; i < max; i++) {
			bledy[i] = abs(xPrzyblizone[i] - rozwiazanie_analityczne(punktPomiarowy));
			plik << punktPomiarowy << "\t" << xPrzyblizone[i] << "\t" << rozwiazanie_analityczne(punktPomiarowy) << endl;
			punktPomiarowy += h;
		}
		plik.close();
	}
	max_blad = maksymalny_blad(bledy, max);

	return max_blad;
}

//algorytmThomasa
void eliminacja(double *u, double *d, double *l, const int rozmiar){
	for(int i = 1; i < rozmiar; i++) {
		d[i] = d[i] - ((l[i - 1] / d[i - 1]) * u[i - 1]);
	}
}
void operacja_rx(double *u, double *d, double *l, double *x, double*b, const int rozmiar){
	for(int i = 1; i < rozmiar; i++) {
		b[i] = b[i] - ((l[i - 1] / d[i - 1]) * b[i - 1]);
	}

	x[rozmiar - 1] = b[rozmiar - 1] / d[rozmiar - 1];
	for(int i = rozmiar - 2; i >= 0; i--) {
		x[i] = (b[i] - u[i] * x[i + 1]) / d[i];
	}
}

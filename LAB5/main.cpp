#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;
//Macierz

class Macierz{
    int liczba_wierszy,liczba_kolumn;
    double **tablica;

public:
    Macierz(int,int);
    virtual ~Macierz();


    //czy kwadratowa
    bool kwadratowa(){
        return this->liczba_kolumn==this->liczba_wierszy;
        }
    void show()
        {
        for(int i = 0; i < this->liczba_wierszy; i++) {
            for(int j = 0; j < this->liczba_kolumn; j++) {
                cout.precision(3);
                cout.width(11);
                cout.setf(ios::fixed);
                cout << this->tablica[i][j];
            }
            cout << endl;
        }
        cout << endl;
        }

	void showL()
	{
            for(int i = 0; i < this->liczba_wierszy; i++) {
            for(int j = 0; j < i + 1; j++) {
                if(j == i) {
                    cout.precision(3);
                    cout.width(11);
                    cout.setf(ios::fixed);
                    cout << 1.0;
                }
                else{
                    cout.precision(3);
                    cout.width(11);
                    cout.setf(ios::fixed);
                    cout << this->tablica[i][j];
                }
            }
            for(int j = i + 1; j < this->liczba_kolumn; j++) {
                cout.precision(3);
                cout.width(11);
                cout.setf(ios::fixed);
                cout << 0.0;
            }
            cout << endl;
        }
        cout << endl;
	}
void showU()
	{
            for(int i = 0; i < this->liczba_wierszy; i++) {
            for(int j = 0; j < i; j++) {
                cout.precision(3);
                cout.width(11);
                cout.setf(ios::fixed);
                cout << 0.0;
            }
            for(int j = i; j < this->liczba_kolumn; j++) {
                cout.precision(3);
                cout.width(11);
                cout.setf(ios::fixed);
                cout << this->tablica[i][j];
            }
            cout << endl;
        }
        cout << endl;
    }
    double& operator ()(int x, int y){
		return this->tablica[x][y];
	}
	int getW(){
		return this->liczba_wierszy;
	}
	int getK(){
		return this->liczba_kolumn;
	}
     void wypelnij(double macierz[4][4]){

        for(int i = 0; i < this->liczba_wierszy; i++)
		for(int j = 0; j < this->liczba_kolumn; j++)
			this->tablica[i][j] = macierz[i][j];

    }
};
//konstruktor
Macierz::Macierz(int l_w = 4, int l_k = 4){
	if(l_k > 0 and l_w > 0) {
		this->liczba_wierszy = l_w;
		this->liczba_kolumn = l_k;
	}
	else{
		cout << "Blad wymiarow" << endl;
		exit(-1);
	}
	tablica = new double*[this->liczba_wierszy];
	for(int i = 0; i < this->liczba_wierszy; i++)
		tablica[i] = new double[this->liczba_kolumn];
}

//destruktor
Macierz::~Macierz(){
	for(int i = 0; i < this->liczba_wierszy; i++)
		delete[] tablica[i];
	delete[] tablica;
}
///////////////////////////////

class Wektor {
private:
	int rozmiar;
	double *tablica;
public:
	Wektor(int);
	virtual ~Wektor(){
		delete[] tablica;
	};
	void wypelnij(double A[4])
	{
                for(int i = 0; i < this->rozmiar; i++)
            this->tablica[i] = A[i];
	}
	void wyzeruj()
	{
            for(int i = 0; i < this->rozmiar; i++)
            this->tablica[i] = 0.0;
	}
	void show()
	{
          for(int i = 0; i < this->rozmiar; i++) {
            cout.precision(3);
            cout.width(11);
            cout.setf(ios::fixed);
            cout << this->tablica[i];
            cout << endl;
        }
        cout << endl;
        }

	//przeciazenie operatora dla latwiejszego zwracania
	double& operator()(int x){
		return this->tablica[x];
	}
};

Wektor::Wektor(const int rozmiar = 4){
	if(rozmiar > 0)
		this->rozmiar = rozmiar;
	else{
		cout << "B³êdny wymiar wektora!" << endl;
		exit(-1);
	}
	this->tablica = new double[this->rozmiar];
}


int wyborCzesciowy(Macierz& macierz, int j, int n, int *indeksy){
	cout << "Wybor czesciowy elementu podstawowego!" << endl;
	//j=1,n=3
	int wiersz;
	//dla j=1 sprawdzamy wiersz kolejny 2 w danej kolumnie i jezeli jest on wiekszy zamieniamy i tak do konca zjezdzajac w dol macierzy
	for(int i = j + 1; i < n; i++) {
		if(fabs(macierz(indeksy[i], j)) < fabs(macierz(indeksy[i + 1], j))) {
			wiersz = indeksy[i + 1];
		}
		else{
			wiersz = indeksy[i];
		}
	}
	return wiersz;
}



void dekompozycjaLU(Macierz& macierz){
	//sprawdzenie czy jest kwadratowa
	if(macierz.kwadratowa() == false) {
		cout << "Macierz nie kwadratowa" << endl;
		exit(-1);
	}
	int rozmiar = macierz.getW(); //wydobycie rozmiaru(wierszy)
	int* indeksy = new int[rozmiar]; //wektor przechowywujący kolejność wierszy

	for(int i = 0; i < rozmiar; i++)
		indeksy[i] = i; //numeracja wierszy macierzy po kolei

	//k oznacza k etap
	for(int k = 1; k < rozmiar; k++) {
		double element = macierz(indeksy[k - 1], k - 1);

		if(element == 0.0) {
			//bierzemy biezacy (dla nas to bedzie 1 czyli kolejny) wiersz i iterujemy do max rozmiaru
			int wiersz = wyborCzesciowy(macierz, indeksy[k], 3, indeksy);
			indeksy[wiersz] = indeksy[k];
			indeksy[k] = wiersz;
			element = macierz(indeksy[k - 1], k - 1);
		}
		//rozpoczecie iterowania dla tablicy zaczynajacej sie od przekatnej
		for(int i = k; i < rozmiar; i++) {
			//zaczynamy od 'wiersza nizej' bo i=k
			double mnoznik = macierz(indeksy[i], k - 1) / element;
			macierz(indeksy[i], k - 1) = mnoznik;
			for(int j = k; j < rozmiar; j++) {
				macierz(indeksy[i], j) -= macierz(indeksy[k - 1], j) * mnoznik;
				//wyswietlanie informacji po kazdym kroku
				cout << "Krok k(" << k - 1 << ")" << endl;
				cout << "Wiersz i(" << i << ")" << endl;
				cout << "Kolumna j(" << j << ")" << endl;
				macierz.show();
			}
		}
	}
	cout << "Macierz A po dekompozycji:" << endl;
	macierz.show();
	cout << "Macierz L:" << endl;
	macierz.showL();
	cout << "Macierz U:" << endl;
	macierz.showU();
}

void rozwiaz(Macierz& macierz, Wektor& wektor, int rozmiar){
	Wektor x, y;
	cout << "Wektor Pomocniczny Y......................." << endl;
	y.show();
	cout << "Wektor Rozwiazan X......................." << endl;
	x.show();

	//obliczamy wektor Y
	y(0) = wektor(0); //dla 0 0 ze wzoru
	for(int i = 1; i < rozmiar; i++) {
		y(i) = wektor(i);
		for(int j = 0; j < i; j++) {
			y(i) = y(i) - macierz(i, j) * y(j);
			cout << "Wektor Pomocniczny" << endl;
			y.show();
		}
	}

	//obliczamy wektor X
	x(rozmiar - 1) = y(rozmiar - 1) / macierz(rozmiar - 1, rozmiar - 1);
	for(int i = rozmiar - 1; i >= 0; i--) {
		x(i) = y(i);
		for(int j = i + 1; j < rozmiar; j++) {
			//od przekatnej do konca
			x(i) = x(i) - macierz(i, j) * x(j);
		}
		//koncowe podzielenie przez wartosc na przekatnej
		x(i) = x(i) / macierz(i, i);
		cout << "Wektor Rozwiazan" << endl;
		x.show();
	}
}



int main()
{
    const int N = 4;
	double macierz1 [N][N] = {{1.0, -20.0, 30.0, -4.0}, {2.0, -40.0, -6.0, 50.0}, {9.0, -180.0, 11.0, -12.0}, {-16.0, 15.0, 140.0, 13.0}};
	double wektor_w_w [4] = {35.0, 104.0, -366.0, -354.0};

	Macierz A = Macierz(N, N);
	A.wypelnij(macierz1);
	cout << "Macierz A:" << endl;
	A.show();

	Wektor b = Wektor(N);
	b.wypelnij(wektor_w_w);
	cout << "Wektor b:" << endl;
	b.show();

	dekompozycjaLU(A);
	rozwiaz(A, b, N);

    return 0;
}

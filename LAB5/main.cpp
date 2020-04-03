#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;
//Macierz

class Macierz{
    int liczba_wierszy,liczba_kolumn;
    double **tab;

public:
    Macierz(int,int);
    virtual ~Macierz();

    void wypelnij(double macierz[4][4]){

        for(int i = 0; i < this->liczba_wierszy; i++)
		for(int j = 0; j < this->liczba_kolumn; j++)
			this->tablica[i][j] = macierz[i][j];

    }
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





int main()
{
    cout << "Hello world!" << endl;
    return 0;
}

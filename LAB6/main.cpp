#include <iostream>
#include <cmath>

using namespace std;
class Wektor {

int rozmiar;
double* tablica;
public:
Wektor(const int rozmiar = 6){
	if(rozmiar > 0)
		this->rozmiar = rozmiar;
	else{
		cout << "B³êdny wymiar wektora!" << endl;
		exit(-1);
	}
	this->tablica = new double[this->rozmiar];
}

virtual ~Wektor(){
	delete [] tablica;
};
void wypelnij(double* A){
	for(int i = 0; i < this->rozmiar; i++)
		this->tablica[i] = A[i];
}
void wyzeruj(){
	for(int i = 0; i < this->rozmiar; i++)
		this->tablica[i] = 0.0;
}
void show(){
	for(int i = 0; i < this->rozmiar; i++){
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


class Macierz_kwadratowa {
int rozmiar;
Wektor* w1;
Wektor* w2;
Wektor* w3;
public:
Macierz_kwadratowa(const int rozmiar = 6){
        if(rozmiar > 0)
            this->rozmiar = rozmiar;
        else{
            cout << "B³êdny wymiar macierzy!" << endl;
            exit(-1);
        }
        this->w1 = new Wektor(rozmiar);
        this->w2 = new Wektor(rozmiar);
        this->w3 = new Wektor(rozmiar);
    }

virtual ~Macierz_kwadratowa(){
        if(this->w1 != nullptr){
            delete this->w1;
        }
        if(this->w2 != nullptr){
            delete this->w2;
        }
        if(this->w3 != nullptr){
            delete this->w3;
        }
    }


double zwroc(int i, int j){
        if(i == j)
            return (*this->w2)(i);
        else if(i + 1 == j)
            return (*this->w3)(j);
        else if(i - 1 == j)
            return (*this->w1)(j);
        else
            return 0.0;
    }

//zwracanie calych wektorow
    Wektor* zwrocU(){
        return this->w1;
    }
    Wektor* zwrocD(){
        return this->w2;
    }
    Wektor* zwrocL(){
        return this->w3;
    }
};

//functions

void eliminacja(Wektor& L, Wektor& D, Wektor& U, int rozmiar){
	//wzory dla kolejnych od 1 nizej dla obliczenia wspolczynnikow przy eliminacji
	for(int i = 1; i < rozmiar; i++)
		D(i) = D(i) - (U(i - 1) * L(i) / D(i - 1));
}

void operacja_r(Wektor& L, Wektor& D, Wektor& b, int rozmiar){
	//ze wzoru na obliczanie wektora b jako osobne polecenie
	//posiadajac to co mamy obliczamy wyrazy wolne
	for(int i = 1; i < rozmiar; i++)
		b(i) = b(i) - (L(i) * b(i - 1) / D(i - 1));
}

void operacja_x(Wektor& x, Wektor& b, Wektor& U, Wektor& D, int rozmiar){
	//robimy to od konca do poczatku macierzy ze wzoru
	x(rozmiar - 1) = b(rozmiar - 1) / D(rozmiar - 1);
	for(int i = rozmiar - 2; i >= 0; i--)
		x(i) = (b(i) - U(i) * x(i + 1)) / D(i);
}

void algorytmThomasa(Wektor& L, Wektor& D, Wektor& U, Wektor& b, int rozmiar){
	cout << "Algorytm Thomasa" << endl;
	Wektor x(rozmiar);
	x.wyzeruj();

	eliminacja(L, D, U, rozmiar);
	operacja_r(L, D, b, rozmiar);
	operacja_x(x, b, U, D, rozmiar);

	cout << "Wektor Rozwi¹zañ" << endl;
	x.wyswietl();
}



int main()
{

   const int rozmiar = 6;
	Macierz_kwadratowa macierz = Macierz_kwadratowa(rozmiar);
	Wektor wektorWyrazowWolnych = Wektor(rozmiar);

	double b[rozmiar] = {
		31.0,
		165.0 / 4.0,
		917.0 / 30.0,
		851.0 / 28.0,
		3637.0 / 90.0,
		332.0 / 11.0
	};

	double wektorL[rozmiar] = {
		0.0,
		3.0 / 4.0,
		7.0 / 8.0,
		11.0 / 12.0,
		15.0 / 16.0,
		19.0 / 20.0
	};

	double wektorD[rozmiar] = {
		30.0,
		20.0,
		10.0,
		10.0,
		20.0,
		30.0
	};

	double wektorU[rozmiar] = {
		2.0 / 3.0,
		5.0 / 6.0,
		9.0 / 10.0,
		13.0 / 14.0,
		17.0 / 18.0,
		0.0
	};

	wektorWyrazowWolnych.wypelnijZadaniem(b);
	cout << "Wektor wolnych wyrazow" << endl;
	wektorWyrazowWolnych.wyswietl();

	macierz.zwrocL()->wypelnijZadaniem(wektorL);
	cout << "Wektor L" << endl;
	macierz.zwrocL()->wyswietl();

	macierz.zwrocD()->wypelnijZadaniem(wektorD);
	cout << "Wektor D" << endl;
	macierz.zwrocD()->wyswietl();

	macierz.zwrocU()->wypelnijZadaniem(wektorU);
	cout << "Wektor U" << endl;
	macierz.zwrocU()->wyswietl();

	algorytmThomasa(*(macierz.zwrocL()), *(macierz.zwrocD()), *(macierz.zwrocU()), wektorWyrazowWolnych, rozmiar);

    return 0;
}

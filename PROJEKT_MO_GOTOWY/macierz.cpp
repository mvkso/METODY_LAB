

//#include "StdAfx.h"



#include "macierz.h"

void macierz::allokuj(){

	try{
		tab = new  long double*[n];
			for(int i=0;i<n;i++)
					tab[i]= new  long double[m];
	}
	catch(bad_alloc a){
		cout<<"blad alokacji pamieci!"<<endl;
		system("pause");
	}
	//czyszczenie:
	for(int i=0; i<n; i++)
		for(int j=0; j<m; j++)
			tab[i][j]=0.0;

}
macierz::macierz(){
		n=2;
		m=2;
		allokuj();

	}
macierz::~macierz()
{
	for(int i=0;i<n;i++)
		delete[] tab[i];
	delete[] tab;
}
macierz::macierz(int nn,int mm){
		if(nn<0 || mm<0) {cout<<"macierz::macierz(int,int) zly rozmiar macierzy!"<<endl;}
			n=abs(nn);
			m=abs(mm);
			allokuj();

		}
macierz::macierz(int nn){
		if(nn<0) {cout<<"macierz:: zly rozmiar macierzy!"<<endl; }
			n=abs(nn);
			m=n;
			allokuj();

		}
macierz::macierz(const macierz &m1){
		n=m1.n;
		m=m1.m;
		allokuj();

			for(int i=0; i<n; i++)
				for(int j=0; j<m; j++)
					tab[i][j]=m1(i,j);

	}
macierz::macierz( long double **t , int nn,int mm){
		if(nn<0 || mm<0) {cout<<"macierz::macierz(double**,int,int) zly rozmiar macierzy!"<<endl;}
			n=abs(nn);
			m=abs(mm);
			allokuj();

		for(int i=0 ; i<n ;i++)
			for(int j=0 ; j<m ;j++)
				tab[i][j]=t[i][j];
	}
 long double &macierz::operator() (unsigned row, unsigned col){

	return tab[row][col];
}
 long double  macierz::operator() (unsigned row, unsigned col) const {

	return tab[row][col];
}

 void macierz::wstaw(wektor w,int row )
 {
	 if(m==w.n){
		 for(int i=0;i<m;i++){
			 tab[row][i] = w[i];
		 }
	 }
	 else{
		 cout<<"rozne wymiary macierza i wektora, nie wstawiono wektora do macierzy !! "<<endl;
		 system("pause");
	 }
 }

 wektor macierz::operator*(wektor w1) const{
		if(w1.n != n)
		{
			cout<<"wektor macierz::operator*(wektor w1): rozne wymiary macierzy i wektora"<<endl;
			return wektor(n);
		}

		wektor temp(n);
		for(int j=0; j<n; j++)
				temp[j] = 0.0;

		for(int i=0 ; i< n ; i++)
		{
			for(int j=0 ; j<m ; j++)
				temp[i] += w1[i]* tab[i][j];

		}

		return temp;
	}
void macierz::wypisz() const {
		for(int i=0; i<n; i++){

			for(int j=0; j<m; j++)
			{
				cout.width(7);
				cout<<setprecision(3)<<tab[i][j];}
				cout<<endl;
			}
	}
macierz macierz::operator*(macierz m1) const{
		if( m!= m1.n){
			cout<<"macierz operator*: zle wymiary macierzy! zwracam pusta macierz!"<<endl;
			return macierz(n,m);
		}
		int mm=n,pp=m,nn=m1.m;
		macierz temp(mm,nn);

		 long double s;
		for(int i = 0; i < mm ; i++){
			for(int j = 0; j < nn ; j++)
			{

				s = 0.0;
				for(int k = 0; k < pp ; k++)
					s += tab[i][k] * m1(k,j);

				temp(i,j) = s;
			}
		}
		return temp;
	}
macierz macierz::operator*( long double x) const {
		macierz temp(n,m);

		for(int i=0;i<n;i++)
			for(int j=0 ; j<m ; j++)
				temp(i,j)=tab[i][j]*x;

		return temp;
	}
macierz & macierz::operator=(macierz ob){
		if(ob.tab)
		{
			if(tab) delete[] tab;
			  tab=NULL;


			n=ob.n;
			m=ob.m;

			allokuj();

			for(int i=0 ; i<n ; i++)
				for(int j=0 ; j<m ; j++)
					tab[i][j] = ob(i,j);


			}
			else{
				 if(tab)
					delete[] tab;
				 tab=NULL;
				}

		return *this;
	}
macierz macierz::operator +(macierz m1) const {
		if(n != m1.n || m!= m1.m){
			cout<<"zle wymiary macierzy! zwracam pusta macierz!"<<endl;
			return macierz(n,m);
		}
		macierz temp(n,n);



		for(int i=0; i<n; i++)
			for(int j=0; j<n; j++)
				temp(i,j) = tab[i][j] + m1(i,j);

		return temp;
	}
macierz macierz::odwroc(){
	if(n!=m){cout<<" n != m =>nie  da sie odwrocic!"; return *this;	}

	macierz I(n,n);
	macierz A(*this);


	for(int i=0 ; i<n ; i++)
		for(int j=0 ; j<n ; j++)
		{													// Tworzenie macierzy jednostkowej
			if(i==j)										// | 1  0  0 |
				I(i,j)=( long double) 1.0;					// | 0  1  0 |
			else											// | 0  0  1 |
				I(i,j)=( long double) 0.0;


		}


	int k,j,i;
	long double temp=0.0;

	for(k=0;k<n;k++)
	{

		if(A(k,k) == 0.0){

			//JEZELI NA PRZEKATNEJ JEST 0.0:
			//PROCEDURA WYBIERANIA CZESCIOWEGO ELEMENTU PODSTAWOWEGO:
			int max = k+1;
			for(int mm=max+1 ; mm<n; mm++){
				if( fabs(A(mm,k)) > fabs(A(max,k)))
					max=mm;

				}
			if (max < n && k<n && A(max,k)!=0.0)
			{
				cout<<" nie da sie odwrocic\n";//system("pause");
				getchar();
				cout<<" Uwaga! bledna macierz :"<<endl;

				return macierz(n);
			}

			//zamiana wierszy w macierzy A i I  o indeksie 'max'
			long double *kopia;
			kopia = new long  double[n];

			for(int z=0; z<n ;z++)
				kopia[z] = A(k,z);

			for(int z=0; z<n ;z++)
				A(k,z) = A(max,z);

			for(int z=0; z<n ;z++)
				A(max,z)=kopia[z];
			//////////////////////////////////
			for(int z=0; z<n ;z++)
				kopia[z] = I(k,z);

			for(int z=0; z<n ;z++)
				I(k,z) = I(max,z);

			for(int z=0; z<n ;z++)
				I(max,z)=kopia[z];
			delete[] kopia;

		}


			temp=A(k,k);

		for(j=0;j<n;j++)
		{										//DZIELENIE CA£EGO WIERSZA PRZEZ WSPÓ£CZYNNIK ABY UZYSKAC 1 NA PRZEKATNEJ
			A(k,j)/=temp;
			I(k,j)/=temp;
		}
		/////////////////////////////////////
		for(i=0;i<n;i++)						//przemieszczanie wzd³u¿ kolumny, w dó³
		{
			temp=A(i,k);
			for(j=0;j<n;j++)
			{
				if(i==k)						// operacje dzielenia przez wspolczynnik( taki zeby uzyskac 0 pod przekatna) wykonujemy dla czesci wiersza pod przekatna macierzy A i I
					break;						// W ELEMENCIE pod przekatna maja byc 0.0
				A(i,j) -= A(k,j)*temp;
				I(i,j) -= I(k,j)*temp;
			}
		}
	}


	return I;
}
macierz operator *( long double x,macierz m1) {
	macierz temp(m1.n,m1.m);

		for(int i=0;i<m1.n;i++)
			for(int j=0 ; j<m1.m ; j++)
				temp(i,j)=m1(i,j)*x;

		return temp;



}
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
	}
void macierz::rozkladLU(macierz &l,macierz &u) const {
	if(n != m ||n != l.n || l.n != l.m || l.n != u.n || u.n!= u.m){
		cout<<"Rozne wymiary macierzy A,L,U !!!"<<endl;
		return;
	}

	//rozklad  LU :
	for(int k=0; k<n; k++){
		//przyjmujemy l[k][k] =1.0
		l(k,k) =1.0;

		 long double suma=0.0;
		for(int s=0; s < k -1;s++)
			suma += l(k,s)*u(s,k);

		u(k,k) = tab[k][k] - suma;


		for(int j=k; j<n ; j++){

			suma=0.0;
			for(int s=0; s < k ;s++)
				suma += l(k,s)*u(s,j);

			u(k,j) = tab[k][j] - suma;
			u(k,j) /= l(k,k);

			suma=0.0;
			for(int s=0; s < k ;s++)
				suma += l.tab[j][s]*u.tab[s][k];

			l(j,k) = tab[j][k] - suma;
			l(j,k) /= u(k,k);

		}

	}

	for(int i=0 ; i<n ;i++)
		for(int j=0 ;j<n;j++)
		{
			if(i<j)
				l(i,j)=0.0;
			else if(i>j)
				u(i,j)=0.0;

		}


}
ostream & operator<<(ostream &strm, const macierz & a) {

	for(int i=0; i<a.n; i++){
			strm<<"|";
			for(int j=0; j<a.m; j++)
			{
				strm.width(9);
				strm<<setprecision(2)<<a(i,j);
			}
			strm<<"|";
			strm<<endl;
	}
	return strm;
}
wektor macierz::algThomasa(wektor b){

		wektor x(n);
		wektor a(n-1),d(n),c(n-1);

		//przepisanie wartosci "trzech przekatnych" macierzy do oddzielnych wektorow
		int i;
		for( i=0 ; i<n-1; i++){
			a[i] = tab[i+1][i];
			d[i] = tab[i][i];
			c[i] = tab[i][i+1];

		}
		d[n-1]=tab[n-1][n-1];
		//cout<<"a: \n"<<a<<"d: \n"<<d<<"c: \n"<<c<<endl;

		//czy diagonalnie dominujaca:
		for(int i=1; i< n;i++)
			if(fabs(d[i]) <= fabs(a[i])+fabs(c[i])){
				cout<<" Uwaga! macierz nie jest diagonalnie dominujaca! brak rozwiazania!"<<endl;
				return b;
			}

		long double mnoznik;
		for(i=1;i<n;i++){


			if(d[i-1] == 0.0){
				cout<<"procedura czesciowego wyboru elementu podstawowego"<<endl;
				int k=i ;
				int max=k;

				for(k++; k<n-1 ;k++)
					if(fabs(d[k]) > fabs(d[max]))
						max=k;

				//cout<<"max: "<<max<<endl;
				 long double temp;
				temp=d[i-1];
				d[i-1]=d[max];
				d[max]=temp;

				//cout<<"d: "<<d;

				temp=a[i-1];
				a[i-1]=a[max];
				a[max]=temp;
				//cout<<"a: "<<a;

				temp=b[i-1];
				b[i-1]=b[max];
				b[max]=temp;
				//cout<<"b: "<<b;

				temp=c[i-1];
				c[i-1]=c[max];
				c[max]=temp;
				//cout<<"c: "<<c<<endl;

			}

			mnoznik = (a[i-1]/d[i-1]);
			d[i] -= mnoznik * c[i-1];
			b[i] -= mnoznik * b[i-1];
		}
		x[n-1]=b[n-1] / d[n-1];

		for(i=n-2;i>-1;i--){
			x[i] = ( b[i] - c[i]*x[i+1])/d[i]; //obliczanie wartosci wektoraX od n-2 do x[0], od do³u0 od góry
		}



	return x;

}

void macierz::zapisz( string nazwa )
{

	fstream fl;
	fl.open(nazwa.c_str(), ios::out);
	if (!fl.is_open())
	{
		cout<<"Blad otwierania pliku! exit(-1)"<<endl;
		system("pause");
		exit(-1);
	}

	for (int i=0;i<n;i++)
	{

		//fl<<setw(15)<<l2s(t).c_str();

		for (int j=0 ; j<m;j++)
		{
			fl<<setw(15)<<l2s(tab[i][j]).c_str();
		}

		fl<<endl;
	}



	fl.close();
	cout<<"Zapisano macierz do pliku :  \""<<nazwa.c_str()<<"\""<<endl;

}

wektor macierz::metodaGS( wektor b ) const
{
	int Iteracje=50; //ilosc iteracji algorytmu
	double eps=1e-06;//dokladnosc rozwiazania
	macierz A=*this;
	if(A.n != b.n){ cout<<"metoda Gaussa - Seidla: rozne wymiary macierzy i wektora! brak rozwiazan!"<<endl; return wektor(b.n);}

	int n=A.n;
	wektor xn(n),x_n(n); //wektory rozwiazan

	//wartosci poczatkowe wektora x(n-1):
	for(int i=0 ; i<n ; i++){
		x_n[i]=0.0;
		xn[i] =0.0;
	}
	for( int k=0 ; k<Iteracje ; k++){
		for (int i=0; i<n;i++){

			double sum1 = 0;
			double sum2 = 0;
			for (int j = 0 ; j < i ; j++)
				sum1 += A(i,j) * xn[j];

			for (int j = i + 1 ; j < n ; j++)
				sum2 += A(i,j) * x_n[j];

			xn[i] = (b[i] - sum1 - sum2) / A(i,i);

		}

		double max_diff = abs(xn[0] - x_n[0]);
		for (int i = 1 ; i < n ; i++){
			double diff = abs(xn[0] - x_n[0]);
			if (diff > max_diff)
				max_diff = diff;
		}

		double max_x = x_n[0];
		for (int i = 1 ; i < n ; i++){
			if (x_n[i] > max_x)
				max_x = x_n[i];
		}

		if (max_diff < max_x * eps){
			return xn;
		}
		x_n = xn;
	}



	return xn;
}


string l2s(double liczba){

	using std::replace;
	char str[50];
	sprintf_s(str,"%lf",(double)liczba,sizeof(double));
	string s(str);

	replace(s.begin(), s.end(), '.', ',');
	return s;

}

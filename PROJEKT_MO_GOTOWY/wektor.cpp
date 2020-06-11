
//#include "StdAfx.h"

#include "wektor.h"


//wektor::~wektor(){	delete[] tab;}
void wektor::wypisz() const {
		
		for(int i=0 ; i<n ; i++)
		{cout<<" ";cout.width(2);cout<<tab[i];}
		cout<<endl;
		
}
long double & wektor::operator[](int x){
		return tab[x];
	}
long  double  wektor::operator[](int x)const{
		return tab[x];
	}
wektor wektor::operator+(wektor ob) const {
		
		if(n != ob.n) { cout <<"wektor wektor::operator+(wektor ob): rozne wymiary wektorow!"<<endl; return wektor(n);}
		
		wektor temp(n);
		
		for(int i=0 ; i<n ; i++)
			temp[i]=tab[i]+ob[i];
		
		
		return temp;
		
}
wektor::wektor(){
		n=2;
		tab = new long  double[2];
	}
wektor::wektor(int nn){
		if(nn<0) { cout<<" zly rozmiar wektora!"<<endl; tab=NULL; return;}
		n=nn;
		tab=new long  double[n];
	}
wektor::wektor( long double *w,int nn){
		if(nn<0) {cout<<" zly rozmiar wektora!"<<endl; tab=NULL; return;}
		n=nn;
		tab=new  long double[n];
		for(int i=0 ; i<n ; i++)
			tab[i]=w[i];
	}
wektor::wektor(const wektor &w){
		n=w.n;
		tab=new long  double[n];
		for(int i=0 ; i<n ; i++)
			tab[i]=w[i];
	}
wektor wektor::operator*( long double x) const {
		
		wektor temp(n);
		
				
		for(int i=0 ; i< n ; i++)
		{
			
				tab[i]= tab[i]*x;
			
		}
		
		return temp;
	}
wektor & wektor::operator=(const wektor & ob1){
		if(ob1.tab)
		{
		if(tab) delete[] tab;
		 tab=NULL;
	     n=ob1.n;
		 tab=new  long double[n];

		for(int i=0 ; i<n ; i++)
			tab[i] = ob1.tab[i];


	}
	else{
		 if(tab)
			 delete []tab;
		 tab=NULL;

	}

	return *this;
		
}

wektor algThomasa( wektor &b,wektor low,wektor diag,wektor up )
{
	int n = b.n;
	wektor x(n);		
	
	//cout<<a<<d<<c<<endl;
	//przepisanie wartosci "trzech przekatnych" macierzy do oddzielnych wektorow
	int i;
	
	//czy diagonalnie dominujaca:	
	for(int ii=1; ii< n;ii++)
		if(fabs(diag[ii]) <= fabs(up[ii])+fabs(low[ii])){				
			cout<<" Uwaga! macierz nie jest diagonalnie dominujaca! brak rozwiazania!"<<endl;
			return b;
		}
	double mnoznik;

	for(i=1;i<n;i++){

			
			/*if(diag[i-1] == 0.0){
			cout<<"procedura czesciowego wyboru elementu podstawowego"<<endl;
			int k=i ;
			int max=k;

			for(k++; k<n-1 ;k++)
			if(fabs(diag[k]) > fabs(diag[max]))
			max=k;

			//cout<<"max: "<<max<<endl;
			long double temp;
			temp=diag[i-1];
			diag[i-1]=d[max];
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
			/**/

		//mnoznik = (up[i-1]/diag[i-1]);
		//diag[i] -= mnoznik * low[i-1];
		//b[i] -= mnoznik * b[i-1];
			b[i] = b[i] - (low[i-1]/diag[i-1])*b[i-1];
	}
	x[n-1]=b[n-1] / diag[n-1];		

	for( i=n-2;i>-1;i--){
		x[i] = ( b[i] - up[i]*x[i+1]  )/diag[i]; 
	}



		return x;

}


ostream & operator << (ostream &strm, const wektor &ob){
	strm<<"[";	
	for(int i=0 ; i<ob.n ; i++){
		strm.width(6);
		strm<<setprecision(3)<<ob.tab[i];
		if(i!=ob.n-1)
			strm<<",";
	}
		strm<<"]"<<endl;

	return strm;
}

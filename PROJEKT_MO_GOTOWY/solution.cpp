//#include "StdAfx.h"
#include "solution.h"



solution::solution(int ile_xx)
{
	D=1.0;
	t_max=2.0;
	alfa=1.0;
	beta=0.0;
	// 1.
	x_min =-6.*sqrt(D*t_max)-1.0;
	x_max = -x_min; 
	
	
	ile_x  = ile_xx;
	
	// 2. obliczenie kroku h i dt
	h  =  (x_max -x_min)/ (ile_x-1);  
	
	dt = h * h;
	//dt = t_max/(ile_t -1); 


	ile_t =(int)( t_max/dt +1);
	
	dt = (t_max/(ile_t-1));

	lambda = D * dt / (h*h);

	cout<<"ile_x = "<<ile_x<<"\tile_t = "<<ile_t<<"\nt_max = "<<t_max<<endl;
	cout<<"h = "<<h<<"\tkrok dt = "<<dt<<endl;
	cout<<"lambda = "<<lambda<<endl;
	

	rozwiazanieT  = macierz(ile_t,ile_x);
	rozwiazanieA  = macierz(ile_t,ile_x);
	rozwiazanieGS = macierz(ile_t,ile_x);

	blad_GS=macierz(ile_t,ile_x);
	blad_T=macierz(ile_t,ile_x);
}

macierz solution::rozwiaz_analityczne()
{
	double x,t;

	warunekA();

	//cout<<rozwiazanieA;
	
	
	x=x_min+h;
	t=dt;

	for (int i=1;i<ile_t;i++)
	{
		for (int j=1 ; j<ile_x-1;j++)
		{

			rozwiazanieA(i,j) = 0.5 * erfc( x / (2.0 * sqrt(D*t)));

			x+=h;
		}
		x=x_min+h;
		t+=dt;


	}
	return rozwiazanieA;
}

macierz solution::rozwiaz_laasonen_thomasa(){
	cout<<"rozwiaz_laasonen_thomasa():"<<endl;
	wektor up(ile_x-1);
	wektor diag(ile_x);
	wektor low(ile_x-1);

	wektor b(ile_x);
	wektor wyn(ile_x);


	warunekT();

	up[0]  =lambda;
	diag[0]=1.0;
	low[0]=lambda;



	for(int i=1 ; i<ile_x-1; i++){
		up[i]   = lambda;
		diag[i] =-( 1 + (2*lambda) );
		low[i] = lambda;
	}

	diag[ile_x-1]=1.0;


	//pierwszy etap alg Thomassa - zerowanie podprzekatnej:
	for(int i=1; i< ile_x;i++)
		diag[i] = diag[i] - (low[i-1]/diag[i-1])*up[i-1];

	// przebieg czasowy:
	for( int k = 1; k < ile_t; k++ ) {
		
		b[0] = alfa*2; //<===MAGICZNA DWOJKA !!
			
		for(int i=1 ; i<ile_x; i++)		
			b[i] = - rozwiazanieT(k-1,i);	

		b[ile_x-1] = beta ;
		
		

		wyn = algThomasa(b,low,diag,up);

	
		//kopiowanie wynikow z chwili "k'atej" do macierzy wynikow:
		for( int i = 1; i < ile_x-1; i++ )
			rozwiazanieT(k,i) = wyn[i]; 


	}
	return rozwiazanieT;
}

macierz solution::rozwiaz_laasonen_GS(){
	double x;
	wektor b(ile_x);
	wektor wyn(ile_x);
	
	x = x_min;
	for( int i = 0; i < ile_x; i++ ){ //warunek poczatkowy

		if(x<0)			  
			rozwiazanieGS(0,i) = 1.0;
		else
			rozwiazanieGS(0,i) = 0.0;
		x=x+h;
	}


	for( int i = 0; i < ile_t; i++ ) {// warunki brzegowe 
		rozwiazanieGS(i,0) = alfa;
		rozwiazanieGS(i,ile_x-1) = beta;
	}
	
	macierz A(ile_x,ile_x);//macierz trojprzekatniowa
	A(1,0)=lambda;
	A(0,1)=lambda;
	A(0,0)=1.0;		
	for( int i = 1; i < ile_x-1; i++ )//wypelnianie macierzy 
	{
		A(i+1,i) = lambda;//down
		A(i,i) = -( 1 + (2*lambda) );//diag
		A(i,i+1) = lambda;//up		
	}
	A(ile_x-1,ile_x-1) = 1.0
	
	//przebieg czasowy
	for( int k = 1; k < ile_t; k++ )
	{		
		
	
		b[0] =  alfa ;
		for( int i = 1; i < ile_x-1; i++ )//wypelnianie macierzy i 
			b[i] = -rozwiazanieGS(k-1,i);
		b[ile_x-1] = beta ;

		wyn = A.metodaGS(b);
		
		for( int i = 1; i < ile_x-1; i++ )//kopiowanie wynikow 
			rozwiazanieGS(k,i) = wyn[i]; 
	}	
	return rozwiazanieGS;
}



void solution::saveT( string nazwa )
{
	
	fstream fl;
	fl.open(nazwa.c_str(), ios::out);
	if (!fl.is_open())
	{
		cout<<"Blad otwierania pliku!"<<endl;
		system("pause");
		return;
	}
	//fl<<"Alg Thomassa:"<<endl;


	//double t,x=x_min;
	//w formie macierzy:

/*
	fl<<setw(15)<<"\t"<<setw(15)<<l2s(x).c_str();
	x+=h;
	for (int j=1 ; j<ile_x;j++,x+=h)
	{

		fl<<setw(15)<<l2s(x).c_str();
	}
	*/
	//fl<<endl<<setw(9)<<"t:"<<endl;
	//t=0.0;
	//fl<<endl;
	
	
	for (int i=0;i<ile_t;i++)
	{

		//fl<<setw(15)<<l2s(t).c_str();

		for (int j=0 ; j<ile_x;j++)
		{

			fl<<setw(15)<<l2s(rozwiazanieT(i,j)).c_str();
		}

		fl<<endl;
	}

	

	fl.close();
	cout<<"Zapisano do pliku :  \""<<nazwa.c_str()<<"\""<<endl;

}
void solution::saveGS( string nazwa )
{

	fstream fl;
	fl.open(nazwa.c_str(), ios::out);
	if (!fl.is_open())
	{
		cout<<"Blad otwierania pliku!"<<endl;
		system("pause");
		return;
	}
//	fl<<"Alg Gaussa-Seidla:"<<endl;


	//double t,x=x_min;
	//w formie macierzy:


	/*fl<<setw(15)<<"\t"<<setw(15)<<l2s(x).c_str();
	x+=h;
	for (int j=1 ; j<ile_x;j++,x+=h)
	{

		fl<<setw(15)<<l2s(x).c_str();
	}
	//fl<<endl<<setw(9)<<"t:"<<endl;
	*/
	//t=0.0;
	//fl<<endl;
	for (int i=0;i<ile_t;i++)
	{

		//fl<<setw(15)<<l2s(t).c_str();

		for (int j=0 ; j<ile_x;j++)
		{

			fl<<setw(15)<<l2s(rozwiazanieGS(i,j)).c_str();
		}

		fl<<endl;
	}



	fl.close();
	cout<<"Zapisano do pliku :  \""<<nazwa.c_str()<<"\""<<endl;
}
void solution::saveA( string nazwa )
{
	fstream fl;
	fl.open(nazwa.c_str(), ios::out);
	if (!fl.is_open())
	{
		cout<<"Blad otwierania pliku!"<<endl;
		system("pause");
		return;
	}
	//double t,x=x_min;
	//fl<<"Rozw Analityczne:"<<endl;


	
	//w formie macierzy:

/*
	fl<<setw(15)<<"\t"<<setw(15)<<l2s(x).c_str();
	x+=h;
	for (int j=1 ; j<ile_x;j++,x+=h)
	{

		fl<<setw(15)<<l2s(x).c_str();
	}
	*/
	//fl<<endl<<setw(9)<<"t:"<<endl;
	//t=0.0;
	//fl<<endl;
	for (int i=0;i<ile_t;i++)
	{

		//fl<<setw(15)<<l2s(t).c_str();

		for (int j=0 ; j<ile_x;j++)
		{

			fl<<setw(15)<<l2s(rozwiazanieA(i,j)).c_str();
		}

		fl<<endl;
	}



	fl.close();
	cout<<"Zapisano do pliku :  \""<<nazwa.c_str()<<"\""<<endl;
}

macierz solution::blad_bezwgl_T()
{
	blad_max_T=0.0;
	//double x = x_min,t = 0.0;

	

	
	for (int j = 0; j < ile_x; j++)
	{
		for (int i = 0; i < ile_t ; i++)
		{
			blad_T(i,j)=fabs(rozwiazanieA(i,j) - rozwiazanieT(i,j));
			
		}
		if(blad_T(ile_t-1,j)>blad_max_T)
				blad_max_T = blad_T(ile_t-1,j); //blad bezwzgledny
		
	}  
	return blad_T;
}
macierz solution::blad_bezwgl_GS()
{
	//double x = x_min,t = 0;

	
	blad_max_GS=0.0;
	
	for (int j = 0; j < ile_x; j++)
	{
		for (int i = 0; i < ile_t ; i++)
		{
			blad_GS(i,j)=fabs(rozwiazanieA(i,j) - rozwiazanieGS(i,j));
			
		}
		if(blad_GS(ile_t-1,j)>blad_max_GS)
				blad_max_GS = blad_GS(ile_t-1,j); //blad bezwzgledny
		
	}  
	return blad_GS;
}

void solution::warunekA()
{
	double x ;
	for( int i = 0; i < ile_x; i++ ){ //warunek poczatkowy

		x=h*(double)i+x_min;

	
		if (x>=0)
			rozwiazanieA(0,i) = 0.0;

		else if   (x<0)
			rozwiazanieA(0,i) = 1.0;

		


	}

	for( int i = 0; i < ile_t; i++ ) {// warunki brzegowe , dla kzdej chwili T - konkretna wartosc 1. i ostatniej zmiennej X
		rozwiazanieA(i,0) = alfa;
		rozwiazanieA(i,ile_x-1) = beta;
	}

}

void solution::warunekGS(){

	//warunek poczatkowe: U(x,0) = 1 if(x<0) else if(x>=0) 0.0 

	//warunek brzegowy :  U(x_min,t)=1  U(x_max,t)=0;

	double x ;
	for( int i = 0; i < ile_x; i++ ){ //warunek poczatkowy

		x=h*(double)i+x_min;


		if (x>=0)
			rozwiazanieGS(0,i) = 0.0;

		else if   (x<0)
			rozwiazanieGS(0,i) = 1.0;

		


	}



	for( int i = 0; i < ile_t; i++ ) {// warunki brzegowe 
		rozwiazanieGS(i,0) = alfa;
		rozwiazanieGS(i,ile_x-1) = beta;
	}


}

void solution::warunekT(){

	//warunki poczatkowe: U(x,0) = 1 if(x<0) else  :=0.0 

	//warunki brzegowy :  U(x_min,t)=alfa,   U(x_max,t)=beta;

	//warunek poczatkowy dla t=0, rozwiazanieT => macierz(ile_t,ile_x);
	double x ;
	for( int i = 0; i < ile_x; i++ ){ //warunek poczatkowy
		
		x=h*(double)i+x_min;
		
		if      (x<0)
			rozwiazanieT(0,i) = 1.0;

		else if (x>=0)
			rozwiazanieT(0,i) = 0.0;

		
		
	}

	for( int i = 0; i < ile_t; i++ ) {// warunki brzegowe 
		rozwiazanieT(i,0) = alfa;
		rozwiazanieT(i,ile_x-1) = beta;
	}
	

}

void solution::saveT_gnuplot( string nazwa )
{

	fstream fl;
	fl.open(nazwa.c_str(), ios::out);
	if (!fl.is_open())
	{
		cout<<"Blad otwierania pliku!"<<endl;
		system("pause");
		return;
	}


	double t,x;
	//w formie macierzy:



	x=x_min;	
	t=0.0;
	for (int i=0;i<ile_t;i++,t=dt*i)
	{

		for(int j=0 ;  j<ile_x ; j++,x=x_min + h*j ){		

			fl<<t<< "\t" <<  x  <<"\t"<<rozwiazanieT(i,j)<<endl;									
		}
		fl<<endl;
		x=x_min;	

	}



	fl.close();
	cout<<"Zapisano do pliku :  \""<<nazwa.c_str()<<"\""<<endl;
}


void solution::saveGS_gnuplot( string nazwa )
{

	

		fstream fl;
		fl.open(nazwa.c_str(), ios::out);
		if (!fl.is_open())
		{
			cout<<"Blad otwierania pliku!"<<endl;
			system("pause");
			return;
		}


		double t,x;
		//w formie macierzy:



		x=x_min;	
		t=0.0;
		for (int i=0;i<ile_t;i++,t=dt*i)
		{

			for(int j=0 ;  j<ile_x ; j++,x=x_min + h*j ){		

				fl<<t<< "\t" <<  x  <<"\t"<<rozwiazanieGS(i,j)<<endl;									
			}
			fl<<endl;
			x=x_min;	

		}



		fl.close();
		cout<<"Zapisano do pliku :  \""<<nazwa.c_str()<<"\""<<endl;
	
}

void solution::saveA_gnuplot( string nazwa )
{

	
		fstream fl;
		fl.open(nazwa.c_str(), ios::out);
		if (!fl.is_open())
		{
			cout<<"Blad otwierania pliku!"<<endl;
			system("pause");
			return;
		}
	

		double t,x;
		//w formie macierzy:


		
		x=x_min;	
		t=0.0;
		for (int i=0;i<ile_t;i++,t=dt*i)
		{
			
			for(int j=0 ;  j<ile_x ; j++,x=x_min + h*j ){		
			
					fl<<t<< "\t" <<  x  <<"\t"<<rozwiazanieA(i,j)<<endl;									
			}
			fl<<endl;
			x=x_min;	
		
		}



	fl.close();
	cout<<"Zapisano do pliku :  \""<<nazwa.c_str()<<"\""<<endl;

}



/*
Uklad::redukujDoThomasa() {
	for(int i=1; i<A.n; i++)  
		A.wsk[i][1] = A.wsk[i][1] - (A.wsk[i-1][0]/A.wsk[i-1][1])*A.wsk[i-1][2] ; 
}
void Uklad::rozwiazThomas() { 
	for(int i=1; i<A.n; i++)  
		b[i] = b[i] - (A.wsk[i-1][0]/A.wsk[i-1][1])* b[i-1];  

	x[A.n-1] = b[A.n-1]/A.wsk[A.n-1][1];
	
	for(int i=A.n-2; i>=0; i--) 
		x[i] = (b[i]-A.wsk[i][2]*x[i+1])/A.wsk[i][1]; //return x; 
}
*/


#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;


void wybor_czesciowy(double macierz[][4],double macierzL[][4], double *wektor, int i)
{
    double tmp[4]={0,0,0,0};
    double y;
    int a;
    double x,maximum=fabs(macierz[i+1][i]);
    for(int j=i+1 ; j<4 ; j++)
    {
        x=fabs(macierz[j][i]);
        if(maximum<x)
        {
            maximum=x;
            a=j;
        }
    }
    for(int k=0; k<4 ; k++)
    {
        tmp[k]=macierz[i][k];
        macierz[i][k]=macierz[a][k];
        macierz[a][k]=tmp[k];
    }
    for(int k=0; k<4 ; k++)
    {
        tmp[k]=macierzL[i][k];
        macierzL[i][k]=macierzL[a][k];
        macierzL[a][k]=tmp[k];
    }
        y=wektor[i];
        wektor[i]=wektor[a];
        wektor[a]=y;

}
void showMacierz(double macierz[][4])
{
    cout<<"_______________________________________________"<<endl;
    for(int i=0 ; i<4 ; i++)
    {
        cout<<"| ";
        for(int j=0 ; j<4 ; j++)
        {
            cout.width(10);
            cout<<macierz[i][j];
        }
        cout<<endl;
    }
    cout<<"_______________________________________________"<<endl;
}



void showVec(double *Vec)
{
    cout<<"\t_________"<<endl;
    for(int i=0 ; i<4 ; i++)
    {
        cout<<"\t| "<<Vec[i]<<endl;
    }
    cout<<"\t_________"<<endl;
}

void dekompozycjaLU(double macierz[][4],double macierzL[][4],double *vec)
{
    cout<<"\tDEKOMPOZYCJA LU"<<endl;
    int i,j,k;
    double a;

    showMacierz(macierz);
    cout<<endl;
    for (i=0 ; i<=3 ; i++)
    {
        if(macierz[i][i]==0)
        {
        cout<<"******** WYBOR CZESCIOWY ********"<<endl;
        wybor_czesciowy(macierz,macierzL,vec,i);
        }
        macierzL[i][i]=1;
        for (j=i+1 ; j<4 ; j++)
        {
            a=macierz[j][i];
            macierzL[j][i]=a/macierz[i][i];
            for (k=i ; k<4 ; k++)
            {
                //cout<<"Wiersz: "<<j+1<<endl;
                //cout<<"Kolumna: "<<k+1<<endl;
                macierz[j][k]-=macierz[i][k]*a/macierz[i][i];
                cout<<"\tMacierzU:"<<endl;
                showMacierz(macierz);
                cout<<"\tMacierz L:"<<endl;
                showMacierz(macierzL);
                cout<<endl<<endl;

            }
        }
    }
}

void wektorY(double *vec, double macierz[][4])
{
    double y[4]={0},suma;
    y[0]=vec[0];
    for(int i=1; i<4 ; i++)
    {
        suma=0;
        for(int j=0 ; j<i ; j++)
        {
            suma-=macierz[i][j]*y[j];
        }
        y[i]=suma+vec[i];
        vec[i]=y[i];
    }
}

void solve(double *vec, double macierz[][4])
{

    double x[4]={0},suma=0;
    x[3]=vec[3]/macierz[3][3];
    for(int i=2; i>=0 ; i--)
    {
        suma=0;
        for(int j=3 ; j>=i ; j--)
        {
            suma=suma-macierz[i][j]*x[j];
        }
        x[i]=(suma+vec[i])/macierz[i][i];
    }
    for (int j=0 ; j<4 ; j++)
        vec[j]=x[j];
}

int main()
{
    cout.precision(5);
    const int N = 4;
	double macierz1 [N][N] =
	{{1.0, -20.0, 30.0, -4.0},
	 {2.0, -40.0, -6.0, 50.0},
	  {9.0, -180.0, 11.0, -12.0},
	   {-16.0, 15.0, 140.0, 13.0}};
	double macierzU[N][N];
	double macierzL[N][N];
	double wektor_w_w[N] =
	{35.0,
	 104.0,
	  -366.0,
	   -354.0};
      for(int i=0 ; i<4 ; i++)
        {
            for(int j=0 ; j<4 ; j++)
            {
                macierzL[i][j]=0.0;
                macierzU[i][j]=macierz1[i][j];

            }
        }

    dekompozycjaLU(macierzU,macierzL,wektor_w_w);


    cout<<"\tmacierz U:"<<endl;

    showMacierz(macierzU);

    cout<<"\tmacierz L:"<<endl;

    showMacierz(macierzL);

    cout<<"\n\twektor B:"<<endl;

    showVec(wektor_w_w);

    cout<<"\n\twektor Y:"<<endl;

    wektorY(wektor_w_w,macierzL);

    showVec(wektor_w_w);

    cout<<"\n\twektor X"<<endl;

    solve(wektor_w_w,macierzU);
    showVec(wektor_w_w);





    return 0;
}

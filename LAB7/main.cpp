#include <iostream>
#include <iomanip>
#include <math.h>

#define N 4
#define MAX_ITERATION 50
#define EPSILON 1e-12
#define szerokosc 15
using namespace std;
void wypelnijM(double m[][N]);
void wypelnijV(double v[]);
void wypelnijx(double x[]);
void showMacierz(double m[][N]);
void showV(double v[]);
void print_step(int i, double x[],double est, double res);
double norma(double vec[]);
void metodaJacobiego(double A[][N],double b[], double x[]);
void metodaGaussaSeidela(double A[][N],double b[], double x[]);
void metodaSOR(double A[][N], double b[], double x[], double omega);

int main()
{
    cout.precision(5);
    double omega=0.5;
    double A[N][N],b[N],x[N];


    //cout<<"\n\tmacierz A:\n"<<endl;
    wypelnijM(A);
    //showMacierz(A);

    //cout<<"\n\tvector b:\n"<<endl;
    wypelnijV(b);
    //showV(b);
    //cout<<endl;

    //cout<<"\n\tvector x:\n"<<endl;
    wypelnijx(x);
    //showV(x);
    //cout<<endl;
    //cout<<endl;

    cout<<"\t Jacobiego:\n"<<endl;
    cout<<"step \tx0\t\tx1\t\tx2\t\tx3\t\t\t\tEst\t\tRes"<<endl;
    print_step(0,x,0,0);
    metodaJacobiego(A,b,x);
    cout<<endl;

    cout<<"\t Gaussa-Seidela:\n"<<endl;
    wypelnijM(A);
    wypelnijV(b);
    wypelnijx(x);
    cout<<"step \tx0\t\tx1\t\tx2\t\tx3\t\t\t\tEst\t\tRes"<<endl;
    cout.width(5);
    cout<<left<<0;
    showV(x);
    cout<<endl;
    metodaGaussaSeidela(A,b,x);

    cout<<endl;
    cout<<"\t SOR:\n"<<endl;
    wypelnijM(A);
    wypelnijV(b);
    wypelnijx(x);
    cout<<"step \tx0\t\tx1\t\tx2\t\tx3\t\t\t\tEst\t\tRes"<<endl;
    cout.width(5);
    cout<<left<<0;
    showV(x);
    cout<<endl;
    metodaSOR(A,b,x,omega);
    return 0;
}

void wypelnijM(double a[][N])
{
    double tab[N][N]={ {100.0,-1.0,2.0,-3.0},
                       {1.0,200.0,-4.0,5.0},
                       {-2.0,4.0,300.0,6.0},
                       {3.0, -5.0, 6.0, 400.0} };
    for (int i=0 ; i<N ; i++)
        for (int j=0 ; j<N ; j++)
            a[i][j]=tab[i][j];
}
void showMacierz(double m[][N])
{
    for (int i=0 ; i<N ; i++)
    {
        for (int j=0 ; j<N ; j++)
        {
            cout.width(szerokosc);
            cout<<m[i][j];
        }
        cout<<endl;
    }
}

void wypelnijV(double v[])
{
    double tab[N]={116.0,-226.0,912.0,-1174.0};
    for (int i=0 ; i<N ; i++)
        v[i]=tab[i];
}


void wypelnijx(double x[])
{
    double tab[N]={2.0,2.0,2.0,2.0};
    for (int i=0 ; i<N ; i++)
        x[i]=tab[i];
}
void showV(double v[])
{
    for (int i=0 ; i<N ; i++)
    {
        cout.width(szerokosc);
        cout<<v[i];
    }

}

void print_step(int i, double tab[],double est, double res)
{
    cout.width(5);
    cout<<left<<i;
    showV(tab);
    cout.width(szerokosc+szerokosc);
    cout<<right<<est;
    cout.width(szerokosc);
    cout<<res<<endl;
}

double norma(double vec[])
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

void metodaJacobiego(double A[][N],double b[], double x[])
{
    int i=0,j=0,k=0;
    double res[N],x_new[N];
    double temp=0,estymator,residuum;
    double roznica[N];
    for(int it=0 ; it<MAX_ITERATION ; it++)
    {
        for(i=0 ; i<N ; i++)
        {
            for(j=0 ; j<N ; j++)
            {
                if(i != j)
                    temp += A[i][j]*x[j];
            }
            x_new[i] = (b[i] - temp) * (1/ A[i][i]);
            temp=0;
        }


        for (i=0 ; i<N ; i++)
            roznica[i]=x_new[i]-x[i];

        estymator=norma(roznica);

        for(i=0 ; i<N ; i++)
        {
            for(j=0 ; j<N ; j++)
            {
                    temp += A[i][j]*x_new[j];
            }
            res[i] = temp-b[i];
            temp=0;
        }
        residuum=norma(res);

        if(estymator<EPSILON || residuum<EPSILON)
            break;
        for(k=0 ; k<N ; k++)
            x[k]=x_new[k];

        print_step(it+1,x_new,estymator,residuum);
    }
}

void metodaGaussaSeidela(double A[][N],double b[], double x[])
{
    double x_poprzedni[N];
    int i=0,j=0;
    double res[N];
    double tmp=0;
    double estymator=0,residuum=0;
    double roznica[N];
    for(int it=0 ; it<MAX_ITERATION ; it++)
    {
        for(i=0 ; i<N ; i++)
        {
            x_poprzedni[i]=x[i];
            for(j=0 ; j<N ; j++)
            {
                if(i != j)
                    tmp += A[i][j]*x[j];
            }
            x[i] = (b[i] - tmp) * (1/ A[i][i]);
            tmp=0;
        }
        for (i=0 ; i<N ; i++)
            roznica[i]=x_poprzedni[i]-x[i];
        estymator=norma(roznica);
        for(i=0 ; i<N ; i++)
        {
            for(j=0 ; j<N ; j++)
            {
                    tmp += A[i][j]*x[j];
            }
            res[i] = tmp-b[i];
            tmp=0;
        }
        residuum=norma(res);

        if(estymator<EPSILON || residuum<EPSILON)
            break;

        print_step(it+1,x,estymator,residuum);
    }
}

void metodaSOR(double A[][N], double b[], double x[], double omega)
{
    double x_poprzedni[N];
    int i=0,j=0;
    double res[N];
    double tmp=0;
    double estymator=0,residuum=0;
    double roznica[N];

    for(int iter=1 ; iter<=MAX_ITERATION ; iter++)
    {
        for (i=0 ; i<N ; i++)
        {
            x_poprzedni[i]=x[i];
            for (j=0 ; j<N ; j++)
                if(i==j)
                    tmp += (1 - (1.0 / omega))*A[i][i] * x[j];
                else
                    tmp += A[i][j] * x[j];
            x[i] = (b[i] - tmp)*(1.0 / A[i][i] * omega);
            tmp=0;
        }

        //Estymator
        for (i=0 ; i<N ; i++)
            roznica[i]=x_poprzedni[i]-x[i];
        estymator=norma(roznica);

        //Residuum
        for(i=0 ; i<N ; i++)
        {
            for(j=0 ; j<N ; j++)
            {
                    tmp += A[i][j]*x[j];
            }
            res[i] = tmp-b[i];
            tmp=0;
        }
        residuum=norma(res);

        if(estymator<EPSILON || residuum<EPSILON)
            break;

        print_step(iter,x,estymator,residuum);
    }
}



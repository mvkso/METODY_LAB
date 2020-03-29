#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <limits>

using namespace std;

double taylorSeries(double x)
{
    double value = 1;
    double result = 1;

    for (int i = 1; i < 30; i++)
    {
        value *= -(x / ((double)i + 1.0));
        result += value;
    }
    return result;
}

double function(double x)
{
    return (1.0 - exp(-x)) / x;
}

int main()
{
    double x, result, absoluteError, log, exactValue, relativeError;

    ifstream f_in;
    f_in.open("dane.txt");
    if (!f_in.good())
        return 1;

    ofstream f_out;
    f_out.open("wynik.txt", ios::out);
    if (!f_out.good())
        return 1;

    fstream wykres;
    wykres.open("rysuj.dat", ios::out);
    if (!wykres.good())
        return 1;

    cout   << scientific << setprecision(20);
    f_in   >> scientific >> setprecision(20);
    f_out  << scientific << setprecision(20);

    cout << setw(27) << "x ";
    cout << setw(27) << "Wartosc dokladna ";
    cout << setw(27) << "Wynik ";
    cout << setw(27) << "Roznica";
    cout << setw(27) << "Blad wzgledny ";
    cout << endl;

    while (!f_in.eof())
    {
        f_in >> log >> x >> exactValue;

        result        = function(x);
        absoluteError = fabs(exactValue - result);


        if (absoluteError > numeric_limits<double>::epsilon())
        {
            result        = taylorSeries(x);
            absoluteError = fabs(exactValue - result);
        }

        relativeError = fabs(absoluteError / exactValue);

        f_out << log << " " << log10(relativeError) << endl;

        cout << setw(27) << x<< "  ";
        cout << setw(27) << exactValue<< "  ";
        cout << setw(27) << result<< "  ";
        cout << setw(27) << absoluteError<< "  ";
        cout << setw(27) << relativeError<< "  ";
        cout << endl;
    }

    f_in.close();
    f_out.close();

    return 0;
}

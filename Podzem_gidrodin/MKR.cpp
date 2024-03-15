#include "MKR.h"
#include <iostream>
#include<fstream>

using namespace std;

void print1(vector<double>& r, vector<double>& p) {
    for (int i = 0; i < p.size(); i++)
    {
        cout << "r ==" << r[i] << "    p== " << p[i] << endl;
    }
}
void print1(vector<double>& tetta) {
    for (int i = 0; i < tetta.size(); i++)
    {
        cout << "tetta== " << tetta[i] << endl;
    }
}

void out_file1(vector<double>& r, vector<double>& p) {
    ofstream fout("P.dat");
    fout << "r p" << endl;
    for (int i = 0; i < r.size(); i++)
    {
        fout << r[i] << " " << p[i] << endl;
    }
    fout.close();
}
static std::vector<double> solveTridiagonalMatrix(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& d) {
    int n = d.size();
    std::vector<double> alpha(n), beta(n);
    std::vector<double> x(n);

    alpha[0] = b[0];
    beta[0] = d[0] / alpha[0];

    for (int i = 1; i < n; ++i) {
        alpha[i] = b[i] - a[i] * c[i - 1] / alpha[i - 1];
        beta[i] = (d[i] - a[i] * beta[i - 1]) / alpha[i];
    }

    x[n - 1] = beta[n - 1];

    for (int i = n - 2; i >= 0; --i) {
        x[i] = beta[i] - c[i] * x[i + 1] / alpha[i];
    }

    return x;
}
void MKR(const double R, const double r, const double p_0, const double p_1, const int N) {
    vector<double> ri;
    const double rc = r / R;
    ri.resize(N + 1);
    for (double i = 0; i < N + 1; i++)
    {
        //r_i.push_back(rc + ((1 - rc) / N) * i);
        ri[i] = rc + ((1 - rc) / N) * i;
        //r_i[i] = rc + (1 - rc)* pow(1.1, i- (double)N);
    }
   
    vector<double> k(N, 1); // к-констант
    vector<double> a{ 0 };
    vector<double> b{ (ri[1] * k[0]) / (ri[1] - ri[0]) + (ri[2] * k[1]) / (ri[2] - ri[1]) };
    vector<double> c{ -(ri[2] * k[1]) / (ri[2] - ri[1]) };
    vector<double> d(N);


    for (int i = 1; i < N ; i++)
    {
        a.push_back(-(ri[i] * k[i-1]) / (ri[i] - ri[i - 1]));
        b.push_back((ri[i] * k[i-1]) / (ri[i] - ri[i - 1]) + (ri[i + 1] * k[i]) / (ri[i+1] - ri[i ]));
        c.push_back(-(ri[i + 1] * k[i]) / (ri[i+1] - ri[i]));
    }
    a.push_back(-(ri[N - 1] * k[N - 2]) / (ri[N - 1] - ri[N - 2]));
    b.push_back((ri[N - 1] * k[N - 2]) / (ri[N - 1] - ri[N - 2]) + (k[N-1]) / (1 - ri[N - 1]));
    c.push_back(0);
    d.push_back(k[N-1] / (1 - ri[N - 1]));
    vector<double> p = solveTridiagonalMatrix(a, b, c, d);

    //print(ri, p);
    out_file1(ri, p);
}


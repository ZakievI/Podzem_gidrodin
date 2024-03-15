#include "Zadanie#1.h"
#include "iostream"
#include "fstream"
#include "vector"
#include "string"
#include <iomanip>
#include <limits>
#include <numbers>

struct point {
    point(double x, double y): x{x},y{y}{}
    point(double x):x{x}{}
    point(){}
    point operator*(double a) {
        point b;
        b.x = this->x * a;
        b.y = this->y * a;
        return b;
    }
    point operator/(double a) {
        point b;
        b.x = this->x / a;
        b.y = this->y / a;
        return b;
    }
    
    double x;
    double y;
};
template <typename T, typename N>
std::vector<T> operator *(const std::vector<T>& lhs, const N& rhs)
{
    std::vector<T> result;
    result.resize(lhs.size());
    std::copy(lhs.cbegin(), lhs.cend(), result.begin());
    for (auto& t : result)
        t *= rhs;
    return result;
}
template <typename T, typename N>
std::vector<T> operator /(const std::vector<T>& lhs, const N& rhs)
{
    std::vector<T> result;
    result.resize(lhs.size());
    std::copy(lhs.cbegin(), lhs.cend(), result.begin());
    for (auto& t : result)
        t = t/rhs;
    return result;
}
template <typename T, typename N>
std::vector<T> operator +(const std::vector<T>& lhs, const N& rhs)
{
    std::vector<T> result;
    result.resize(lhs.size());
    std::copy(lhs.cbegin(), lhs.cend(), result.begin());
    for (auto& t : result)
        t = t + rhs;
    return result;
}
template <typename T, typename N>
std::vector<T> operator -(const std::vector<T>& lhs, const N& rhs)
{
    std::vector<T> result;
    result.resize(lhs.size());
    std::copy(lhs.cbegin(), lhs.cend(), result.begin());
    for (auto& t : result)
        t = t - rhs;
    return result;
}
// точое значение для давления для задачи 1.2 (2-ой подпункт)
double exact_P_1(double r){
    return 1-r;
}
// точое значение для давления для задачи 1.2 (2-ой подпункт)
double exact_P_12_2(double r, double r_0) {
    if (r <= 0.5)
    {
        return (1 - 0.1 * (r)-(1-0.1)* (0.5)) / (1 - 0.1 * r_0 - (1 - 0.1) * (0.5));
    }
    else {
        return (r - 1) / (0.1 * (r_0-0.5) - 0.5); 
    }
}
// точое значение для давления для задачи 1.2 (3-ой подпункт)
double exact_P_12_3(double r, double r_0) {
    if (r <= 0.75)
    {
        return (-1 + 0.0526316 * r + 1 * 0.75 - 0.0526316 * 0.75) / (-1 + 0.0526316 * r_0 + 0.75 - 0.0526316 * 0.75);
    }
    else {
        return (-1 + r) / (-1 + 0.0526316 * r_0 + 0.75 - 0.0526316 * 0.75);
    }
}
double compute_norm(std::vector<long double> p, std::vector<point> Mh,int z) {
    double Norm=0.0;
    
    switch (z)
    {
    case(1):
    {
        for (int i = 0; i < p.size(); i++)
        {
            Norm += (p[i] -exact_P_1(Mh[i].x)) * (p[i] - exact_P_1(Mh[i].x));
        }
        break;
    }
    case(2):
    {
        for (int i = 0; i < p.size(); i++)
        {
            Norm += (p[i] - exact_P_12_2(Mh[i].x, Mh[0].x)) * (p[i] - exact_P_12_2(Mh[i].x, Mh[0].x));
        }
        break;
    }
    case(3):
    {
        for (int i = 0; i < p.size(); i++)
        {
            Norm += (p[i] - exact_P_12_3(Mh[i].x, Mh[0].x))* (p[i] - exact_P_12_3(Mh[i].x, Mh[0].x));
        }
        break;
    }
    default:
        break;
    }
  
    return sqrt(Norm / p.size());
}
double compute_taim(std::vector<long double> v, std::vector<point> Mh) {
    double taim = 0.0;
    for (int i = 0; i < v.size(); i++)
    {
        taim += (Mh[i + 1].x - Mh[i].x) / abs(v[i]);
    }
    return taim;
}
std::vector<double> compute_k(std::vector<point> Mh, int i)
{
    std::vector<double> k1;
    switch (i)
    {
    case(1): {
        std::vector<double> k(Mh.size() - 1, 1.0);
        k1 = k;
        break;
    }
    case(2): {
        for (int i = 1; i < Mh.size(); i++)
        {
            if (Mh[i].x <= 0.5)
            {
                k1.push_back(1.0);
            }
            else {
                k1.push_back(0.1);
            }
        }
        break;
    }
    case(3): {
        for (int i = 1; i < Mh.size(); i++)
        {
            if (Mh[i].x <= 0.75)
            {
                k1.push_back(1.0);
            }
            else {
                k1.push_back(0.0526316);
            }
        }
        break;
    }
    default:
        break;
    }
    return k1;
}
template <typename T>
void out_file(std::vector<point>& r, std::vector<T>& p, std::vector<T>& v,int z) {
    std::ofstream fout("P.dat");
    fout << "x p_n v p_e" << std::endl;
    switch (z)
    {
    case(1):
    {
        for (int i = 0; i < r.size() - 1; i++)
        {
            fout << std::setprecision(10) << r[i].x << " " << p[i] << " " << v[i]<< " " <<exact_P_1(r[i].x) << std::endl;
        }
        fout << std::setprecision(10) << r[r.size() - 1].x << " " << p[r.size() - 1] << " " << exact_P_1(r[r.size() - 1].x) << std::endl;
        break;
    }
    case(2):
    {
        for (int i = 0; i < r.size() - 1; i++)
        {
            fout << std::setprecision(10) << r[i].x << " " << p[i] << " " << v[i] << " " << exact_P_12_2(r[i].x, r[0].x) << std::endl;
        }
        fout << std::setprecision(10) << r[r.size() - 1].x << " " << p[r.size() - 1] << " "<<exact_P_12_2(r[r.size() - 1].x, r[0].x) << std::endl;
        break;
    }
    case(3):
    {
        for (int i = 0; i < r.size() - 1; i++)
        {
            fout << std::setprecision(10) << r[i].x << " " << p[i] << " " << v[i] << " " << exact_P_12_3(r[i].x, r[0].x) << std::endl;
        }
        fout << std::setprecision(10) << r[r.size() - 1].x << " " << p[r.size() - 1] << " " << exact_P_12_3(r[r.size() - 1].x, r[0].x) << std::endl;
        break;
    }
    default:
        break;
    }

    fout.close();
}
void create_mesh(int N) {
    std::ofstream fout("mesh.dat");
    double h_0 = 0.0;
    double h_1 = 100.0;
    for (double i = 0; i < N + 1; i++)
    {
        fout << h_0+(h_1-h_0)*i/N<< std::endl;
        
    };
    fout.close();
}
std::vector<long double> velocyte(std::vector<long double> p,std::vector<point> Mh, std::vector<double> k) {
    std::vector<long double> v;
    v.resize(p.size()-1);
    for (int i = 0; i < v.size(); i++)
    {
        v[i] = -(p[i + 1] - p[i]) / (Mh[i + 1].x - Mh[i].x);
    }
    return v;
}


static std::vector<long double> solveTridiagonalMatrix(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& d) {
    int n = d.size();
    std::vector<double> alpha(n), beta(n);
    std::vector<long double> x(n);

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

void Zadanie_1(const double p_0, const double p_1, std::string mesh) {
    create_mesh(100);
    int zd = 1;
    static std::vector<point> Mh;
    std::ifstream mesh_(mesh);
    if (mesh_.is_open())
    {
        double x;
        while (mesh_>>x)
        {
            Mh.push_back(point(x));
        }
    }
    int N = Mh.size();
    double L = (Mh[N - 1].x - Mh[0].x);
    Mh = Mh / (Mh[N - 1].x - Mh[0].x);
    std::vector<double> k=compute_k(Mh,zd); // к-констант
    std::vector<double> a{ 0 };
    std::vector<double> b{ -(k[0]) / (Mh[1].x - Mh[0].x) - (k[1]) / (Mh[2].x - Mh[1].x) };
    std::vector<double> c{ (k[1]) / (Mh[2].x - Mh[1].x) };
    std::vector<double> d(N-3);
    d[0] = -k[0] / ((Mh[1].x - Mh[0].x));

    for (int i = 1; i < N-3; i++)
    {
        a.push_back((k[i]) / (Mh[i+1].x - Mh[i].x));
        b.push_back(-(k[i + 1]) / (Mh[i+2].x - Mh[i + 1].x) - (k[i]) / (Mh[i + 1].x - Mh[i].x));
        c.push_back((k[i+1]) / (Mh[i + 2].x - Mh[i+1].x));
    }

    a.push_back((k[N - 3]) / (Mh[N - 2].x - Mh[N -3].x));
    b.push_back(-(k[N - 3]) / (Mh[N - 2].x - Mh[N - 3].x) - (k[N - 2]) / (Mh[N-1].x - Mh[N - 2].x));
    c.push_back(0);
    d.push_back(0);
    std::vector<long double> p = solveTridiagonalMatrix(a, b, c, d);
    auto iter = p.cbegin();
    p.emplace(iter, 1);
    p.push_back(0);
    //p = p * p_1;
    //print1(Mh, p);
    mesh_.close();
    double norma = compute_norm(p, Mh,zd);
    bool pr=0;
    if (pr)
    {
        for (int i = 0; i < N; i++)
        {
            p[i] = (p_0 - p_1) * p[i];
            Mh[i].x = Mh[i].x * L;
        }
    }
    std::vector<long double> v = velocyte(p, Mh,k);
    double t_n=0.0;
    double t_a=0.0;
    double mu = 1e-3;
    double m = 0.2;
    t_a = (0.2 * pow(10,-3) * L * L) / (pow(10,-12)* (p_0 - p_1));
    t_n = compute_taim(v, Mh)*k[0]*m*abs(p_0-p_1)/mu/L;
    out_file(Mh, p,v,2);

    
    double tau = 0.001;

    std::ofstream fout("C.dat");
    for (int i = 0; i < Mh.size(); i++)
    {
        fout << Mh[i].x << "  ";

    }
    int l=0;
    setlocale(LC_ALL, "Russian");
    std::cout << "---------------------------------------Welcome-----------------------------------------" << std::endl;
    std::cout << "---------------------------------------------------------------------------------------" << std::endl;
    std::cout << "Norma:=   " << norma << std::endl;
    std::cout << "---------------------------------------------------------------------------------------" << std::endl;
    std::cout << "Taim_a:=   " << t_a << " " << "Taim_n:=   " << t_n << std::endl;
    bool stop=0;
    while (!stop)
    {
        std::cout << "---------------------------------------------------------------------------------------" << std::endl;
        std::cout << "#0 (выйти) #1 (график для разных итераций) #2 (график для различных коэффициентов дифузии)" << std::endl;
        std::cout << "---------------------------------------------------------------------------------------" << std::endl;
        std::cin >> l;
        switch (l)
        {
        case 0:
        {
            stop = 1;
            break; 
        }
        case 1: {
            fout << std::endl;
            std::vector<double> c1(N - 1, 0);
            auto iter_c = c1.cbegin();
            c1.emplace(iter_c, 1);
            double D = 1e-5;
            for (int k = 1; k < 80; k++)
            {
                for (int i = 1; i < c1.size() - 1; i++)
                {
                    c1[i] = c1[i] + (tau)*D * ((c1[i + 1] - c1[i]) / (Mh[i + 1].x - Mh[i].x) - (c1[i] - c1[i - 1]) / (Mh[i].x - Mh[i - 1].x)) / (m * (Mh[i].x - Mh[i - 1].x)) - tau * (c1[i] * v[i] - c1[i - 1] * v[i - 1]) / (m * (Mh[i].x - Mh[i - 1].x));
                }
                c1[c1.size() - 1] = c1[c1.size() - 2];
                if (k % 5 == 0)
                {
                    double err = 0.0;
                    std::cout << "D__" << D << "__iter__" << k << "__taim__" << tau * k  << std::endl;
                    std::cout << "---------------------------------------------------------------------------------------" << std::endl;

                    for (int i = 0; i < c1.size(); i++)
                    {
                        std::cout << c1[i] << "  ";
                        fout << c1[i] << "  ";
                        //err += c1[i] - m * exp(-pow(Mh[i].x - v[i] * tau * k, 2) / 4 / D / tau / k)/2/sqrt(3.1415*D*tau*k);
                    }
                    std::cout << std::endl;
                  

                    //std::cout << c1[c1.size()-1] << "  ";
                    //std::cout << std::endl;
                    //std::cout << "err: "<<err;
                    std::cout << std::endl;
                    fout << std::endl;
                }

            }
            break;
        }
        case 2: {
            int N_;
            std::cout << "введи число итераций!!!" << std::endl;
            std::cin >> N_;
            for (double D : { 1e-1, 1e-2, 1e-5, 0.0, })
            {
                std::vector<double> c1(N - 1, 0);
                auto iter_c = c1.cbegin();
                c1.emplace(iter_c, 1);
                for (int k = 1; k < N_ + 1; k++)
                {
                    for (int i = 1; i < c1.size() - 1; i++)
                    {
                        c1[i] = c1[i] + (tau)*D * ((c1[i + 1] - c1[i]) / (Mh[i + 1].x - Mh[i].x) - (c1[i] - c1[i - 1]) / (Mh[i].x - Mh[i - 1].x)) / (m * (Mh[i].x - Mh[i - 1].x)) - tau * (c1[i] * v[i] - c1[i - 1] * v[i - 1]) / (m * (Mh[i].x - Mh[i - 1].x));
                    }
                    c1[c1.size() - 1] = c1[c1.size() - 2];

                }
                double err = 0.0;
                std::cout << "---------------------------------------------------------------------------------------" << std::endl;
                std::cout << "D__" << D << "__iter__" << N_ << "__taim__" << tau * N_ << std::endl;
                std::cout << "---------------------------------------------------------------------------------------" << std::endl;
                for (int i = 0; i < c1.size(); i++)
                {
                    std::cout << c1[i] << "  ";
                    fout << c1[i] << "  ";
                    //err += c1[i] - m * exp(-pow(Mh[i].x - v[i] * tau * k, 2) / 4 / D / tau / k)/2/sqrt(3.1415*D*tau*k);
                }
                std::cout << std::endl;
                

                //std::cout << c1[c1.size()-1] << "  ";
                //std::cout << std::endl;
                //std::cout << "err: "<<err;
                std::cout << std::endl;
                fout << std::endl;
            }
            break;
        }
        default:
        {std::cout << "error" << std::endl;
        break;
        }
        }
    }
    
    
    fout.close();


}


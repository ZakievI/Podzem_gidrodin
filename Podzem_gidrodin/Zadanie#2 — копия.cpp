#include "Zadanie#1.h"
#include "iostream"
#include "fstream"
#include "vector"
#include "string"
#include <iomanip>
#include <limits>
#include <numbers>

constexpr auto PI =3.14159265359;
constexpr auto mu=1e-3;
constexpr auto m = 0.2;
constexpr auto k_0 = 1e-12;
bool poprav_kof=1; //включать или отключать поправочный коэфициент 

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
namespace zd2 {
    double compute_norm(std::vector<long double> p, std::vector<point> Mh) {
        double Norm = 0;
        double exi = 0.0;
        for (int i = 0; i < p.size(); i++)
        {
            exi = -(log(Mh[i].x) - log(Mh[0].x)) / (log(Mh[0].x));
            Norm += (p[i] - exi) * (p[i] - exi);
        }
        return sqrt(Norm / p.size());
    }
    double compute_taim(std::vector<long double> v, std::vector<point> Mh, std::vector<double> k) {
        double taim = 0.0;
        for (int i = 0; i < v.size(); i++)
        {
            taim += k[i] * (Mh[i + 1].x - Mh[i].x) / abs(v[i]);
        }
        return taim;
    }
    double exact_P_1(double r, double r_0) {
        return -(log(r) - log(r_0)) / (log(r_0));
    }
    double exact_Q_1(double r, double r_0) {
        return  -2 * PI / log(r_0);
    }
    double velocity_u_0(double r_0) {
        return -1 / (r_0*log(r_0));
    }
    double exact_P_2(double r, double r_0) {
        return 1+exact_Q_1(r,r_0)*r_0*log(r);
    }
    double exact_Q_2(double r, double r_0) {
        return  exact_Q_1(r,r_0)*2 * PI *(r_0);
    }
    std::vector<long double> compout_debit(std::vector<long double> p, std::vector<point> Mh, std::vector<double> k) {
        std::vector<long double> Q;
        Q.resize(p.size() - 1);
        for (int i = 0; i < p.size() - 1; i++)
        {
            Q[i] = k[i] * (p[i + 1] - p[i]) * 2 * PI * Mh[i].x / ((Mh[i + 1].x - Mh[i].x));
        }
        return Q;
    }
    template <typename T>
    void out_file(std::vector<point>& r, std::vector<T>& p, std::vector<T>& v, std::vector<double> k) {
        std::ofstream fout("P.dat");
        fout << "x p_n v p_e Q Q_A" << std::endl;
        for (int i = 0; i < r.size() - 1; i++)
        {
            fout << std::uppercase << std::scientific << std::setprecision(10) << r[i].x << " " << p[i] << " " << v[i] << " " << exact_P_1(r[i].x,r[0].x) << " " << compout_debit(p, r, k)[i]<< " " << exact_Q_1(r[i].x, r[0].x) << std::endl;
        }
        fout << std::uppercase << std::scientific << std::setprecision(10) << r[r.size() - 1].x << " " << p[r.size() - 1] << " . " << -(log(r[r.size() - 1].x) - log(r[0].x)) / (log(r[0].x)) << std::endl;
        fout.close();
    }
    void create_mesh(int N) {
        std::ofstream fout("mesh.dat");
        double h_0 = 0.1;
        double h_1 = 100.0;
        double m;
        for (double i = 0; i < N + 1; i++)
        {
           // fout << h_0 + (h_1 - h_0) * i / N << std::endl;
            m = h_0 * pow(h_1 / h_0, i / N);
            fout << m << " ";

        };
        fout.close();
    }
    std::vector<long double> velocyte(std::vector<long double> p, std::vector<point> Mh, std::vector<double> k) {
        std::vector<long double> v;
        v.resize(p.size() - 1);
        for (int i = 0; i < v.size(); i++)
        {
            v[i] = -k[i] * (p[i + 1] - p[i]) / (Mh[i + 1].x - Mh[i].x);
        }
        return v;
    }
    static std::vector<long double> solveTridiagonalMatrix(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& d) {
        size_t n = d.size();
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
    static std::vector<long double> compute_usl_1_roda(std::vector<point> Mh, std::vector<double> k) {
        int N = Mh.size();
        std::vector<double> Q;
        if (poprav_kof) {
            for (int i = 0; i < k.size()-1; i++)
            {
                Q.push_back((Mh[i + 1].x - Mh[i].x)/ (Mh[i].x * log(Mh[i + 1].x/ Mh[i].x)));
            }
        }
        else {
            std::vector<double> Q1(N, 1.0);
            Q = Q1;
        }
        std::vector<double> a{ 0 };
        std::vector<double> b{ -(Mh[0].x * k[0]* Q[0]) / (Mh[1].x - Mh[0].x) - (Mh[1].x * k[1] * Q[1]) / (Mh[2].x - Mh[1].x) };
        std::vector<double> c{ (Mh[1].x * k[1] * Q[1]) / (Mh[2].x - Mh[1].x) };
        std::vector<double> d(N - 3);

        for (int i = 2; i < N - 2; i++)
        {
            a.push_back((Mh[i-1].x * k[i - 1] * Q[i-1]) / (Mh[i].x - Mh[i - 1].x));
            b.push_back(-(Mh[i-1].x * k[i - 1] * Q[i-1]) / (Mh[i].x - Mh[i - 1].x) - (Mh[i].x * k[i] * Q[i]) / (Mh[i + 1].x - Mh[i].x));
            c.push_back((Mh[i].x * k[i] * Q[i]) / (Mh[i + 1].x - Mh[i].x));
        }

        a.push_back((Mh[N - 3].x * k[N - 3] * Q[N - 3]) / (Mh[N - 2].x - Mh[N - 3].x));
        b.push_back(-(Mh[N - 3].x * k[N - 3] * Q[N - 3]) / (Mh[N - 2].x - Mh[N - 3].x) - (Mh[N - 2].x * k[N - 2] * Q[N - 2]) / (Mh[N - 1].x - Mh[N - 2].x));
        c.push_back(0);
        d.push_back(-Mh[N - 2].x * k[N - 2] * Q[N - 2] / (Mh[N - 1].x - Mh[N - 2].x));

        std::vector<long double> p = solveTridiagonalMatrix(a, b, c, d);
        auto iter = p.cbegin();
        p.emplace(iter, 0);
        p.push_back(1);
        return p;
    }
    static std::vector<long double> compute_usl_1_2_roda(std::vector<point> Mh, std::vector<double> k) {
        int N = Mh.size();
        std::vector<double> Q;
        if (poprav_kof) {
            for (int i = 0; i < k.size() - 1; i++)
            {
                Q.push_back((Mh[i + 1].x - Mh[i].x) / (Mh[i].x * log(Mh[i + 1].x / Mh[i].x)));
            }
        }
        else {
            std::vector<double> Q1(N, 1.0);
            Q = Q1;
        }
        std::vector<double> a{ 0 };
        std::vector<double> b{- (Mh[1].x * k[1] * Q[1]) / (Mh[2].x - Mh[1].x)};
        std::vector<double> c{ (Mh[1].x * k[1] * Q[1]) / (Mh[2].x - Mh[1].x) };
        std::vector<double> d(N - 3);
        d[0] = Mh[0].x * k[0] * Q[0] * velocity_u_0(Mh[0].x);

        for (int i = 2; i < N - 2; i++)
        {
            a.push_back((Mh[i - 1].x * k[i - 1] * Q[i - 1]) / (Mh[i].x - Mh[i - 1].x));
            b.push_back(-(Mh[i - 1].x * k[i - 1] * Q[i - 1]) / (Mh[i].x - Mh[i - 1].x) - (Mh[i].x * k[i] * Q[i]) / (Mh[i + 1].x - Mh[i].x));
            c.push_back((Mh[i].x * k[i] * Q[i]) / (Mh[i + 1].x - Mh[i].x));
        }

        a.push_back((Mh[N - 3].x * k[N - 3] * Q[N - 3]) / (Mh[N - 2].x - Mh[N - 3].x));
        b.push_back(-(Mh[N - 3].x * k[N - 3] * Q[N - 3]) / (Mh[N - 2].x - Mh[N - 3].x) - (Mh[N - 2].x * k[N - 2] * Q[N - 2]) / (Mh[N - 1].x - Mh[N - 2].x));
        c.push_back(0);
        d.push_back(-Mh[N - 2].x * k[N - 2] * Q[N - 2] / (Mh[N - 1].x - Mh[N - 2].x));
        std::vector<long double> p = solveTridiagonalMatrix(a, b, c, d);
        auto iter = p.cbegin();
        p.emplace(iter, p[0]-(Mh[1].x-Mh[0].x)*(velocity_u_0(Mh[0].x)) );
        p.push_back(1);
        return p;
    }
}
using namespace zd2;
void Zadanie_2(const double p_0, const double p_1, std::string mesh) {
    setlocale(LC_ALL, "Russian");
    std::cout << "---------------------------------------Welcome-----------------------------------------" << std::endl;
    std::cout << "---------------------------------------------------------------------------------------" << std::endl;
    static std::vector<point> Mh;
    std::ifstream mesh_(mesh);
    for (int i : {10, 100, 1000, 10000})
    {
        std::ifstream mesh_(mesh);
        zd2::create_mesh(i);
        if (mesh_.is_open())
        {
            double x;
            while (mesh_ >> x)
            {
                Mh.push_back(point(x));
            }
        }
        mesh_.close();
        int N = Mh.size();
        double L = (Mh[N - 1].x);
        Mh = Mh / L;
        std::vector<double> k(N, 1); // к-констант
        std::vector<long double> p = compute_usl_1_2_roda(Mh, k);
        std::cout << "norma " << i << " " << compute_norm(p, Mh) << std::endl;
        Mh.clear();
    }
    zd2::create_mesh(100);
    if (mesh_.is_open())
    {
        double x;
        while (mesh_>>x)
        {
            Mh.push_back(point(x));
        }
    }
    mesh_.close();
    int N = Mh.size();
    double L = (Mh[N - 1].x);
    Mh = Mh / L;
    std::vector<double> k(N, 1); // к-констант
    std::vector<long double> p = compute_usl_1_2_roda(Mh, k);
    double norma = compute_norm(p, Mh);
    std::vector<long double> v = velocyte(p, Mh,k);
    double t_n=0.0;
    double t_a=0.0;
    t_a = (m * mu * log(Mh[N-1].x/Mh[0].x)* L* L*(Mh[N - 1].x - Mh[0].x) * (Mh[N - 1].x + Mh[0].x)) / (2*k_0 * abs(p_0 - p_1));
    t_n = compute_taim(v, Mh,k)*m*abs(p_0-p_1)/mu/L;
    zd2::out_file(Mh, p,v,k);
    double tau = 0.001;
    std::ofstream fout("C.dat");
    for (int i = 0; i < Mh.size(); i++)
    {
        fout << Mh[i].x << "  ";

    }
    int l=0;

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


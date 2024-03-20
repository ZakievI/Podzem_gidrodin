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
    // подсчет времени численно 
    double compute_taim(std::vector<long double> v, std::vector<point> Mh, std::vector<double> k) {
        double taim = 0.0;
        for (int i = 0; i < v.size(); i++)
        {
            taim += (Mh[i + 1].x - Mh[i].x) / abs(v[i]);
        }
        return taim;
    }
    // точоне решение для давления задача 2.1 (когда граничныеусловия 1-ого рода) , 2.2(1-ая сетка)
    double exact_P_1(double r, double r_0) {
        return -(log(r) - log(r_0)) / (log(r_0));
    }
    // дебит скважины (пока хз правильно подсчитан или нет)
    double exact_Q_1(double r, double r_0) {
        return  -2 * PI / log(r_0);
    }
    // точное значение скорости для задачи 2.1, 2.2(1-ая сетка)
    double velocity_u_0(double r_0) {
        return 1 / (r_0*log(r_0));
    }
    // точное значение давления для задачи 2.1 (когда граничные условия 1-ого и 2-ого рода)
    double exact_P_2(double r, double r_0) {
        return 1+velocity_u_0(r_0)* r_0*log(r);
    }
    // дебит для постановки выше 
    double exact_Q_2(double r, double r_0) {
        return  exact_Q_1(r,r_0)*2 * PI *(r_0);
    }
    // точое значение для давления для задачи 2.2 (2-ой подпункт)
    double exact_P_22_2(double r,double r_0){
        if (r<=0.5)
        {
            return 0.1 * (-log(r) + log(r_0)) / (0.1 * log(r_0) + 0.9 * log(0.5));
        }
        else {
            return (-log(r)+0.1*log(r_0)+0.9*log(0.5)) / (0.1 * log(r_0) + 0.9 * log(0.5));
        }
    }
    // точое значение для давления для задачи 2.2 (3-ой подпункт)
    double exact_P_22_3(double r, double r_0) {
        if (r <= 0.75)
        {
            return 0.0440824 * (-log(r) + log(r_0)) / (0.0440824 * log(r_0) + 1.0 * log(0.75)- 0.0440824* log(0.75));
        }
        else {
            return (-log(r) + 0.0440824 * log(r_0) + log(0.75)- 0.0440824*log(0.75)) / (0.0440824 * log(r_0) + (1- 0.0440824)*log(0.75));
        }
    }
    // vr отвечает за вариант поля проницаемости 
    
    template <typename T>
    void out_file(std::vector<point>& r, std::vector<T>& p, std::vector<T>& v, std::vector<double> k, int zd) {
        std::ofstream fout("P.dat");
        fout << "x p_n v p_e Q Q_A" << std::endl;
        switch (zd)
        {
        case(1):
        {
            for (int i = 0; i < r.size() - 1; i++)
            {
                fout << std::setprecision(10) << r[i].x << " " << p[i] << " " << v[i] << " " << exact_P_1(r[i].x, r[0].x) /*<< " " << compout_debit()*/ << " " << exact_Q_1(r[i].x, r[0].x) << std::endl;
            }
            fout << std::setprecision(10) << r[r.size() - 1].x << " " << p[r.size() - 1] << " . " << -(log(r[r.size() - 1].x) - log(r[0].x)) / (log(r[0].x)) << std::endl;
            break;
        }
        case(2):
        {
            for (int i = 0; i < r.size() - 1; i++)
            {
                fout << std::setprecision(10) << r[i].x << " " << p[i] << " " << v[i] << " " << exact_P_22_2(r[i].x, r[0].x)<< std::endl;
            }
            fout << std::setprecision(10) << r[r.size() - 1].x << " " << p[r.size() - 1] << " . " << -(log(r[r.size() - 1].x) - log(r[0].x)) / (log(r[0].x)) << std::endl;
            break;
        }
        case(3):
        {
            for (int i = 0; i < r.size() - 1; i++)
            {
                fout << std::setprecision(10) << r[i].x << " " << p[i] << " " << v[i] << " " << exact_P_22_3(r[i].x, r[0].x)<< std::endl;
            }
            fout << std::setprecision(10) << r[r.size() - 1].x << " " << p[r.size() - 1] << " . " << -(log(r[r.size() - 1].x) - log(r[0].x)) / (log(r[0].x)) << std::endl;
            break;
        }
        case(4):
        {
            for (int i = 0; i < r.size() - 1; i++)
            {
                fout << std::setprecision(10) << r[i].x << " " << p[i] << " " << v[i] << " " << exact_P_2(r[i].x, r[0].x) << std::endl;
            }
            fout << std::setprecision(10) << r[r.size() - 1].x << " " << p[r.size() - 1] << " . " << -(log(r[r.size() - 1].x) - log(r[0].x)) / (log(r[0].x)) << std::endl;
            break;
        }
        default:
            break;
        }
        fout.close();
        
    }
    void create_mesh(int N) {
        std::ofstream fout("mesh.dat");
        double h_0 = 0.1;
        double h_1 = 100.0;
        double m;
        for (double i = 0; i < N; i++)
        {
            fout << h_0 + (h_1 - h_0) * i / N << std::endl;//-регулярная сетка
            //fout << h_0 * pow(h_1 / h_0, i / (N - 1)) << " ";//-логарифмическая сетка

        };
        fout.close();
    }
    std::vector<long double> velocyte(std::vector<long double> p, std::vector<point> Mh, std::vector<double> k, double L) {
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
    static std::vector<long double> compute_usl_1_roda(double p_0, double p_1, std::vector<point> Mh, std::vector<double> k) {
        int N = Mh.size();
        std::vector<double> Q;
        if (poprav_kof) {
            for (int i = 0; i < Mh.size()-1; i++)
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
        d[0] = -p_0*(Mh[0].x * k[0] * Q[0]) / ((Mh[1].x - Mh[0].x));

        for (int i = 2; i < N - 2; i++)
        {
            a.push_back((Mh[i-1].x * k[i - 1] * Q[i-1]) / (Mh[i].x - Mh[i - 1].x));
            b.push_back(-(Mh[i-1].x * k[i - 1] * Q[i-1]) / (Mh[i].x - Mh[i - 1].x) - (Mh[i].x * k[i] * Q[i]) / (Mh[i + 1].x - Mh[i].x));
            c.push_back((Mh[i].x * k[i] * Q[i]) / (Mh[i + 1].x - Mh[i].x));
        }

        a.push_back((Mh[N - 3].x * k[N - 3] * Q[N - 3]) / (Mh[N - 2].x - Mh[N - 3].x));
        b.push_back(-(Mh[N - 3].x * k[N - 3] * Q[N - 3]) / (Mh[N - 2].x - Mh[N - 3].x) - (Mh[N - 2].x * k[N - 2] * Q[N - 2]) / (Mh[N - 1].x - Mh[N - 2].x));
        c.push_back(0);
        d.push_back(-p_1*Mh[N - 2].x * k[N - 2] * Q[N - 2] / (Mh[N - 1].x - Mh[N - 2].x));

        std::vector<long double> p = solveTridiagonalMatrix(a, b, c, d);
        auto iter = p.cbegin();
        p.emplace(iter, p_0);
        p.push_back(p_1);
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

        d[0] = -Mh[0].x * velocity_u_0(Mh[0].x)*2* Mh[0].x*k[0]*Q[0]/(Mh[1].x+ Mh[0].x);

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
        p.emplace(iter, p[0]-Mh[0].x*(velocity_u_0(Mh[0].x))* 2*(Mh[1].x- Mh[0].x)/ (Mh[1].x + Mh[0].x));
        p.push_back(1);
        return p;
    }
    static std::vector<double> compute_c(std::vector<double> c1,std::vector<long double> u, std::vector<point> Mh,double tau, double D)
    {
        double r_1;
        double r_2;
        for (int i = 1; i < c1.size() - 1; i++)
        {
            r_1 = (Mh[i].x + Mh[i - 1].x)/2;
            r_2 = (Mh[i+1].x + Mh[i].x)/2;
            c1[i] = c1[i] - (tau / m/ Mh[i].x)  * ((c1[i] * u[i]*r_2 - c1[i-1] * u[i - 1] * r_1)/(Mh[i].x- Mh[i-1].x)) + (tau* D / m / Mh[i].x) * (Mh[i + 1].x * (c1[i + 1] - c1[i]) / (pow(Mh[i + 1].x - Mh[i].x, 2)) - Mh[i].x*(c1[i] - c1[i - 1]) / (Mh[i].x - Mh[i - 1].x) / (Mh[i + 1].x - Mh[i].x));
            //c1[i] = c1[i] - (tau / m / Mh[i].x) * ((c1[i + 1] * u[i] * r_2 - c1[i] * u[i - 1] * r_1) / (Mh[i + 1].x - Mh[i].x)) + (tau * D / m ) * ((c1[i + 1] - c1[i]) /Mh[i].x - (c1[i+1] - 2*c1[i - 1]+c1[i-1]) / (Mh[i].x - Mh[i - 1].x) / (Mh[i + 1].x - Mh[i].x));
            c1[c1.size() - 1] = c1[c1.size() - 2];
        }
        return c1;
    }
    std::vector<double> compute_k(std::vector<point> Mh, int i)
    {
        std::vector<double> k1;
        switch (i)
        {
        case(1||4): {
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
                    k1.push_back(0.0440824);
                }
            }
            break; 
        }
        default:
            break;
        }
        return k1;
    }
    double compute_norm(std::vector<long double> p, std::vector<point> Mh, int k) {
        double Norm = 0;
        double exi = 0.0;
        switch (k)
        {
        case(1): {
            for (int i = 0; i < p.size(); i++)
            {
                exi = exact_P_1(Mh[i].x,Mh[0].x);
                Norm += (p[i] - exi) * (p[i] - exi);
            }
            return sqrt(Norm / p.size());
            break;
        }
        case(4): {
            for (int i = 0; i < p.size(); i++)
            {
                exi = exact_P_2(Mh[i].x, Mh[0].x);
                Norm += (p[i] - exi) * (p[i] - exi);
            }
            return sqrt(Norm / p.size());
            break;
        }
        case(2): {
            for (int i = 0; i < p.size(); i++)
            {
                exi = exact_P_22_2(Mh[i].x, Mh[0].x);
                Norm += (p[i] - exi) * (p[i] - exi);
            }
            return sqrt(Norm / p.size());
            break;
        }  
        case(3): {
            for (int i = 0; i < p.size(); i++)
            {
                exi = exact_P_22_3(Mh[i].x, Mh[0].x);
                Norm += (p[i] - exi) * (p[i] - exi);
            }
            return sqrt(Norm / p.size());
            break;
        }
        default:
            break;
        }

    }
    std::vector<double> compout_debit(std::vector<point> Mh, double r, double k, int vr) {
     std::vector<long double> p;
     std::vector<double> k1;
     for (int i = 1; i < Mh.size(); i++)
     {
         if (Mh[i].x <= r)
         {
             k1.push_back(1.0);
         }
         else {
             k1.push_back(k);
         }
     }
     p = compute_usl_1_roda(0,1,Mh, k1);
     std::vector<double> Q;
     for (int i = 1; i < Mh.size(); i++)
     {
         Q.push_back(-2 * PI * Mh[i].x * ( k1[i - 1] * (p[i] - p[i - 1]) / (Mh[i].x - Mh[i - 1].x)));
     }
     return Q;
         }
    void out_file_2_2_5(std::vector<point> Mh,int vr){
        std::ofstream fout("2.2.5.dat");
        fout << "r" << std::endl;
        for (double r:{0.01,0.02,0.05,0.1,0.2,0.5,0.75})
        {
            for (double k : {0.1, 0.2, 0.3,0.4, 0.5, 0.6,0.7,0.8,0.9,1.0})
            {
                fout << r << " "<< k << " _ ";
                for (size_t i = 0; i < Mh.size()-1; i++)
                {
                    fout << compout_debit(Mh, r, k, vr)[i] << " ";
                }
                fout << std::endl;
            }
            fout << std::endl;
        }
        fout.close();
    }
}
using namespace zd2;
void Zadanie_2(const double p_0, const double p_1, std::string mesh) {
    setlocale(LC_ALL, "Russian");
    std::cout << "---------------------------------------Welcome-----------------------------------------" << std::endl;
    std::cout << "---------------------------------------------------------------------------------------" << std::endl;
    static std::vector<point> Mh;
    std::ifstream mesh_(mesh);
    for (int i : {100, 1000, 10000})// i отвечает за количество разбиений
    {
        ///////////////////////////////////////Loookk___sudaaaaa_and_read/////////////////////////////////////////////////////////////////////////
        int zd=2;//1-для задачи 2.1, 2.2(1 сетка),  2-для задачи 2.2 (2 сетка), 3-для задачи 2.2(3 сетка), 4-для задачи 2.1(пункт с граничным условием 2ого рода
        int punkt=1;//1-граничные условия 1ого рода, 2-граничное условие 1-ого и 2-ого рода
        ///////////////////////////////////////Loookk___sudaaaaa_and_read/////////////////////////////////////////////////////////////////////////
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
        std::vector<double> k = compute_k(Mh, zd);
        std::vector<long double> p;
        switch (punkt)
        {
        case(1):
        {
            p = compute_usl_1_roda(0,1,Mh, k);
            break;
        }
        case(2):
        {
            p = compute_usl_1_2_roda(Mh, k);
            break;
        }
        }
        
        double norma=compute_norm(p, Mh,zd);
        std::vector<long double> v = velocyte(p, Mh, k,L);
        out_file(Mh, p, v, k,zd);
        double tau = (Mh[1].x-Mh[0].x)/abs(v[0]);
        std::ofstream fout("C.dat");
        for (int i = 0; i < Mh.size(); i++)
        {
            fout << Mh[i].x << "  ";

        }
        int l = 0;

        std::cout << "Mesh: N=" << i << " " << "Norma:=   " << norma << std::endl;
        std::cout << "---------------------------------------------------------------------------------------" << std::endl;
        //std::cout << exact_P_22_2(0.5, 0.001)*p_1- p_0 << std::endl;
        bool stop = 0;
        bool bre_ak = 0;
        while (!stop)
        {
            std::cout << "---------------------------------------------------------------------------------------" << std::endl;
            std::cout << "#0 (выйти) #1 (график для разных итераций) #2 (график для различных коэффициентов дифузии)"<<std::endl<<
                "#3(время) #4(значение давления на забойне) #5(следующая сетка) #6(задание 2.2.5)" << std::endl;
            std::cout << "---------------------------------------------------------------------------------------" << std::endl;
            std::cin >> l;
            switch (l)
            {
            case 5:
            {
                stop = 1;
                break;
            }
            case 6:
            {
                out_file_2_2_5(Mh,2);
                break;
            }
            case 0:
            {
                stop = 1;
                bre_ak = 1;
                break;
            }
            case 3:
            {
                double t_n = 0.0;
                double t_a = 0.0;
                switch (zd)
                {
                case(1): {
                    t_a = (m * mu * log(Mh[N - 1].x / Mh[0].x) * L * L * (Mh[N - 1].x - Mh[0].x) * (Mh[N - 1].x + Mh[0].x)) / (2 * k_0 * abs(p_0 - p_1));
                    std::cout << "Taim_a:=   " << t_a << " " << "Taim_n:=   " << t_n << std::endl;
                    break;
                }
                case(2): {
                    t_a = (m * mu * log(0.5 / Mh[0].x) * L * L * (0.5 - Mh[0].x) * (0.5 + Mh[0].x)) / (2 * k[0] * k_0 * abs(p_0 - p_1*exact_P_22_2(0.5, Mh[0].x))) + (m * mu * log(Mh[N - 1].x / 0.5) * L * L * (Mh[N - 1].x - 0.5) * (Mh[N - 1].x + 0.5)) / (2 * k[N - 2] * k_0 * abs(p_1 * exact_P_22_2(0.5, Mh[0].x) - p_1));
                    t_n = compute_taim(v, Mh, k) * m * abs(p_0 - p_1) / mu / L;
                    std::cout << "Taim_a:=   " << t_a << " " << "Taim_n:=   " << t_n << std::endl;
                    break;
                }
                case(3): {
                    t_a = (m * mu * log(0.75 / Mh[0].x) * L * L * (0.75 - Mh[0].x) * (0.75 + Mh[0].x)) / (2 * k[0] * k_0 * abs(p_0 - p_1 * exact_P_22_3(0.75, Mh[0].x))) + (m * mu * log(Mh[N - 1].x / 0.75) * L * L * (Mh[N - 1].x - 0.75) * (Mh[N - 1].x + 0.75)) / (2 * k[N - 2] * k_0 * abs(p_1 * exact_P_22_3(0.75, Mh[0].x) - p_1));
                    t_n = compute_taim(v, Mh, k) * m * abs(p_0 - p_1) / mu / L;
                    std::cout << "Taim_a:=   " << t_a << " " << "Taim_n:=   " << t_n << std::endl;
                    break;
                }
                }
                break;  
            }
            case 4:
            {
                std::cout << "p(rw)=" << p[0] << std::endl;
                break;
            }
            case 1: {
                double D = 0.0;
                fout << std::endl;
                std::vector<double> c1(N - 1, 0);
                auto iter_c = c1.cbegin();
                c1.emplace(iter_c, 1);
                for(int k1=1; k1<10;k1++)
                {
                std::vector<double> c1_new = compute_c(c1, velocyte(compute_usl_1_roda(1, 0, Mh, k), Mh, k, L), Mh, tau, D);
                c1 = c1_new;
                    if (k1 % 2 == 0)
                    {
                        double err = 0.0;
                        std::cout << "D__" << D << "__iter__" << k1 << "__taim__" << tau * k1 << std::endl;
                        std::cout << "---------------------------------------------------------------------------------------" << std::endl;

                        for (int i = 0; i < c1.size(); i++)
                        {
                            std::cout << c1_new[i] << "  ";
                            fout << c1_new[i] << "  ";
                        }
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
                    for (int k = 1; k < N_ + 1; k++) {
                        std::vector<double> c1_new = compute_c(c1, v, Mh, tau, D);
                        c1 = c1_new;
                    }
                    double err = 0.0;
                    std::cout << "---------------------------------------------------------------------------------------" << std::endl;
                    std::cout << "D__" << D << "__iter__" << N_ << "__taim__" << tau * N_ << std::endl;
                    std::cout << "---------------------------------------------------------------------------------------" << std::endl;
                    for (int i = 0; i < c1.size(); i++)
                    {
                        std::cout << c1[i] << "  ";
                        fout << c1[i] << "  ";
                    }
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
        Mh.clear();
        fout.close();
        if (bre_ak) break;
    }
}
#include<stdio.h>
#include<math.h>
#include<iostream>
#include<vector>
#include <algorithm>
using namespace std;

const double Mu = 0.01;
const double C = 1;
const double X = 1;
const double T = 1;
const int N = 10;//по t
const int M = 10;// по x
const double t = T / N;
const double h = X / M; 

double f0 (double x, double t);
double f (double x, double t);

class Matrix
{
    public:
        vector<double> _main;
        vector<double> _up;
        vector<double> _down;
        int _len;

        Matrix (int len): _len(len){}

        void three_diag_meth (const vector<double>& f, vector<double>& x)
        {  
            //В начале диагонали _down должен быть 0
            auto it = _down.begin();
            _down.insert(it, 0);
             // _up = B 
            // _main = C
            // _down = A
            //WIKIPEDIA
            vector<double> alpha;
            vector<double> beta;
            x.clear();

            alpha.push_back(0);//номер начинается с 2-х на вики
            beta.push_back(0);//номер начинается с 2-х на вики

            alpha.push_back (-_up[0] / _main[0]);
            beta.push_back (f[0] / _main[0]);

            for (int i = 2; i < _len; i++)
            {
                alpha.push_back (-_up[i - 1] / 
                        (_down[i - 1] * alpha[i - 1] + _main[i - 1]));
                
                beta.push_back ( (f[i - 1] - _down[i - 1] * beta.back()) / 
                        (_down[i - 1] * alpha[i - 1] + _main[i - 1]));
            }

            x.push_back ((f.back() - _down.back() * beta.back()) / 
                        (_main.back() + _down.back() * alpha.back())
            );

            for (int i = _len - 1; i >= 1; i--)
            {
                x.push_back (alpha[i] * x.back() + beta[i]);
            }

            reverse(begin(x), end(x));
            _down.erase(_down.begin());//Удаляем начальный
            

        }
};

double f0 (double x, double t)
{
    double Ro = exp(t) * (cos (M_PI * x / 10.) + 1.5);
    double u = cos(2 * M_PI * t) * sin (M_PI * x * x / 100);
    double res = 0;
    res = Ro;
    res += exp(t) * cos(2 * M_PI * t) *
            (
                -M_PI / 10. * sin(M_PI * x / 10.) * sin(M_PI * x * x / 100.) 
                + 
                 M_PI * x / 50. * cos(M_PI * x * x / 100.) * (cos(M_PI * x / 10.) + 1.5)
            );
    return res;
}


double f (double x, double t)
{
    double Ro = exp(t) * (cos (M_PI * x / 10.) + 1.5);//+

    double u = cos(2 * M_PI * t) * sin (M_PI * x * x / 100);//+

    double du_dt = sin(M_PI * x * x / 100.) * (-2 * M_PI * sin(2 * M_PI * t));//+

    double du_dx = cos(2 * M_PI * t) * cos(M_PI * x * x / 100.) * 2 * M_PI * x / 100.;//+

    double dRo_dx = exp(t) * (-M_PI / 10. * sin(M_PI * x / 10.));//+

    double d2u_dx2 = -M_PI * cos(2 * M_PI * t) / 2500.;//+
    
    d2u_dx2 *= M_PI * x * x * sin(M_PI * x * x / 100.)
               -
               50 * cos(M_PI * x * x / 100);//+
    
    
    double res = Ro * du_dt 
                + Ro * u * du_dx 
                + dRo_dx
                - Mu * d2u_dx2;//+
    res /= Ro;
    return res;

}

double Ro_0 (double x, double t)
{
    double Ro = exp(t) * (cos (M_PI * x / 10.) + 1.5);//+
    return Ro;
}

double u_0(double x, double t)
{
    double u = cos(2 * M_PI * t) * sin (M_PI * x * x / 100);//+
    return u;
}

double count_residual(vector<double> H, vector<double> V, int num_sloy)
{
    return 0.;
}

void solve10Task()
{
    vector<double> H_n;
    vector<double> H_n_1;
    vector<double> V_n;
    vector<double> V_n_1;

    for (int m = 0; m < M + 1; m++)//Начальные условия
    {
        H_n_1.push_back (Ro_0 (h * m, 0));
        V_n_1.push_back (u_0 (h * m, 0));
    }

    for (int i = 1; i < N + 1; i++)
    {
        H_n.clear();
        V_n.clear();

        for (auto i : H_n_1)
            H_n.push_back(i);

        for (auto i : V_n_1)
            V_n.push_back(i);

        Matrix mat(M + 1);
        vector<double> b;

        mat._main.push_back (1. / t 
                            - V_n[0] / (2 * h)
        );

        mat._up.push_back (V_n[1] / (2 * h)); 

        b.push_back(
                    H_n[0] / t
                    - H_n[0] * (V_n[1] - V_n[0]) / (2 * h)
                    + 1/2. * h * (
                                        (H_n[2] * V_n[2] - 2 * H_n[1] * V_n[1] + H_n[0] * V_n[0]) / (h * h)
                                        - 1/2. * (H_n[3] * V_n[3] - 2 * H_n[2] * V_n[2] + H_n[1] * V_n[1]) / (h * h)
                                        + H_n[0] * (
                                                        (V_n[2] - 2 * V_n[1] + V_n[0]) / (h * h)
                                                        - 1/2. * (V_n[3] - 2 * V_n[2] + V_n[1]) / (h * h)
                                                   )
                                   )

        );

        for (int m = 1; m < M; m++)
        {
            mat._main.push_back (1. / t);

            mat._down.push_back (
                                 -V_n[m] / (4 * h)
                                 - V_n[m - 1] / (4 * h)
            );

            mat._up.push_back(
                               V_n[m] / (4 * h)
                               + V_n[m + 1] / (4 * h) 
            );

            b.push_back(
                        H_n[m] / t
                        - H_n[m] * (V_n[m + 1] - V_n[m - 1]) / (4 * h)
            );
        }

        mat._main.push_back(
                            1 / t
                            + V_n[M] / (2 * h)
        );

        mat._down.push_back(
                  V_n[M - 1] / (2 * h)
        );

        b.push_back(
                    H_n[M] / t
                    - 1/2. * H_n[M] * (V_n[M] - V_n[M - 1]) / h
                    + 1/2. * h * (
                                    (H_n[M] * V_n[M] - 2 * H_n[M - 1] * V_n[M - 1] + H_n[M - 2] * V_n[M - 2]) / (h * h)
                                    - 1/2. * (H_n[M - 1] * V_n[M - 1] - 2 * H_n[M - 2] * V_n[M - 2] + H_n[M - 3] * V_n[M - 3]) / (h * h)
                                    + H_n[M] * (
                                                    (V_n[M] - 2 * V_n[M - 1] + V_n[M - 2]) / (h * h)
                                                    - 1/2. * (V_n[M - 1] - 2 * V_n[M - 2] + V_n[M - 3]) / (h * h)
                                               )
                               )
        );
        mat.three_diag_meth(b, H_n_1);//Посчитали значение H на n+1 слое
        b.clear();

        V_n.push_back(1);
        b.push_back(0);
        for (int m = 1; m < M; m++)
        {
            mat._main.push_back (
                                 H_n_1[m] / t
                                 + 2 * Mu / (h * h)
            );

            mat._down.push_back (
                                 -Mu / (h * h)
                                 - H_n_1[m - 1] * V_n[m - 1] / (3 * h)
                                 - H_n_1[m] * V_n[m] / (3 * h)
            );

            mat._up.push_back(
                               -Mu / (h * h)
                               + H_n_1[m + 1] * V_n[m + 1] / (3 * h)
                               + H_n_1[m] * V_n[m] / (3 * h)
            );

            b.push_back(
                        H_n_1[m] * f(h * m, t * (i + 1))
                        + H_n[m] * V_n[m] / t
                        - 1/3. * V_n[m] * V_n[m] * (H_n_1[m + 1] - H_n_1[m - 1]) / (2 * h)
                        - C * (H_n_1[m + 1] - H_n_1[m - 1]) / (2 * h)
            );
        }

        V_n.push_back(1);
        b.push_back(0);
        mat.three_diag_meth(b, V_n_1);//Посчитали значение V на n+1 слое
        count_residual(H_n_1, V_n_1, i);
    }
    
}


int main()
{
    vector<double> main, up, down, b, x;
    Matrix a(6);
    a._up = {0,0,0,0,0};
    a._down = {0,0,0,0,0};
    a._main = {1,1,1,12,4,5};
    b = {1,1,1,12,4,5};
    
    
    a.three_diag_meth(b, x);
    
    for (auto i: x)
        cout << i << "\n";  
    
    
    return 0;
}
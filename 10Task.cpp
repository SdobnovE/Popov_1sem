#include<stdio.h>
#include<math.h>
#include<iostream>
#include<vector>
#include <algorithm>
#include <utility>
using namespace std;
typedef  pair<double, double> residual;
residual final_residual;
double Mu = 0.01;
double C = 1.0;
const double X = 10;
const double T = 1;
string norm;
int N = 400;//по t
int M = 400;// по x
double t = T / N;
double h = X / M; 
const double EPS = 1e-16;

double f0 (double x, double t);
double f (double x, double t);
double Ro_0 (double x, double t);
double u_0(double x, double t);
residual count_residual(const vector<double>& H, const vector<double>& V);
void solve10Task();


class Matrix
{
    public:
        vector<double> _main;
        vector<double> _up;
        vector<double> _down;
        int _len;

        Matrix (int len): _len(len){}
        void print()
        {
            printf("mid:\n");
            for (int i = 0; i < _len; i++)
                printf ("%e ", _main[i]);
            printf("\n");

            printf("up:\n");
            for (int i = 0; i < _len - 1; i++)
                printf ("%e ", _up[i]);
            printf("\n");

            printf("down:\n");
            for (int i = 0; i < _len - 1; i++)
                printf ("%e ", _down[i]);
            printf("\n");
            
        }

        double residual(const vector<double>& f, vector<double>& x)
        {
            vector<double> res;
            this->prod_vec(x, res);
            double resu = 0;
            for (int i = 0; i < _len; i++)
                resu += (res[i] - f[i]) * (res[i] - f[i]);
            return sqrt (resu);

        }

        void prod_vec (const vector<double>& f, vector<double>& x)
        {
            x.clear();
            x.push_back (_main[0] * f[0] + _up[0] * f[1]);
            
            for (int i = 1; i < _len; i++)
            {
                x.push_back (_main[i] * f[i] + _up[i] * f[i + 1] + _down[i - 1] * f[i - 1]);
            }
            x.push_back(_main.back() * f[_len - 1] + _down.back() * f[_len - 2]);
        }

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

double f0 (int i, int j)
{
    double x = j * h;
    double Ro = exp(t * i) * (cos (M_PI * x / 10.) + 1.5);
    
    double res = 0;
    res = Ro;
    res += exp(t * i) * cos(2 * M_PI * t * i) *
            (
                -M_PI / 10. * sin(M_PI * x / 10.) * sin(M_PI * x * x / 100.) 
                + 
                 M_PI * x / 50. * cos(M_PI * x * x / 100.) * (cos(M_PI * x / 10.) + 1.5)
            );
    return res;
}


double f (int i, int j)
{
    
    return (
        -exp(i * t) 
            * (cos(M_PI * h * j / 10.) + 1.5) 
            * 2. * M_PI * sin(2. * M_PI * i * t) 
            * sin(M_PI * j * j * h * h / 100.)
        + exp(i * t) 
            * (cos(M_PI * j * h / 10.) + 1.5) 
            * sin(M_PI * j * j * h * h / 100.) 
            * cos(2. * M_PI * i * t) 
            * cos(2. * M_PI * i * t) 
            * (M_PI * h  * j / 50.) 
            * cos(M_PI * j * j * h * h / 100.)
        - C * (M_PI / 10.) 
            * exp(i * t)
            * sin(M_PI * h * j / 10.)
        - Mu * cos(2. * M_PI * i * t)
            * (
                (M_PI / 50.) * cos(M_PI * j * h * j * h / 100.)
                - (M_PI *M_PI * j * j * h * h / 2500)
                    * sin(M_PI * j * h * j * h / 100.)
               )
            ) / (
                    exp(i * t) 
                    * (cos(M_PI * j * h / 10.) + 1.5)
                );

}

double Ro_0 (double x, double t)
{
    double Ro = exp(t) * (cos (M_PI * x / 10.) + 1.5);//+
    return Ro;
}

double u_0 (double x, double t)
{
    double u = cos(2 * M_PI * t) * sin (M_PI * x * x / 100);//+
    return u;
}

residual count_residual(const vector<double>& H, const vector<double>& V)
{
    residual resid;
    resid.first = 0;
    resid.second = 0;
    
   if (norm == "C")
    {
        for (int i = 0; i < M + 1; i++)
        {
            if (fabs (Ro_0(h * i, 1) - H[i]) > resid.first)
                resid.first = fabs (Ro_0(h * i, 1) - H[i]);
        }
        
        for (int i = 0; i < M + 1; i++)
        {
            if (fabs (u_0 (h * i, 1) - V[i]) > resid.second)
                resid.second = fabs (u_0 (h * i, 1) - V[i]);
        }
    }

    if (norm == "L2")
    {
        double scal_V = 0.0;
        double scal_H = 0.0;

        


        for (int i = 1; i < M; i++)//1...N
            scal_H += (Ro_0(h * i, 1) - H[i]) * (Ro_0(h * i, 1) - H[i]);
            
        for (int i = 1; i < M; i++)//1...N
            scal_V += (u_0 (h * i, 1) - V[i]) * (u_0 (h * i, 1) - V[i]);
        
        scal_H *= h;
        scal_V *= h;
        
        scal_H += 0.5 * h * ((Ro_0(0, 1) - H[0]) * (Ro_0(0, 1) - H[0])
                             + (Ro_0(h * M, 1) - H[M]) * (Ro_0(h * M, 1) - H[M])
                            );
        scal_V += 0.5 * h * ((u_0 (h * 0, 1) - V[0]) * (u_0 (h * 0, 1) - V[0])
                             + (u_0 (h * M, 1) - V[M]) * (u_0 (h * M, 1) - V[M]));

        resid.first = sqrt (scal_H);
        resid.second = sqrt (scal_V);
    }

    if (norm == "W")
    {
        double scal_V = 0.0;
        double scal_H = 0.0;

        


        for (int i = 1; i < M; i++)//1...N
            scal_H += (Ro_0(h * i, 1) - H[i]) * (Ro_0(h * i, 1) - H[i]);
            
        for (int i = 1; i < M; i++)//1...N
            scal_V += (u_0 (h * i, 1) - V[i]) * (u_0 (h * i, 1) - V[i]);
        
        scal_H *= h;
        scal_V *= h;
        
        scal_H += 0.5 * h * ((Ro_0(0, 1) - H[0]) * (Ro_0(0, 1) - H[0])
                             + (Ro_0(h * M, 1) - H[M]) * (Ro_0(h * M, 1) - H[M])
                            );
        scal_V += 0.5 * h * ((u_0 (h * 0, 1) - V[0]) * (u_0 (h * 0, 1) - V[0])
                             + (u_0 (h * M, 1) - V[M]) * (u_0 (h * M, 1) - V[M]));

        resid.first = scal_H;
        resid.second = scal_V;
        
        double temp1 = 0.;
        for (int i = 0; i < M; i++)
        {
            temp1 += (u_0 (h * (i), 1) - V[i] - u_0 (h * (i + 1), 1) + V[i + 1]) / h * 
                     (u_0 (h * (i), 1) - V[i] - u_0 (h * (i + 1), 1) + V[i + 1]) / h;
        }
        temp1 *= h;

        double temp2 = 0.;
        for (int i = 0; i < M; i++)
        {
            temp2 += (Ro_0 (h * (i), 1) - H[i] + H[i + 1] - Ro_0 (h * (i + 1), 1)) / h * 
                     (Ro_0 (h * (i), 1) - H[i] + H[i + 1] - Ro_0 (h * (i + 1), 1)) / h;
        }
        temp2 *= h;
        resid.first = sqrt(fabs (resid.first + temp2));
        resid.second = sqrt(fabs (resid.second + temp1));
        


    }

    
    
    return resid;
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
    }//Ok

    

    for (int i = 1; i < N + 1; i++)
    {
        H_n.clear();
        V_n.clear();
        for (auto i1 : H_n_1)
            H_n.push_back(i1);

        for (auto i1 : V_n_1)
            V_n.push_back(i1);

        Matrix mat(M + 1);
        vector<double> b;

        mat._main.push_back (
                            1. / t 
                            - V_n[0] / (2 * h)
        );
        
        mat._up.push_back (V_n[1] / (2 * h)); 

        b.push_back(
                    H_n[0] / t
                    - H_n[0] * (V_n[1] - V_n[0]) / (2.0 * h)
                    + (1 / (2.0 * h)) * (
                                            H_n[2] * V_n[2] - 2 * H_n[1] * V_n[1] + H_n[0] * V_n[0]
                                            - 0.5 * (H_n[3] * V_n[3] - 2 * H_n[2] * V_n[2] + H_n[1] * V_n[1])
                                            + H_n[0] * (
                                                        V_n[2] - 2 * V_n[1] + V_n[0] 
                                                        - 0.5 * (V_n[3] - 2 * V_n[2] + V_n[1])
                                                    )
                                        )
                    
                    + f0(i, 0)
        );


        

        for (int m = 1; m < M; m++)////Сто проц верный цикл
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
                        + f0(i, m)
            );
        }

        mat._main.push_back(
                            1 / t
                            + V_n[M] / (2. * h)
        );

        mat._down.push_back(
                  -V_n[M - 1] / (2. * h)
        );//Возможно минус

        b.push_back(
                    H_n[M] / t
                    - 1/2. * H_n[M] * (V_n[M] - V_n[M - 1]) / h
                    - 1/2. * h * (
                                    (H_n[M] * V_n[M] - 2 * H_n[M - 1] * V_n[M - 1] + H_n[M - 2] * V_n[M - 2]) / (h * h)
                                    - 1/2. * (H_n[M - 1] * V_n[M - 1] - 2 * H_n[M - 2] * V_n[M - 2] + H_n[M - 3] * V_n[M - 3]) / (h * h)
                                    + H_n[M] * (
                                                    (V_n[M] - 2 * V_n[M - 1] + V_n[M - 2]) / (h * h)
                                                    - 1/2. * (V_n[M - 1] - 2 * V_n[M - 2] + V_n[M - 3]) / (h * h)
                                               )
                               )
                    + f0(i, M)
        );
        

        mat.three_diag_meth(b, H_n_1);//Посчитали значение H на n+1 слое
        //printf("res %e\n", mat.residual(b, H_n_1));
        mat._main.clear();
        mat._up.clear();
        mat._down.clear();
        b.clear();

        mat._main.push_back(1);
        mat._up.push_back(0);
        b.push_back(0);

       double mm = 1. / H_n_1[0];
        
        for (int j = 1; j <= M; j++)
        {
            if (1. / H_n_1[j] > mm)
            {
                mm = 1. / H_n_1[j];
            }
        }
        mm = Mu * mm;

        for (int j = 1; j < M; j++)
        {
            mat._down.push_back(
                -(V_n[j] + V_n[j - 1]) / (6 * h)
                - mm / (h * h)
            );

            mat._main.push_back(
                1. / t 
                + 2. * mm / (h * h)
            );
            
            mat._up.push_back(
                (V_n[j] + V_n[j + 1]) / (6*h)
                - mm / (h * h)
            );

            b.push_back(
                V_n[j] / t
                - C * (H_n_1[j + 1] - H_n_1[j - 1]) / (2. * h * H_n_1[j])
                - (mm - Mu / H_n_1[j]) * (V_n[j + 1] - 2 * V_n[j] + V_n[j - 1]) / (h * h)
                + f(i, j)
            );
            
        }

        mat._down.push_back(0);
        mat._main.push_back(1);
        b.push_back(0);

        mat.three_diag_meth(b, V_n_1);//Посчитали значение V на n+1 слое    
        
        //count_residual(H_n_1, V_n_1, i);//Посчитали невязку
        
        
    }
    final_residual = count_residual(H_n_1, V_n_1);
    
}


int main(int argc, char* argv[])
{
    if (argc != 6)
    {
        printf("USAGE: Mu, C, N_t, M_x, norm\n");
        return -1;
    }

    sscanf(argv[1], "%lf", &Mu);
    sscanf(argv[2], "%lf", &C);
    sscanf(argv[3], "%d", &N);
    //printf("argv[3] %s\n", argv[3]);
    sscanf(argv[4], "%d", &M);
    
    norm = argv[5];
    //printf("HA %f %f %d %d\n", Mu, C, N, M);
    t = T / N;
    h = X / M; 

    solve10Task();
    vector<double> b;
    vector<double> a;
    printf ("%e %e\n", final_residual.first, final_residual.second);

    
    
    return 0;
}
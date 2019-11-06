#include<stdio.h>
#include<iostream>
#include<vector>
#include <algorithm>
using std::vector;
using std::cout;

class Matrix
{
    public:
        vector<double> _main;
        vector<double> _up;
        vector<double> _down;
        int64_t _len;

        Matrix (uint64_t len): _len(len){}

        void three_diag_meth (const vector<double>& f, vector<double>& x) const
        {  
            //В начале диагонали _down должен быть 0
             // _up = B 
            // _main = C
            // _down = A
            //WIKIPEDIA
            vector<double> alpha;
            vector<double> beta;
            x.clear();

            alpha.push_back(0);
            beta.push_back(0);

            alpha.push_back (-_up[0] / _main[0]);
            beta.push_back (f[0] / _main[0]);

            for (auto i = 2; i < _len; i++)
            {
                alpha.push_back (-_up[i - 1] / 
                        (_down[i - 1] * alpha[i - 1] + _main[i - 1]));
                
                beta.push_back ( (f[i - 1] - _down[i - 1] * beta.back()) / 
                        (_down[i - 1] * alpha[i - 1] + _main[i - 1]));
            }

            x.push_back ((f.back() - _down.back() * beta.back()) / 
                        (_main.back() + _down.back() * alpha.back())
            );

            cout <<"\n\tHERE:    " << x.back() << "\n\n";
            cout << _len << " " << alpha.size() << "\n";

            for (int64_t i = _len - 1; i >= 1; i--)
            {
                x.push_back (alpha[i] * x.back() + beta[i]);
            }

            std::reverse(std::begin(x), std::end(x));
            return;

        }
};



int main()
{
    vector<double> main, up, down, b, x;
    Matrix a(6);
    a._up = {2,3,2,4,2};
    a._down = {0,-2,-2,12,3,4};
    a._main = {1,1,1,12,4,5};
    b = {1,2,33,24,4,7};
    // a._up = {2};
    // a._down = {0,2};
    // a._main = {10,10};
    // b = {1,1};
    
    a.three_diag_meth(b, x);
    
    for (auto i: x)
        cout << i << "\n";
    return 0;
}
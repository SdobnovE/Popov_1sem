#ifndef MATRIX_3DIAGONAL_H
#define MATRIX_3DIAGONAL_H
#include <math.h>
#include <vector>

#define eps 1e-20
#define EQ(x) (fabs(x) < eps)


class matrix_3diagonal
{
private:
  std::vector<double> RD;
  std::vector<double> MD;
  std::vector<double> LD;
  int n;
  std::vector<double> rhs;
public:
  matrix_3diagonal (int n_)
  {
    n = n_;
    resize ();
  }

  void fill_matrix_H (double tau, double h, std::vector<double> &V, std::vector<double> &H, std::vector<double> &debug_vector);
  void fill_matrix_V (double gamma, double tau, double mu, double h, std::vector<double> &V, std::vector<double> &H, std::vector<double> &H_old, std::vector<double> &debug_vector);
  void solve (std::vector<double> &ans);
  void resize ()
  {
    RD.resize (n - 1);
    MD.resize (n);
    LD.resize (n);
    rhs.resize (n);
  }
  ~matrix_3diagonal()
  {
    RD.clear ();
    MD.clear ();
    LD.clear ();
    rhs.clear ();
  }
};

#endif // MATRIX_3DIAGONAL_H

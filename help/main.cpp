#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "matrix_3diagonal.h"
#include "functions.h"
#include <string.h>
#include <string>

#define NMAX 11

class data
{
public:
  double h;
  double tau;
  int m;
  int n;
  double x0, x1;
  double t0, t1;

  double exponent;
  double viscosity;
  data ()
  {
    m = 1000;
    n = 1000;
    x0 = 0;
    x1 = 10;
    t0 = 0;
    t1 = 10;
    h = (x1 - x0) / m;
    tau = (t1 - t0) / n;
    exponent = 1.4;
    viscosity = 0.1;
  }
  void print_properties ()
  {
    printf ("All properties will be defaulted\nX = [%f, %f]\nT = [%f, %f]\nn = %d, m = %d\nh = %f, tau = %f\nexponent (gamma) = %f\nviscosity (mu) = %f\n",
            x0, x1,
            t0, t1,
            n, m,
            h, tau,
            exponent,
            viscosity);
  }

  void compute_residual_v (std::vector<double> &v, double t, FILE * fp = stdout)
  {
    int k = v.size ();
    double s = 0;
    double max = 0;
    double res = 0;

    for (int i = 0; i < k; i++)
      {
        s += (velocity(x0 + i * h, t) - v[i]) * (velocity(x0 + i * h, t) - v[i]);
//        printf ("{%f} %f \t %f\n", x0 + i * h, v[i], velocity(x0 + i * h, t));

        if (i < k - 1)
          res += (velocity(x0 + i * h, t) - v[i]) * (velocity(x0 + i * h, t) - v[i]);

        if (fabs (v[i] - velocity(x0 + i * h, t)) > max)
          max = fabs (v[i] - velocity(x0 + i * h, t));
      }
    s *= h;

    s += h / 2 * ((velocity(x0, t) - v[0]) * (velocity(x0, t) - v[0]) + (velocity(x1, t) - v[k - 1]) * (velocity(x1, t) - v[k - 1]));
    res *= h;

    fprintf (fp, "V L2h-residual = %e\n", sqrt (s));

    fprintf (fp, "V Ch-residual  = %e\n", max);

    fprintf (fp, "V 2-residual   = %e\n", sqrt (s + res));
  }
  void compute_residual_h (std::vector<double> &v, double t, FILE * fp = stdout)
  {
    int k = v.size ();
    double s = 0;
    double max = 0;
    double res = 0;

    for (int i = 0; i < k; i++)
      {
        s += (rho(x0 + i * h + h / 2, t) - v[i]) * (rho(x0 + i * h + h / 2, t) - v[i]);
//        printf ("{%f} %f \t %f\n", x0 + i * h + h / 2, v[i], rho(x0 + i * h + h / 2, t));

        if (i < k - 1)
          res += (rho(x0 + i * h + h / 2, t) - v[i]) * (rho(x0 + i * h + h / 2, t) - v[i]);

        if (fabs (v[i] - rho(x0 + i * h + h / 2, t)) > max)
          max = fabs (v[i] - rho(x0 + i * h + h / 2, t));
      }
    s *= h;

    s += h / 2 * ((rho(x0 + h / 2, t) - v[0]) * (rho(x0 + h / 2, t) - v[0]) + (rho(x1 - h / 2, t) - v[k - 1]) * (rho(x1 - h / 2, t) - v[k - 1]));

    fprintf (fp, "H L2h-residual = %e\n", sqrt (s));

    fprintf (fp, "H Ch-residual  = %e\n", max);

    fprintf (fp, "H 2-residual   = %e\n", sqrt (s + res));
  }
};

void fill_h0_v0_f0 (
    data &dat,
    std::vector<double> &h,
    std::vector<double> &v,
    double (*v0)(double),
    double (*h0)(double))
{
  int m = dat.m;
  double left_edge = dat.x0;
  double step = dat.h;
  h.clear ();
  v.clear ();

#if 0
  for (int i = 0; i <= m; i++)
    {
      if (i < m)
        h.push_back (h0(left_edge + i * step + step / 2));

      v.push_back (v0(left_edge + i * step));
    }
#else

  (void) h0;
  (void) v0;
  double x = left_edge;
  for (int i = 0; i <= m; i++)
    {
      if (i < m)
        {
          if (x + step / 2 < 4.5 || x + step / 2 > 5.5)
            {
              h.push_back(1);
            }
          else
            {
              h.push_back(2);
            }
        }
      if (x < 4.5 || x > 5.5)
        {
          v.push_back(0);
        }
      else
        {
          v.push_back(0);
        }
      x += step;
    }
#endif
}
void print_vector (std::vector<double> &v, std::string a)
{
  int n = v.size ();

  int m = n > NMAX ? NMAX : n;
  printf ("%s", a.c_str ());
  for (int i = 0; i < m; i++)
    {
      printf ("[%d] %f\n", i, v[i]);
    }
  if (n > m)
    printf ("...\n");
  printf ("\n");
}
void fill_rhf (double t, data &dat,std::vector<double> &f, std::vector<double> &g, double mu, double gamma)
{
  int m = dat.m;
  f.resize (m);
  g.resize (m + 1);
  for (int i = 0; i <= m; i++)
    {
      if (i < m)
        f[i] = f1 (dat.x0 + i * dat.h + dat.h / 2, t);
      g[i] = f2 (dat.x0 + i * dat.h, t, mu, gamma);
    }
}

double compute_mas (data &dat, std::vector<double> &rh)
{
  double s = 0;
  for (int i = 0; i < dat.m; i++)
    {
      s += rh[i];
    }
  return s * dat.h;
}

double norm (std::vector<double> &v, double h)
{
    double sum = 0;

    for (unsigned int i = 1; i < v.size () - 1; i++)
      {
        if (EQ(v[i]))
          v[i] = 0;
        sum += v[i] * v[i];
      }
    sum += (v.front() * v.front() + v.back() * v.back()) / 2;
    sum *= h;
    return sqrt (sum);
}

double conts (std::vector<double> &h, double d)
{
  double min = 1e+300;
  double max = 0;

  for (unsigned int i = 0; i < h.size (); i++)
    {
      min = h[i] < min ? h[i] : min;
      max = h[i] > max ? h[i] : max;
    }
  return d * (max - min);
}

int main()
{
  data dat;

  FILE * fh = fopen ("./H.txt", "w");
  FILE * fv = fopen ("./V.txt", "w");
  FILE *res = fopen ("./res.txt", "w");

  dat.print_properties ();

  std::vector<std::string> outH, outV;

  matrix_3diagonal H (dat.m);
  matrix_3diagonal V (dat.m - 1);

  std::vector<double> h;
  std::vector<double> h_old;
  std::vector<double> v;

  std::vector<double> rhf1;
  std::vector<double> rhf2;

  outH.resize (dat.m);
  outV.resize(dat.m + 1);
  for (int i = 0; i <= dat.m; i++)
    {
      if (i < dat.m)
        {
          char buf[100];
          sprintf (buf, "%f ", dat.x0 + i * dat.h + dat.h / 2);
          outH[i] = buf;
        }
      char buf[100];
      sprintf (buf, "%f ", dat.x0 + i * dat.h);
      outV[i] = buf;
    }

  fill_h0_v0_f0 (dat, h, v, v_init, h_init);

  double mas = compute_mas (dat, h);

  printf ("MASSA = %f\n", mas);

  fill_rhf (dat.t0, dat, rhf1, rhf2, dat.viscosity, dat.exponent);
  h_old = h;
  print_vector (h, "H init \n");
  print_vector (v, "V init\n");

  dat.compute_residual_v (v, dat.t0);
  dat.compute_residual_h (h, dat.t0);

  for (int i = 0; i <= dat.m; i++)
    {
      if (i < dat.m)
        {
          char buf[100];
          sprintf (buf, "%f ", h[i]);
          outH[i] += buf;
        }
      char buf[100];
      sprintf (buf, "%f ", v[i]);
      outV[i] += buf;
    }

  printf ("==========step 0===========\n");
  printf ("%f,\t %f\n", (norm (v, 1)), conts (h, 1));

  H.fill_matrix_H (dat.tau, dat.h, v, h, rhf1);
  H.solve (h);

  print_vector (h, "H0 \n");
  dat.compute_residual_h (h, dat.t0 + dat.tau);

  V.fill_matrix_V (dat.exponent, dat.tau, dat.viscosity, dat.h, v, h, h_old, rhf2);
  V.solve (v);

  v.insert (v.begin (), 0);
  v.push_back (0);

  print_vector (v, "V0\n");
  dat.compute_residual_v (v, dat.t0 + dat.tau);

  for (int i = 0; i <= dat.m; i++)
    {
      if (i < dat.m)
        {
          char buf[100];
          sprintf (buf, "%f ", h[i]);
          outH[i] += buf;
        }
      char buf[100];
      sprintf (buf, "%f ", v[i]);
      outV[i] += buf;
    }

  for (int i = 1; !(((norm (v, 1)) < 1e-2) && conts (h, 1) < 1e-2) && i < 1000000; i++)
//  for (int i = 1; i < 2; i++)
    {
//      EQ(norm(h, dat.h));
      printf ("================ITER %d===========\n", i);
      printf ("%f,\t %f\n", (norm (v, 1)), conts (h, 1));
      h_old = h;
      fill_rhf (dat.t0 + i * dat.tau, dat, rhf1, rhf2, dat.viscosity, dat.exponent);
      H.fill_matrix_H (dat.tau, dat.h, v, h, rhf1);
      H.solve (h);

      V.fill_matrix_V (dat.exponent, dat.tau, dat.viscosity, dat.h, v, h, h_old, rhf2);
      V.solve (v);

      v.insert (v.begin (), 0);
      v.push_back (0);

//      for (int q = 0; q <= dat.m; q++)
//        {
//          if (q < dat.m)
//            {
//              char buf[100];
//              sprintf (buf, "%f ", h[q]);
//              outH[q] += buf;
//            }
//          char buf[100];
//          sprintf (buf, "%f ", v[q]);
//          outV[q] += buf;
//        }

//      print_vector (h, "H\n");
//      dat.compute_residual_h (h, dat.t0 + (i + 1) * dat.tau);
//      print_vector (v, "V\n");
//      dat.compute_residual_v (v, dat.t0 + (i + 1) * dat.tau);
    }
//  print_vector (h, "H\n");
//  dat.compute_residual_h (h, dat.t1, res);
//  print_vector (v, "V\n");
//  dat.compute_residual_v (v, dat.t1, res);
//  print_vector (h, "H\n");
//  dat.compute_residual_h (h, dat.t1);
//  print_vector (v, "V\n");
//  dat.compute_residual_v (v, dat.t1);

//  for (int i = 0;  i < dat.m; i++)
//    {
//      fprintf (fh, "%s\n", outH[i].c_str ());
//    }
//  for (int i = 0; i <= dat.m; i++)
//    {
//      fprintf (fv, "%s\n", outV[i].c_str ());
//    }

  printf ("Mat balance = %e\n", fabs (mas - compute_mas (dat, h)));

  fclose (fh);
  fclose (fv);
  fclose (res);

  return 1;
}


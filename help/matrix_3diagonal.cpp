#include "matrix_3diagonal.h"

#define SIGN(X) (EQ(X) ? 0 : X < 0 ? -1 : 1)

void matrix_3diagonal::fill_matrix_H (double tau, double h, std::vector<double> &V, std::vector<double> &H, std::vector<double> &debug_vector)
{
#if 1
  for (int i = 0; i < n; i++)
    {
      if (i < n - 1)
        {
          RD[i] = tau / h / 2 * (V[i + 1] - fabs (V[i + 1]));
        }
      if (i > 0)
        LD[i - 1] = - tau / h / 2 * (V[i] + fabs (V[i]));

      MD[i] = 1 + tau / h / 2 * (V[i + 1] + fabs (V[i + 1]) - V[i] + fabs (V[i]));
      rhs[i] = H[i] + tau * debug_vector[i];
    }
#else
  for (int i = 0; i < n; i++)
    {
      MD[i] = 1;
      rhs[i] = H[i] + tau * debug_vector[i];
      if (V[i] < 0)
        {
          if (i < n - 1)
            RD[i]  = tau / h * V[i + 1];

          MD[i] -= tau / h * V[i];

          if (i > 0)
            LD[i - 1] = 0;
        }
      else
        {
          if (i < n - 1)
            RD[i]  = 0;

          MD[i] += tau / h * V[i + 1];

          if (i > 0)
            LD[i - 1] = -tau / h * V[i];
        }
    }
#endif
}
void matrix_3diagonal::fill_matrix_V (double gamma, double tau, double mu, double h, std::vector<double> &V, std::vector<double> &H, std::vector<double> &H_old, std::vector<double> &debug_vector)
{
#if 1
  for (int i = 1, k = 0; k < n; i++, k++)
    {
      if (EQ(H[i] + H[i-1]))
        {
          LD[k] = 0;
          MD[k] = 1;
          RD[k] = 0;
          rhs[k] = 0;
        }
      else
        {
          if (k < n - 1)
            {
              RD[k] = -tau / h / 4 * ((fabs(V[i]) - V[i]) * H[i] + (fabs (V[i + 1]) - V[i + 1]) * H[i + 1]) - mu * tau / h / h;
            }
          MD[k] = (H[i - 1] + H[i]) / 2 + tau / h / 4 * ((fabs (V[i - 1]) - V[i - 1] + fabs(V[i]) + V[i]) * H[i - 1] + (fabs (V[i + 1]) + V[i + 1] + fabs (V[i]) - V[i]) * H[i]) + mu * 2 * tau / h / h;

          if (k > 0)
            {
              LD[k - 1] = - tau / h / 4 * ((fabs (V[i - 1]) + V[i - 1]) * H[i - 2] + (fabs (V[i]) + V[i]) * H[i - 1]) - mu * tau / h / h;
            }
          rhs[k] = (H_old[i - 1] + H_old[i]) * V[i] / 2 - gamma / (gamma - 1) * tau / h / 2 * (H[i] + H[i - 1]) * (pow((H[i]), gamma - 1) - pow((H[i - 1]), gamma - 1)) + tau * debug_vector[i];
        }
    }
#else
  for (int i = 1, k = 0; k < n; i++, k++)
    {
      if (EQ(H[i] + H[i-1]))
        {
          if (k > 0)
            LD[k - 1] = 0;
          MD[k] = 1;
          RD[k] = 0;
          rhs[k] = 0;
        }
      else
        {
          if (V[i] < 0)
            {
              if (k < n - 1)
                {
                  RD[k] = tau / h / 2 * (H[i] * V[i] + H[i + 1] * V[i + 1]) - mu * tau / h / h;
                }

              MD[k] = (H[i] + H[i - 1]) / 2 - tau / h / 2 * (H[i - 1] * V[i - 1] + H[i] * V[i]) + mu * tau / h / h * 2;

              if (k > 0)
                LD[k - 1] = - mu * tau / h / h;

              rhs[k] = (H_old[i] + H_old[i - 1]) * V[i] / 2 - gamma / (gamma - 1) * tau / h / 2 * (H[i] + H[i - 1]) * (pow(H[i], gamma - 1) - pow(H[i - 1], gamma - 1)) + tau * debug_vector[i];
            }
          else
            {
              if (k < n - 1)
                RD[k] = - mu * tau / h / h;

              MD[k] = (H[i] + H[i - 1]) / 2 + tau / h / 2 * (H[i - 1] * V[i] + H[i] * V[i + 1]) + mu * tau / h / h * 2;

              if (k > 0)
                LD[k - 1] = -tau / h / 2 * (H[i - 2] * V[i - 1] + H[i - 1] * V[i]) - mu * tau / h / h;

              rhs[k] = (H_old[i] + H_old[i - 1]) * V[i] / 2 - gamma / (gamma - 1) * tau / h / 2 * (H[i] + H[i - 1]) * (pow(H[i], gamma - 1) - pow(H[i - 1], gamma - 1)) + tau * debug_vector[i];
            }
        }
    }
#endif
}
void matrix_3diagonal::solve (std::vector<double> &ans)
{
    std::vector<double> a;
    std::vector<double> b;
    a.resize(n);
    b.resize(n);
    ans.resize(n);
    double y = MD[0];

    a[0] = -RD[0] / y;
    b[0] = rhs[0] / y;

    for (int  i = 1; i < n - 1; i++)
    {
        y = MD[i] + LD[i - 1] * a[i - 1];
        a[i] = - RD[i] / y;
        b[i] = (rhs[i] - LD[i - 1] * b[i - 1]) / y;
    }

    ans[n - 1] = (rhs[n - 1] - LD[n - 2] * b[n - 1]) / (MD[n - 1] + LD[n - 2] * a[n - 1]);

    for (int i = n - 2; i >= 0; i--)
    {
        ans[i] = a[i] * ans[i + 1] + b[i];
    }
}

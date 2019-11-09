#include "functions.h"

double rho (double x, double t)
{
  return exp (t) * (x + 1);
  return exp(t) * (cos(M_PI * x / 10) + 1.5);
}
double velocity (double x, double t)
{
  return x * (x - 1);
  return cos (2 * M_PI * t) * sin (M_PI * x * x / 10 / 10);
}

double h_init (double x)
{
  (void) x;
  return rho(x, 0);
}

double v_init (double x)
{
  (void) x;
  return velocity (x, 0);
}

// f1 = d(rho)/dt + d(rho * u) / dx
double f1 (double x, double t)
{
  return  0;
  return exp (t) *(x + 1) + exp (t) * (3 * x * x - 1);
  return rho (x, t) - velocity(x, t) * M_PI / 10 * exp(t) * sin(M_PI * x / 10) + rho(x, t) * cos(2 * M_PI * t) * cos (M_PI * x * x / 100) * M_PI / 50 * x;
}
//f2 = d(rho * u) / dt + d(rho * u * u) / dx + dp/dx - mu d(du/dx)/dx
//p(x) = rho(x) ^ gamma
double f2 (double x, double t, double mu, double gamma)
{
  return 0;
  return exp (t) * x * (x * x - 1) + exp (t) * (2 * x * (x * x - 1) * (x - 1) + 2 * x * x * x * (x - 1) + x * x * (x * x - 1)) + gamma * exp (t) * pow (rho (x, t), gamma - 1) - 2 * mu;
  return velocity (x, t) * rho(x, t) - rho(x, t) * 2 * M_PI * sin (2 * M_PI * t) * sin (M_PI * x * x / 100)
      - M_PI / 10 * exp(t) * sin(M_PI * x / 10) * velocity(x, t) * velocity(x, t) + rho(x, t) * 2 * velocity(x, t) * cos(2 * M_PI * t) * cos(M_PI * x * x / 100) * 2 * M_PI * x / 100
      + gamma * pow (rho(x, t), gamma - 1) * (- exp(t) * sin(M_PI * x / 10) * M_PI / 10)
      - mu * cos (2 * M_PI * t) * (-sin(M_PI * x * x / 100) * M_PI * M_PI * x * x / 50 / 50 + cos(M_PI * x * x / 100) * 2 * M_PI / 100);
}

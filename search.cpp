#include <math.h>
#include <stdio.h>

bool lightnoise(const float qrms, const float mqtq, const float rmsts,
                const float qdiff)
{
  if(mqtq > 0.12 || mqtq < 0) return true;
  if(qdiff > 30e3) return true;
  if(rmsts >= 36 && qrms >= 464 - 8*rmsts) return true;
  return false;
}

__attribute__((unused)) static double
pol6(const double x, const double p0, const double p1, const double p2,
     const double p3, const double p4, const double p5, const double p6)
{
  return p0
       + p1 * x
       + p2 * x*x
       + p3 * x*x*x
       + p4 * x*x*x*x
       + p5 * x*x*x*x*x
       + p6 * x*x*x*x*x*x;
}

__attribute__((unused)) static double
pol5(const double x, const double p0, const double p1, const double p2,
     const double p3, const double p4, const double p5)
{
  return p0
       + p1 * x
       + p2 * x*x
       + p3 * x*x*x
       + p4 * x*x*x*x
       + p5 * x*x*x*x*x;
}


__attribute__((unused)) static double
pol4(const double x, const double p0, const double p1, const double p2,
     const double p3, const double p4)
{
  return p0
       + p1 * x
       + p2 * x*x
       + p3 * x*x*x
       + p4 * x*x*x*x;
}

__attribute__((unused)) static double
pol3(const double x, const double p0, const double p1, const double p2,
     const double p3)
{
  return p0
       + p1 * x
       + p2 * x*x
       + p3 * x*x*x;
}

double fidocorrx(const double x)
{
  return pol5(x, -1.43685, 1-0.102688, -3.53047e-06, 1.35172e-07,
                 -2.12269e-14, -4.2652e-14);
}

double fidocorry(const double y)
{
  return pol5(y, 1.36807e+01, 1-8.65001e-02 , 3.98792e-06, 8.68754e-08, 
                 -2.11371e-12, -2.17610e-14);
}

double fidocorrz(const double z)
{
  return z > 1700? z + 31.5: pol5(z, -6.38165e+01, 1+4.44640e-02, -3.95760e-05, 
           1.05557e-07, 1.21091e-11, -2.26096e-14);
}

__attribute__((unused)) static double
phi(const double fex, const double fey, const double mx,
    const double my)
{
  return atan2(fey-my, fex-mx);
}

__attribute__((unused)) static double
theta(const double fex,const double fey,const double fez,
      const double mx, const double my, const double mz)
{
  return atan2(sqrt(pow(fex-mx,2)+pow(fey-my,2)),fez-mz);
}

double bamacorrz(const double z, const double e)
{
  // Romain's thesis's correction (eq. 7.24):
  return z + 7.466 
       + (0.008475 + 0.01029*e)*z
       - 1.053e-5*z*z
       + 2.694e-8*z*z*z;
}

double bamacorrxy(const double xy, const double e)
{
  return (1.013 - 7.0e-3*e)*xy + 0.0795e-3*xy*fabs(xy);
}

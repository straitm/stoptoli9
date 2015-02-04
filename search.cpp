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

static double pol4(const double gc)
{
  return  94.9593
       + 0.163906    * gc
       - 0.000228265 * gc*gc
       + 7.44587e-08 * gc*gc*gc
       - 1.06138e-11 * gc*gc*gc*gc;
}

static double pol3(const double x, const double p0, const double p1,
                   const double p2, const double p3)
{
  return p0
       + p1 * x
       + p2 * x*x
       + p3 * x*x*x;
}

double fidocorrz(const double z, const double gclen, const double th)
{
  const double shift = gclen <150? -200: pol4(gclen);  
  return z + shift*cos(th)
           - pol3(z, 62.1529, 0.00321, 9.02779e-06, -4.16564e-08);
}

double fidocorry(const double y, const double gclen, const double th,
                 const double phi)
{
  const double shift = gclen <150? -200:pol4(gclen);  
  return y - shift*sin(th)*sin(phi)
           - pol3(y, -9.75229-7.47685, -0.068273+0.0504024, 0, 0);
}

double fidocorrx(const double x, const double gclen, const double th,
                 const double phi)
{
  const double shift = gclen <150? -200:pol4(gclen);  
  return x - shift*sin(th)*cos(phi)
           - pol3(x, -9.75229+9.6452, -0.068273+0.0478563, 0, 0);
}

static double phi(const double fex, const double fey,
                  const double mx,  const double my)
{
  return atan2(fey-my, fex-mx);
}

static double theta(const double fex,const double fey,const double fez,
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
  // Using Romain's z correction also for x and y...
  const double efactor = 1; // 0.013/0.005;
  const double ifactor = 1; // 2;

  return xy + 7.466*ifactor 
       + (0.008475*ifactor + 0.01029*efactor*e)*xy
       - 1.053e-5*ifactor*xy*xy
       + 2.694e-8*ifactor*xy*xy*xy;
}

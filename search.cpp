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

static double fidopol3(const double gc)
{
  return -196.8
       + 0.119601    * gc
       - 2.00327e-05 * gc*gc
       + 1.02096e-08 * gc*gc*gc;
}

static double pol6(const double x, const double p0, const double p1,
                   const double p2, const double p3, const double p4,
                   const double p5, const double p6)
{
  return p0
       + p1 * x
       + p2 * x*x
       + p3 * x*x*x
       + p4 * x*x*x*x
       + p5 * x*x*x*x*x
       + p6 * x*x*x*x*x*x;
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
  if(gclen < 150) return z;
  const double shift = fidopol3(gclen);  
  return z - shift*cos(th)
           + pol6(z, -63.0931, -0.111571, 1.00597e-05,
                  1.04145e-07, -3.09931e-11, -2.0244e-14, 9.87013e-18);
}

double fidocorry(const double y, const double gclen, const double th,
                 const double phi)
{
  if(gclen < 150) return y;
  const double shift = fidopol3(gclen);  
  return y - shift*sin(th)*sin(phi)
           + pol3(y, 9.83751, -0.112057, -4.86037e-06, 3.6581e-08);
}

double fidocorrx(const double x, const double gclen, const double th,
                 const double phi)
{
  if(gclen < 150) return x;
  const double shift = fidopol3(gclen);  
  return x - shift*sin(th)*cos(phi)
           + pol3(x, -6.3569, -0.0897073, -1.72453e-06, 3.09912e-08);
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
  return (1.013 - 7.0e-3*e)*xy + 0.0795e-3*xy*fabs(xy);
}

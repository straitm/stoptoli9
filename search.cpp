#include <math.h>

bool lightnoise(const float qrms, const float mqtq, const float rmsts,
                const float qdiff)
{
  if(mqtq > 0.12 || mqtq < 0) return true;
  if(qdiff > 30e3) return true;
  if(rmsts >= 36 && qrms >= 464 - 8*rmsts) return true;
  return false;
}

double fidocorrz(const double z)
{
  return -86.2      
          +1.0396   *     z
          +4.14e-05 * pow(z, 2)
          -3.32e-08 * pow(z, 3)
          -1.02e-10 * pow(z, 4)
          +7.49e-14 * pow(z, 5)
          +6.17e-17 * pow(z, 6)
          -3.56e-20 * pow(z, 7)
          -1.01e-23 * pow(z, 8)
          +3.84e-27 * pow(z, 9);
}

double fidocorrxy(const double xy)
{
  return  2.76
          +1.0116 *       xy
          -1.64e-06 * pow(xy, 2)
          +1.91e-08 * pow(xy, 3);
}

double bamacorrxy(const double xy, const double e)
{
  return xy * (0.970 + 0.013*e);
}

double bamacorrz(const double z, const double e)
{
  return z  * (0.985 + 0.005*e);
}

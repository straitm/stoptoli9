double deadtime(const float fidoqid)
{
  return fidoqid/8300 * 0.063; // us
}

double targeff(const float fidoqid)
{
  const double ghtime = 5.5;
  if(deadtime(fidoqid) < ghtime) return 1;
  // There's a theralization turnon which more or less lines
  // up with the 5.5us trigger bounce time.  After that, approximate 
  // as exponential
  const double tau = 28.4;
  const double early = 0.0726;
  return (1-early)/exp(-ghtime/tau)*exp(-deadtime(fidoqid)/tau) + early;
}

double gceff(const float fidoqid)
{
  const double ghtime = 5.5;
  if(deadtime(fidoqid) < ghtime) return 1;
  // There is essentially no turn-on time for H capture.
  return exp(-deadtime(fidoqid)/179.) + (1-exp(-ghtime/179));
}

double eff(const float fidoqid, const double x, const double y,
           const double z)
{
  return fabs(z) < 1229 + 0.03*1150 && x*x+y*y < 1160*1160? 
         targeff(fidoqid): gceff(fidoqid);
}

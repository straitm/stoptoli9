
void fcn(int & npar, double * gin, double & chi2, double *par, int flag)
{
  const double a = 1 / 1000. / 2.028e-6 / 1000.;

  const double n1 = par[0], n2 = par[1], n3 = par[2];
  chi2 = pow((n1 + n3 - 1.42*a)/0.33*a , 2)
       + pow((n2 + n3 - 1.43*a)/0.29*a, 2);
}

  TMinuit mn(3);
void foo()
{
  mn.SetFCN(fcn);
  int err;
  mn.mnparm(0, "n1", 0.5, 0.1, 0, 10, err);
  mn.mnparm(1, "n2", 0.5, 0.1, 0, 10, err);
  mn.mnparm(2, "n3", 0.5, 0.1, 0, 10, err);
  mn.Command("set strategy 2");
  mn.Command("MIGRAD");
  mn.Command("MINOS");
}

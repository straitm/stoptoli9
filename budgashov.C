// Reanalysis of Budgashov et al., JETP 31 (1970) 651
// to not assume the 953keV level is zero.

TMinuit mn(6);

void fcn(int & npar, double * gin, double & chi2, double *par, int flag)
{
  const double a = 1 / 1000. / 2.028e-6 / 1000.;

  const double n1 = par[0], n2 = par[1], n3 = par[2];
  const double br_3t0 = par[3], br_2t0 = par[4], br_3t2 = par[5];

  const double eff = 0.03;

  chi2 = pow((n1
            + n2*(1-br_2t0)*(1-eff)
            + n3*(1-br_3t0)*(1-eff)
            - 1.42*a)/0.33*a, 2)

   + pow((n2*br_2t0
        + n3*(1-br_2t0-br_3t2*(1-br_2t0))*(1-eff)
        - 1.43*a)/0.29*a, 2)

   + pow((br_3t0 - 0.06 )/0.01 , 2)
   + pow((br_2t0 - 0.968)/0.004, 2)
   + pow((br_3t2 - 0.14 )/0.03 , 2)
  ;
}

void budgashov()
{
  mn.SetFCN(fcn);
  int err;
  mn.mnparm(0, "n1", 0.5, 0.1, 0, 10, err);
  mn.mnparm(1, "n2", 0.5, 0.1, 0, 10, err);
  mn.mnparm(2, "n3", 0.5, 0.1, 0, 10, err);
  mn.mnparm(3, "br_3t0", 0.06 , 0.01 , 0, 1, err);
  mn.mnparm(4, "br_2t0", 0.968, 0.004, 0, 1, err);
  mn.mnparm(5, "br_3t2", 0.14, 0.03, 0, 1, err);
  mn.Command("set strategy 2");
  mn.Command("MIGRAD");
  mn.Command("MINOS");
  mn.Command("show min");
}

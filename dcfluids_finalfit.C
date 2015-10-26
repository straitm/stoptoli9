/* Terrible name, dcfluids_finalfit is, but what it does is finds
 * the number of nitrogen and oxygen captures like the spreadsheet
 * dcfluids.ods used to.  It depends a lot on the fluid compositions,
 * so maybe it isn't that bad of a name.  Well, it really depends more
 * on the acrylic composition and wild guesses, but whatever.
 */

double E26();
double D23();
double C16();
double C17();
double C18();
double D20();
double D21();
double D22();
double D24();
double D25();
double D26();
double B20();
double B21();
double B22();
double B24();
double B25();
double I18();
double D8();
double C8();
double B11();
double C11();
double E3();
double I3();
double I16();
double I17();
double B23();
double I13();
double I11();
double C10();
double C4();
double D6();
double D5();
double D7();
double F4();
double H4();
double C3();
double D3();
double H3();
double B8();
double A8();
double B10();
double F2();
double H2();
double H5();
double I10();
double I12();
double I14();
double I19();
double J14();
double J19();
double B26();
double D27();
double C26();
double E27();
double G49();
double H49();
double I49();







double E26(){ return D26()*2.*0.85; }

double D23(){ return B23(); }

double C16(){ return 3.14159 * 1708*1708; }

double C17(){ return 3.14159*2.*1708.*1786.*2.; }

double C18(){ return (0.9*C17()+C16())/(2.*C16()+C17()); }

double D20(){ return C18()*814.*(15.999*2./(15.999*2.+12.011*5.+8.*1.008)); }

double D21(){ return D20()/((C10()+C11())*12./14.1); }

double D22(){ return D21()*I11(); }

double D24(){ return D22()*D23(); }

double D25(){ return D24()*0.4; }

double D26(){ return D25()*0.75; }

double B20(){ return 432.*(15.999*2./(15.999*2+12.011*5.+8.*1.008)); } 

double B21(){ return B20()/((C10()+C11())*12./14.1); }

double B22(){ return B21()*I11(); }

double B24(){ return B22()*B23(); }

double B25(){ return B24(); }

double I18(){ return B23()-I13(); }

double D8(){  return 1786-1229; }

double C8(){  return 1708-1150; }

double B11() { return ((A8()+C8())*(A8()+C8())*(B8()+D8())*2*3.14159
             + 2./3. * (A8()+C8())*0.03*(A8()+C8())*(A8()+C8())*3.14159)
             /1000000.-B10(); }

double C11() { return D5()*B11(); }

double E3()  { return 120.-D3(); }

double I3()  { return E3()*C3()*1000.; }

double I16() { return I3()/(C11() * 12./14.1)/1000. * 12./16.; }

double I17() { return I16()*I11(); }

double B23() { return 358.8; }

double I13() { return B23()*129.979327/353.192; }

double I11() { return (102.5/(102.5 + 1./0.002197))/
                        (37.9/(37.9 + 1./0.002197)); }

double C10() { return B10()*D5(); }

double C4() { return 15.999/(15.999 + 8*1.008 + 4.*12.011); }

double D6() { return 0.5; }

double D5() { return 0.804; }

double D7() { return D5()*D6()/100.*1000.; }

double F4() { return D7()*C4(); }

double H4() { return F4()*B10(); }

double C3() { return 15.999/(15.999 + 14.007 + 11.*1.008+15*12.011); }

double D3() { return 75.25; }

double H3() { return D3()*C3()*1000.; }

double B8() { return 1229.; }

double A8() { return 1150.; }

double B10(){ return (A8()*A8()*B8()*2*3.14159
         + 2./3. * A8()*0.03*A8()*A8()*3.14159)/1000000; }

double F2(){ return 0.99 * 6.*15.999/157.25;}

double H2(){ return F2()*B10(); }

double H5(){ return H2() + H3() + H4(); }

double I10(){ return H5()/(C10() * 12./14.1) * 12./16. / 1000.; }

double I12(){ return I10()*I11(); }

double I14(){ return I12() * I13(); }

double I19(){ return I17() * I18(); }

double J14(){ return I14() * 2.; }

double J19(){ return I19()*2.; }

double B26(){ return B25()*0.85; }

double D27(){ return D26()*0.43;}

double C26(){ return B25()*2*0.95;}

double E27(){ return E26()*0.47;}

double G49(){ return I14() + I19() + B26() + D27(); }

double H49(){ return J14() + J19() + C26() + E27(); }

double I49(){ return (G49() + H49())/2.; }

void dcfluids_finalfit()
{
  const double n_o16cap_beta = I49();

  printf("%f\n", C10());

  printf("TECHNOTE 5.3: Gaussian central value of number of O-16 "
    "captures per day: %.8f\n", n_o16cap_beta);
}

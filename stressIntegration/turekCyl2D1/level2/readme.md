Copied from
  Turek only-cylinder validation
  Case 2D1 using half-way bounceback
  level 1: 82
  Result: Cd=5.823, Cl=0.01206

Target: Acheive same using stress integration instead of momemtum-exchange

Successful on 07/12/22 commit(4d99e94d55328fc7d6648126871b7629c01f10f5)
It uses curved BC in contrast to half-way bounceback

H=82, noOfPoints=400, Using all relevant direction
  time      rhoAvg        Cd(Mom-Ex)    Cl(Mom-Ex)    Cd(StressInt) Cl(StressInt)
  100000    0.1004E+01    0.5737E+01    0.1066E-01    0.5844E+01    0.9344E-02

H=82, noOfPoints=400, Using only 1 direction for inerpolation along max dot with normal
  time      rhoAvg        Cd(Mom-Ex)    Cl(Mom-Ex)    Cd(StressInt) Cl(StressInt)
  100000    0.1004E+01    0.5737E+01    0.1066E-01    0.5699E+01    0.8698E-02

H=205, noOfPoints=400, Using only 1 direction
  time              LBM time  rhoAvg        Cd(Mom-Ex)    Cl(Mom-Ex)    Cd(StressInt) Cl(StressInt)
  0.125000E+02      100000    0.1000E+01    0.5627E+01    0.1060E-01    0.5625E+01    0.9928E-02  
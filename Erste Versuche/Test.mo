within ;
model Test
  import SI = Modelica.SIunits;

  PowerSource PS1;
  Elektrolyser E;
  BurnCell B;
  PowerSink PS2;

equation
  connect(PS1.e,E.e);
  connect(B.h,E.h);
  connect(B.e,PS2.e);
end Test;

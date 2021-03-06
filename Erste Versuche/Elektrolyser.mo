within ;
model Elektrolyser
  extends Elektrolyse;
  import SI = Modelica.SIunits;



  parameter Real n = 0.95;

  SI.Current i;
  SI.Voltage u;
  SI.Pressure p;
  SI.MassFlowRate mF;

equation
  i + e.i = 0;
  u = e.u;
  p = h.p;
  mF + h.mF = 0;
  p=10000;

  mF + u*i*n = 0;

end Elektrolyser;

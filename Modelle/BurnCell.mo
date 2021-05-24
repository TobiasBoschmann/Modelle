within ;
model BurnCell
  import SI = Modelica.SIunits;

  HydrogenConnector h;
  ElectricConnector e;

  parameter Real n = 0.9 "Effizienz";
  parameter SI.Voltage u = 230;

  SI.Current i;
  SI.Pressure p;
  SI.MassFlowRate mF;

equation
  i + e.i = 0;
  e.u = u;
  p = h.p;
  mF + h.mF = 0;

  mF*n + u*i = 0;

end BurnCell;

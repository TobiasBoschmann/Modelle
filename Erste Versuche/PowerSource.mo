within ;
model PowerSource
  import SI = Modelica.SIunits;

  ElectricConnector e;

  parameter SI.Current i = 5;
  parameter SI.Voltage u = 230;

  SI.Energy InputPower( start= 0);

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)));

equation
  e.i + i =0;
  e.u = u;

  der(InputPower) = u*i;

end PowerSource;

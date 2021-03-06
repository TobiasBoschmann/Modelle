within ;
model PowerSink
  import SI = Modelica.SIunits;

  ElectricConnector e;

  SI.Current i;
  SI.Voltage u;
  SI.Energy OutputPower( start= 0);


  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)));

equation
  e.i + i =0;
  e.u = u;
  der(OutputPower) = -u * i;

end PowerSink;

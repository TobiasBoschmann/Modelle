within ;
model HydrogenSink
  import SI = Modelica.SIunits;
  HydrogenConnector h;
  SI.Pressure p;
  SI.MassFlowRate mF;
  SI.Mass m( start=0);
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)));

equation
  mF + h.mF = 0;
  p = h.p;
  der(m) = -mF;


end HydrogenSink;

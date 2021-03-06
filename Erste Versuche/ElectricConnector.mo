within ;
connector ElectricConnector
  import SI = Modelica.SIunits;
  flow SI.Current i;
  SI.Voltage u;
               annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end ElectricConnector;

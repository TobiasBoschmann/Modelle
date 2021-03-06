within ;
partial model Elektrolyse
  import SI = Modelica.SIunits;

  //Konnektoren
  ElectricConnector e;
  HydrogenConnector h;


  //Parameter Thermodynamik
  SI.MolarInternalEnergy G = 236483 "freie Gische Energie bei T=298,15K und P=1atm";
  SI.MolarInternalEnergy H = 285830 "Reaktionsenthalpie bei T=298,15K und P=1atm";
  Real n = 2 "Anzahl der ausgetauschten Elektronen pro Elementarreaktion";
  SI.Temprature T;


  //Parameter Elektronik
  SI.Voltage V_rev_0 = G/(n*Modelica.Constants.F) "reversible Cell-Voltage";
  SI.Voltage V_rev;
  SI.Voltage V_th =  H/(n*Modelica.Constants.F) "thermodynamic Cell-Voltage";
  SI.Voltage V_cell "produced/needed Cell-Voltage";
  SI.Voltage V_act "Voltage-losse due to activation processes";
  SI.Voltage V_ohm "ohmic Voltage-losses";
  SI.Current I;

equation
  V_rev = T*Modelica.Constants.R/(n*Modelica.Constants.F)*ln(100);
  V_cell = V_rev + V_act + V_ohm;


end Elektrolyse;

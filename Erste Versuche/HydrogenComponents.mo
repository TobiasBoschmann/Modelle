within ;
package HydrogenComponents
  model Elektrolyser
    extends Elektrolyse;
    import SI = Modelica.SIunits;
    import C = Modelica.Constants;

    parameter Real alpha_m(unit="1") = 0.5;
    SI.CurrentDensity i_0 = 0.05;

  equation

    N_H2 = +I/(z*C.F);

    //p = h.p;
    //p=10000;



    V_cell = V_rev + V_act + V_ohm;


    //Lückenfüller:
    T = 298.15;

    V_ohm = +0.0005*I;
    V_act = +2*C.R*T/(alpha_m*C.F)*asinh(i/(2*i_0));

    p_H2 = p_O2;
    p_H2 + p_O2 = 100000;

    if P_el==0 then
      eff = 0;
    else
      eff = N_H2/P_el;
    end if;

  end Elektrolyser;

  partial model Elektrolyse
    import SI = Modelica.SIunits;
    import C = Modelica.Constants;

    //Konnektoren
    Modelica.Electrical.Analog.Interfaces.PositivePin p;
    Modelica.Electrical.Analog.Interfaces.NegativePin n;
    HydrogenConnector h;

    parameter SI.Area A = 1 "Area of exchange-membrane";

    //Parameter Thermodynamik

    //SI.Pressure p;
    SI.MolarFlowRate N_H2;

    SI.MolarInternalEnergy G = 236483 "freie Gische Energie bei T=298,15K und P=1atm";
    SI.MolarInternalEnergy H = 285830 "Reaktionsenthalpie bei T=298,15K und P=1atm";
    Real z = 2 "Anzahl der ausgetauschten Elektronen pro Elementarreaktion";
    SI.Temperature T_0 = 298.15 "standard temperature";
    SI.Temperature T "temperature of incomming Water";
    SI.Pressure p_H2 "partial pressure H2";
    SI.Pressure p_O2 "partial pressure O2";
    Real a_H2O(unit="PA^1.5") = 1 "water activity (equals 1 for liquid water)";

    //Parameter Elektronik
    SI.Voltage V_cell "produced/needed Cell-Voltage";

    SI.Voltage V_rev_0 = G/(z*C.F) "reversible Cell-Voltage at T=298,15K and P=1atm";
    SI.Voltage V_rev "reversible Cell-Voltage";


    SI.Voltage V_act "Voltage-losse due to activation processes";
    SI.Voltage V_ohm "ohmic Voltage-losses";
    SI.Current I;
    SI.CurrentDensity i;

    SI.Power P_el "electrical power in/output";

    SI.Efficiency eff;

  equation
    //connector equations

    h.N = -N_H2;
    p.i = I "Current from positiv to negative Pin";
    V_cell = p.v - n.v;

    //thermodynamic equations
    V_rev = V_rev_0 + T*C.R/(z*C.F)*log(p_H2*sqrt(p_O2)/a_H2O) "Nernst Equation";

    I=A*i;
    P_el = V_cell*I;

  end Elektrolyse;

  connector HydrogenConnector
    import SI = Modelica.SIunits;
    flow SI.MolarFlowRate N "molar flowrate of Hydrogen";
    //SI.Pressure p;

  end HydrogenConnector;

  model PowerSource
    import SI = Modelica.SIunits;

    Modelica.Electrical.Analog.Interfaces.PositivePin p;
    Modelica.Electrical.Analog.Interfaces.NegativePin n;

    SI.Current i( start=0);
    SI.Voltage v;

    SI.Energy InputPower(start = 0);

  equation
    der(i) = 0.02;

    p.i = -i;
    p.i + n.i = 0;

    n.v = 0;
    p.v - n.v = v;

   der(InputPower) = v*i;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end PowerSource;

  model BurnCell
    extends Elektrolyse;

  equation

    N_H2 = -I/(z*Modelica.Constants.F);

    //p = h.p;
    //p=10000;



    V_cell = V_rev + V_act + V_ohm;


    //Lückenfüller:
    T = 298.15;

    V_ohm = -0.05*I;
    V_act = -0.5*sqrt(abs(I));

    p_H2 = p_O2;
    p_H2 + p_O2 = 100000;

    if N_H2==0 then
      eff = 0;
    else
      eff = P_el/N_H2;
    end if;

  end BurnCell;

  model PowerSink
    import SI = Modelica.SIunits;

    Modelica.Electrical.Analog.Interfaces.PositivePin p;
    Modelica.Electrical.Analog.Interfaces.NegativePin n;

    SI.Current i( start = 0.001);
    SI.Voltage v;
    SI.Energy InputPower( start= 0);

  equation
    p.i = i;
    p.i + n.i = 0;

    n.v = 0;
    p.v - n.v = v;

    der(InputPower) = v * i;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end PowerSink;

  model HydrogenSink
    import SI = Modelica.SIunits;
    HydrogenConnector h;
    //SI.Pressure p;
    SI.MolarFlowRate N_H2;
    SI.AmountOfSubstance n( start=0);

  equation
    N_H2 = h.N;
    //p = h.p;
    der(n) = N_H2;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end HydrogenSink;

  model HydrogenSource
    import SI = Modelica.SIunits;
    HydrogenConnector h;
    //SI.Pressure p;

    SI.MolarFlowRate N_H2( start = 0);
    SI.AmountOfSubstance n( start=0);

  equation
    N_H2 + h.N = 0;
    der(N_H2) = 0.0001;
    //p = h.p;
    der(n) = N_H2;

  end HydrogenSource;

  package Examples
    model ExElOnly
      extends Modelica.Icons.Example;
      import SI = Modelica.SIunits;

      PowerSource PI;
      Elektrolyser E;
      HydrogenSink HS;

    equation
      connect(PI.p,E.p);
      connect(PI.n,E.n);

      connect(HS.h,E.h);

    end ExElOnly;

    model ExBCOnly
      extends Modelica.Icons.Example;
      import SI = Modelica.SIunits;

      HydrogenSource HS;
      BurnCell B;
      PowerSink PS;

    equation
      connect(PS.p,B.p);
      connect(PS.n,B.n);

      connect(HS.h,B.h);

    end ExBCOnly;
  end Examples;
  annotation (uses(Modelica(version="3.2.2")));
end HydrogenComponents;

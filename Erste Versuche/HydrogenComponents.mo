within ;
package HydrogenComponents
  model Elektrolyser
    extends Elektrolyse;
    import SI = Modelica.SIunits;
    import C = Modelica.Constants;

    parameter Real alpha_m(unit="1") = 0.5;
    SI.CurrentDensity i_0 = 21.8*10^(-4) "hat mit Material der Membran zu tun";

    Real alpha = 1 "Tafel equation parameter";

    //PEM-Parameter:
    SI.CurrentDensity i_max = 0.5 "maximal stromdichte der Membran";
    SI.Length delta = 0.35*10^(-3) "membran thickness";
    SI.Conductivity sigma "conductivity of membran";
    Real lambda = 22 "water content of membran (dry: ~0.5; water saturated gas: ~13; liquid water: ~22)";

  equation

    N_H2 = +I/(z*C.F);

    //p = h.p;
    //p=10000;



    V_cell = V_rev + V_act + V_ohm + V_diff;

    V_rev = V_rev_0 + T*C.R/(z*C.F)*log(p_H2/p_0*sqrt(p_O2/p_0)/a_H2O) "Nernst Equation";

    if i == 0 then
      V_act = 0;
    else
      V_act = C.R*T/(alpha*z*C.F)*log(i/i_0) "Tafel equation";
    end if;

    V_ohm = I*delta/sigma;
    sigma = (0.005139*lambda -0.00326)*exp(1268*(1/303-1/T));

    V_diff = 0 "- C.R*T/(z*C.F)*log(1-i/i_max)";//vernachlässigbar?

    //Lückenfüller:
    T = 293.15;

    p_H2 = p_0 + i*24 - p_H2O;
    p_O2 = p_0 + i*28 - p_H2O;

    if P_el==0 then
      eff = 0;
    else
      eff = N_H2*H/P_el;
    end if;


  end Elektrolyser;

  partial model Elektrolyse
    import SI = Modelica.SIunits;
    import C = Modelica.Constants;

    //Konnektoren
    Modelica.Electrical.Analog.Interfaces.PositivePin p;
    Modelica.Electrical.Analog.Interfaces.NegativePin n;
    HydrogenConnector h;

    parameter SI.Area A = 25*10^(-6) "Area of exchange-membrane";

    //Parameter Thermodynamik

    SI.Pressure p_0 = 100000;
    SI.MolarFlowRate N_H2;

    SI.MolarInternalEnergy G = 236483 "freie Gische Energie bei T=298,15K und P=1atm";
    SI.MolarInternalEnergy H = 285830 "Reaktionsenthalpie bei T=298,15K und P=1atm";
    Real z = 2 "Anzahl der ausgetauschten Elektronen pro Elementarreaktion";
    SI.Temperature T_0 = 298.15 "standard temperature";
    SI.Temperature T "temperature of incomming Water";
    SI.Pressure p_H2 "partial pressure H2";
    SI.Pressure p_O2 "partial pressure O2";
    parameter SI.Pressure p_H2O = 46000 "partial pressure H2O";

    Real a_H2O(unit="PA^1.5") = 1 "water activity (equals 1 for liquid water)";

    //Parameter Elektronik
    SI.Voltage V_cell "produced/needed Cell-Voltage";

    SI.Voltage V_rev_0 = G/(z*C.F) "reversible Cell-Voltage at T=298,15K and P=1atm";
    SI.Voltage V_rev "reversible Cell-Voltage";


    SI.Voltage V_act "activation process resistance";
    SI.Voltage V_ohm "ohmic resistance";
    SI.Voltage V_diff "diffusion resistance";
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
    der(i) = 0.2;

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

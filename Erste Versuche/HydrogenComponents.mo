within ;
package HydrogenComponents
  model Elektrolyser
    import SI = Modelica.SIunits;
    import C = Modelica.Constants;

   //Konnektoren
    AixLib.Utilities.MassTransfer.MassPort O2_Output;
    AixLib.Utilities.MassTransfer.MassPort H2_Output;
    AixLib.Utilities.MassTransfer.MassPort H2O_Input;
    //Modelica.Electrical.Analog.Interfaces.PositivePin p;
    //Modelica.Electrical.Analog.Interfaces.NegativePin n;
    Modelica.Blocks.Interfaces.RealInput PowerInput(
      final quantity="Power",
      final unit="W")
    "Input power for the elektrolyser";
    //Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort; // T - Port Temprerature; Q_Flow Heatflow

    //Parameter
    parameter Real alpha_m(unit="1") = 0.5 "für Überspannung";
    parameter SI.Area A = 25*10^(-4) "active area";
    SI.CurrentDensity i_0 = 21.8*10^(-4) "hat mit Material der Membran zu tun";

    Real alpha = 1 "Tafel equation parameter";

    //PEM-Parameter:
    parameter SI.CurrentDensity i_max = 2*10^4 "maximal stromdichte";
    parameter SI.Length delta = 0.175*10^(-2) "membran thickness";
    SI.Conductivity sigma_pem "conductivity of membran";
    parameter Real lambda = 22 "water content of membran (dry: ~0.5; water saturated gas: ~13; liquid water: ~22)";
    Real R_mem;
    Real R_ele = 0.0096;

    //Alkaline parameter
    parameter Real m = 0.2 "molar Concentration of the electrolyt-solution";
    SI.Conductivity sigma_alk "conductivity of elektrolyte";

    //SO Parameter
    SI.Conductivity sigma_so "conductivity of elektrolyte";


    //Parameter Thermodynamik


    SI.MolarFlowRate n_H2;
    constant SI.MolarMass M_H2 = 2.02;
    parameter SI.Pressure p_H2 = 400000 "partial pressure H2";

    SI.MolarFlowRate n_O2;
    constant SI.MolarMass M_O2 = 16;
    parameter SI.Pressure p_O2 = 400000 "partial pressure O2";

    SI.MolarFlowRate n_H2O;
    constant SI.MolarMass M_H2O = 18.02;
    SI.Pressure p_H2O "partial pressure H2O";
    parameter SI.Temperature T_H2O = 298.15 "temperature of incomming Water";
    constant SI.SpecificHeatCapacity cp_H2O = 1004;

    SI.MolarInternalEnergy G_r_0 = 236483 "free gibbs reaktion energy at T=298,15K and P=1atm";
    SI.MolarInternalEnergy H_r "reaktion enthalpy";
    SI.MolarInternalEnergy H_r_0 = 285830 "reaktion enthalpy at T=298,15K and P=1atm";

    Real z = 2 "nuber of exchanged elektrons per elementary reaktion";
    constant SI.Temperature T_0 = 298.15 "standard temperature";
    constant SI.Pressure p_0 = 103105 "standard pressure";
    SI.Temperature T = 298.15 "process temperature";

    Real a_H2O(unit="1") = 1 "water activity (equals 1 for liquid water)";

    //Parameter Elektronik
    SI.Voltage V_cell "needed Cell-Voltage";

    SI.Voltage V_th_0 = H_r_0/(z*C.F) "thermo-neutral Cell-Voltage at T=298,15K and P=1atm";
    SI.Voltage V_rev_0 = G_r_0/(z*C.F) "reversible Cell-Voltage at T=298,15K and P=1atm";
    SI.Voltage V_rev "reversible Cell-Voltage";


    SI.Voltage V_act "activation overpotetial";
    SI.Voltage V_ohm "ohmic resistance";
    SI.Voltage V_diff "diffusion overpotetial";
    SI.Current I;
    SI.CurrentDensity i;

    SI.Power P_el "electrical input power";
    SI.HeatFlowRate Q "cooling power";

    SI.Efficiency eff;

  equation
    //connector equations

    H2_Output.m_flow = M_H2 * n_H2;
    H2_Output.p = p_H2;

    O2_Output.m_flow = M_O2 * n_O2;
    O2_Output.p = p_O2;

    H2O_Input.m_flow = M_H2O * n_H2O;
    H2O_Input.p = p_H2O;

    //p.i = I "Current from positiv to negative Pin";
    //V_cell = p.v - n.v;
    PowerInput = V_cell*I;

    //elektro-chemical equations

    I=A*i;
    P_el = V_cell*I;

    n_H2  = I/(z*C.F);
    n_H2O = n_H2;
    n_O2  = 0.5 * n_H2;


    V_cell = V_rev + V_act + V_ohm + V_diff;

    V_rev = V_rev_0 - 8.5*10^(-4)*(T-298) + T*C.R/(z*C.F)*log(p_H2/p_0*sqrt(p_O2/p_0)/a_H2O) "Temperature dependence and Nernst Equation";

    if i == 0 then
      V_act = 0;
    else
      V_act = C.R*T/(alpha_m*2*C.F)*log(i/i_0) "from Tafel equation";
    end if;

    V_ohm = I*(R_mem + R_ele);
    R_mem=delta/sigma_pem;

    sigma_pem = (0.005139*lambda -0.00326)*exp(1268*(1/303-1/T));
    sigma_alk = -2.04*m -0.0028*m^2  + 0.005332*m*T +207.2*m/T+0.001043*m^3-0.0000003*m^2*T^2;
    sigma_so = 3.34*10^4*exp(-10300/T);

    if i<i_max then
      V_diff = C.R*T/(z*C.F)*log(i_max/(i_max-i));
    else
      V_diff = 0;
    end if;


    //thermodynamic equations
    H_r = H_r_0;

    0 = n_H2O*cp_H2O*M_H2O*(T_H2O - T) - I/(z*C.F)*H_r + P_el - Q;

    //Lückenfüller:
    p_O2 = p_H2O;

    if P_el==0 then
      eff = 0;
    else
      eff = n_H2*H_r_0/P_el;
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

    //Modelica.Electrical.Analog.Interfaces.PositivePin p;
    //Modelica.Electrical.Analog.Interfaces.NegativePin n;
     Modelica.Blocks.Interfaces.RealInput PowerOutput(
      final quantity="Power",
      final unit="W");

    //SI.Current i( start=0);
    //SI.Voltage v;

  equation
    der(PowerOutput) = 100;
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
    AixLib.Utilities.MassTransfer.MassPort O2_Input;
    AixLib.Utilities.MassTransfer.MassPort H2_Input;
    AixLib.Utilities.MassTransfer.MassPort H2O_Output;


    SI.MolarFlowRate n_H2;
    constant SI.MolarMass M_H2 = 2.02;
    SI.Pressure p_H2 "partial pressure H2";

    SI.MolarFlowRate n_O2;
    constant SI.MolarMass M_O2 = 16;
    SI.Pressure p_O2 "partial pressure O2";

    SI.MolarFlowRate n_H2O;
    constant SI.MolarMass M_H2O = 18.02;
    SI.Pressure p_H2O "partial pressure H2O";

    SI.AmountOfSubstance n( start=0);

  equation
    H2_Input.m_flow = M_H2 * n_H2;
    H2_Input.p = p_H2;

    O2_Input.m_flow = M_O2 * n_O2;
    O2_Input.p = p_O2;

    H2O_Output.m_flow = M_H2O * n_H2O;
    H2O_Output.p = p_H2O;




    der(n) = n_H2;

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
      connect(PI.PowerOutput,E.PowerInput);

      connect(HS.H2_Input,E.H2_Output);
      connect(HS.O2_Input,E.O2_Output);
      connect(HS.H2O_Output,E.H2O_Input);
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

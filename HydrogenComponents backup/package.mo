within ;
package HydrogenComponents
  model Elektrolyseur
    extends Elektrolyse;

    import SI = Modelica.SIunits;
    import C = Modelica.Constants;


  equation
    //elektro-chemical equations
    U_zell = U_rev + U_akt + U_ohm + U_konz;

    //thermodynamic equations
    0 = n_H2*(0.5*cp_O2+cp_H2)*(T_u - T) + I*(U_zell-U_th) - Q;


    if P_el==0 then
      eff = 0;
    else
      eff = n_H2*H_r_0/P_el;
    end if;


  end Elektrolyseur;

  model BurnCell
    extends Elektrolyse;

    import SI = Modelica.SIunits;
    import C = Modelica.Constants;


  equation
    //elektro-chemical equations
    U_zell = U_rev - U_akt - U_ohm - U_konz;

    //thermodynamic equations
    0 = -n_H2*(0.5*cp_O2+cp_H2)*(T_u - T) - I*(U_zell-U_th_0) - Q;

    if n_H2==0 then
      eff = 0;
    else
      eff = P_el/(n_H2*H_r_0);
    end if;

  end BurnCell;

  model Power
    import SI = Modelica.SIunits;

     Modelica.Blocks.Interfaces.RealInput Power(start=0,
      final quantity="Power",
      final unit="W");

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end Power;

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

  model Speicher
    AixLib.Utilities.MassTransfer.MassPort Anschluss;

    import SI = Modelica.SIunits;
    import C = Modelica.Constants;

    parameter SI.Volume V=1;
    parameter SI.Temperature T = 298.15;
    parameter SI.MolarMass M = 2.02 "2.02 für H2, 16 für O2";
    SI.pressure p;
    parameter Real eff = 0.5;

    Modelica.Blocks.Interfaces.RealInput Power(
      final quantity="Power",
      final unit="W");

  equation

    if p > Anschluss.p then
      Power = eff*(p-Anschluss.p)*Anschluss.m_flow*C.R*T/ (Anschluss.p*M);
    else
      Power = 0;
    end if;


    der(p) = Anschluss.m_flow*C.R*T / (V*M);

  end Speicher;

  model HydrogenSink
    import SI = Modelica.SIunits;
    AixLib.Utilities.MassTransfer.MassPort O2;
    AixLib.Utilities.MassTransfer.MassPort H2;
    AixLib.Utilities.MassTransfer.MassPort H2O;

  protected
            constant Real H2Verbrauch[:,:] = [0, 1; 1,0; 2, 0.3; 3,4e-14; 4,4e-14; 5,2e-14; 6,1e-14];

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
    H2.m_flow = M_H2 * n_H2;
    H2.p = p_H2;

    O2.m_flow = M_O2 * n_O2;
    O2.p = p_O2;

    H2O.m_flow = M_H2O * n_H2O;
    H2O.p = p_H2O;

    n_H2= -DataFromTable(time,H2Verbrauch);


    der(n) = n_H2;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end HydrogenSink;

  model HydrogenSource
    import SI = Modelica.SIunits;

    AixLib.Utilities.MassTransfer.MassPort O2;
    AixLib.Utilities.MassTransfer.MassPort H2;
    AixLib.Utilities.MassTransfer.MassPort H2O;

  protected
            constant Real H2Verbrauch[:,:] = [0, 0.2; 1,0; 2, 0.3; 3,4e-14; 4,4e-14; 5,2e-14; 6,1e-14];//Vorläufig

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
    H2.m_flow = M_H2 * n_H2;
    H2.p = p_H2;

    O2.m_flow = M_O2 * n_O2;
    O2.p = p_O2;

    H2O.m_flow = M_H2O * n_H2O;
    H2O.p = p_H2O;

    n_H2= -DataFromTable(time,H2Verbrauch);

    der(n) = n_H2;
  end HydrogenSource;

  partial model Elektrolyse
    import SI = Modelica.SIunits;
    import C = Modelica.Constants;

   //Konnektoren
    AixLib.Utilities.MassTransfer.MassPort O2
    annotation(Placement(
    transformation(extent={{-120,-90},{-100,-70}}),
    iconTransformation(extent={{-114,-50},{-94,-30}})));
    AixLib.Utilities.MassTransfer.MassPort H2
    annotation(Placement(
    transformation(extent={{-120,-10},{-100,10}}),
    iconTransformation(extent={{-114,-10},{-94,10}})));
    AixLib.Utilities.MassTransfer.MassPort H2O
    annotation(Placement(
    transformation(extent={{-120,70},{-100,90}}),
    iconTransformation(extent={{-114,30},{-94,50}})));
    Modelica.Blocks.Interfaces.RealInput Power(
      final quantity="Power",
      final unit="W")
    "Input/Output power for the elektrolyser"
    annotation(Placement(
    transformation(extent={{100,30},{120,50}}),
    iconTransformation(extent={{100,30},{120,50}})));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort
    annotation(Placement(
    transformation(extent={{100,-50},{120,-30}}),
    iconTransformation(extent={{100,-50},{120,-30}})));  // T - Port Temprerature; Q_Flow Heatflow


    //Parameter
    parameter SI.Area A = 10 "active area";
    parameter Real n = 60 "Anzahl der in parallelgeschalteten Zellen";
     //PEM-Parameter:

    parameter Real alpha_m(unit="1") = 0.7353 "für Überspannung";
    SI.CurrentDensity i_0_pem = 1.08*10^(-4) "hat mit Material der Membran zu tun";
    parameter SI.CurrentDensity i_max = 2*10^4 "maximal stromdichte";
    parameter SI.Length delta_pem = 0.15*10^(-3) "membran thickness";
    SI.Conductivity sigma_pem "conductivity of membran";
    parameter Real lambda = 22 "water content of membran (dry: ~0.5; water saturated gas: ~13; liquid water: ~22)";
    Real t_min_pem = 0.1 "minimale Teillast";
    constant SI.HeatCapacity Cp_0 = 162116 "Literaturwert";
    constant SI.Area A_0 = 60*290*10^(-4) "Literaturwert";

    //Alkaline parameter
    parameter Real m = 0.2 "molar Concentration of the electrolyt-solution";
    SI.Conductivity sigma_alk "conductivity of elektrolyte";
    Real t_min_alk = 0.2 "minimale Teillast";
    constant Real t_start_alk = 2 "Verhältniss der durchschnittlichen Anlaufzeiten von Alkalischen zu PEM Elektrolyseuren";

    //SO Parameter
    SI.Conductivity sigma_so "conductivity of elektrolyte";
    Real t_min_so = 0.3 "minimale Teillast";



    //Parameter Thermodynamik
    SI.MolarFlowRate n_H2;
    constant SI.MolarMass M_H2 = 2.02;
    parameter SI.Pressure p_H2 = 400000 "partial pressure H2";
    constant SI.SpecificHeatCapacity cp_H2 = 28.9;

    SI.MolarFlowRate n_O2;
    constant SI.MolarMass M_O2 = 16;
    parameter SI.Pressure p_O2 = 400000 "partial pressure O2";
    constant SI.SpecificHeatCapacity cp_O2 = 29.8;

    SI.MolarFlowRate n_H2O;
    constant SI.MolarMass M_H2O = 18.02;
    SI.Pressure p_H2O "partial pressure H2O";
    Real a_H2O(unit="1") = 1 "water activity (equals 1 for liquid water)";



    SI.MolarInternalEnergy G_r_0 = 236483 "free gibbs reaktion energy at T=298,15K and P=1atm";
    SI.MolarInternalEnergy H_r_0 = 285830 "reaktion enthalpy at T=298,15K and P=1atm";

    Real z = 2 "nuber of exchanged elektrons per elementary reaktion";
    constant SI.Temperature T_0 = 298.15 "standard temperature";
    constant SI.Pressure p_0 = 103105 "standard pressure";
    SI.Temperature T = (273.15 + 80) "process temperature";
    parameter SI.Temperature T_u = 293.15 "Umgebungstemperatur";


    //Parameter Elektronik
    SI.Voltage U_zell "needed Cell-Voltage";

    SI.Voltage U_th_0 = H_r_0/(z*C.F) "thermo-neutral Cell-Voltage at T=298,15K and P=1atm";
    SI.Voltage U_rev_0 = G_r_0/(z*C.F) "reversible Cell-Voltage at T=298,15K and P=1atm";
    SI.Voltage U_rev "reversible Cell-Voltage";


    SI.Voltage U_akt "activation overpotetial";
    SI.Voltage U_ohm "ohmic resistance";
    SI.Voltage U_konz "diffusion overpotetial";
    SI.Current I;
    SI.CurrentDensity i;

    SI.Power P_el "electrical input power";
    SI.HeatFlowRate Q "cooling power";

    SI.Efficiency eff;

  equation

    H2.m_flow = M_H2 * n_H2;
    H2.p = p_H2;

    O2.m_flow = M_O2 * n_O2;
    O2.p = p_O2;

    H2O.m_flow = M_H2O * n_H2O;
    H2O.p = p_H2O;

    Power = U_zell*I;

    I=A*i;
    P_el = Power;

    //elektro-chemical equations
    n_H2  = I/(z*C.F);
    n_H2O = n_H2;
    n_O2  = 0.5 * n_H2;

    U_rev = U_rev_0 - 8.5*10^(-4)*(T-298) + T*C.R/(z*C.F)*log(p_H2/p_0*sqrt(p_O2/p_0)/a_H2O) "Temperature dependence and Nernst Equation";

    if i == 0 then
      U_akt = 0;
    else
      U_akt = C.R*T/(alpha_m*2*C.F)*log(i/i_0_pem) "from Tafel equation";
    end if;

    U_ohm = I/n*delta_pem/sigma_pem;

    sigma_pem = (0.005139*lambda -0.00326)*exp(1268*(1/303-1/T));
    sigma_alk = -2.04*m -0.0028*m^2  + 0.005332*m*T +207.2*m/T+0.001043*m^3-0.0000003*m^2*T^2;
    sigma_so = 3.34*10^4*exp(-10300/T);

    if i<i_max then
      U_konz = C.R*T/(z*C.F)*log(i_max/(i_max-i));
    else
      U_konz = 0;
    end if;


    //Lückenfüller:
    p_O2 = p_H2O;

    annotation (El, Icon(graphics={
          Rectangle(
            extent={{-78,92},{76,-92}},
            lineColor={0,0,0},
            fillColor={158,158,158},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-90,100},{-70,-100}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{70,100},{90,-100}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid)}));
  end Elektrolyse;

  package Examples
    model ExElOnly
      extends Modelica.Icons.Example;
      import SI = Modelica.SIunits;

      Power P;
      Elektrolyseur E;
      HydrogenSink HS;

    equation
      connect(P.Power,E.Power);

      connect(HS.H2,E.H2);
      connect(HS.O2,E.O2);
      connect(HS.H2O,E.H2O);
    end ExElOnly;

    model ExBCOnly
      extends Modelica.Icons.Example;
      import SI = Modelica.SIunits;

      HydrogenSource HS;
      BurnCell B;
      Power P;

    equation
      connect(P.Power,B.Power);

      connect(HS.H2,B.H2);
      connect(HS.O2,B.O2);
      connect(HS.H2O,B.H2O);

    end ExBCOnly;

    model PVElektrolyseur
      extends Modelica.Icons.Example;

      Modelica.Blocks.Interfaces.RealOutput Power(
        final quantity="Power",
        final unit="W")
        "Output Power of the PV system including the inverter"
        annotation (Placement(transformation(extent={{-2,58},{18,78}})));


      AixLib.BoundaryConditions.WeatherData.Old.WeatherTRY.Weather Weather(
        Latitude=49.5,
        Longitude=8.5,
        GroundReflection=0.2,
        tableName="wetter",
        extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic,
        Wind_dir=false,
        Air_temp=true,
        Wind_speed=false,
        SOD=AixLib.DataBase.Weather.SurfaceOrientation.SurfaceOrientationData_N_E_S_W_RoofN_Roof_S(),
        fileName=Modelica.Utilities.Files.loadResource(
            "modelica://AixLib/Resources/WeatherData/TRY2010_12_Jahr_Modelica-Library.txt"))
        "Weather data input for simulation of PV power "
        annotation (Placement(transformation(extent={{-93,49},{-68,66}})));

      AixLib.Electrical.PVSystem.PVSystem PVsystem(
        MaxOutputPower=4000,
        NumberOfPanels=5,
        data=AixLib.DataBase.SolarElectric.SymphonyEnergySE6M181())
        "PV system model including the inverter"
        annotation (Placement(transformation(extent={{-42,46},{-22,66}})));


    equation
      connect(Weather.AirTemp, PVsystem.TOutside) annotation (Line(points={{
              -67.1667,60.05},{-62.5834,60.05},{-62.5834,63.6},{-44,63.6}},
                                                                   color={0,0,127}));
      connect(PVsystem.PVPowerW, Power)
        annotation (Line(points={{-21,56},{-8,56},{-8,68},{8,68}},
                                                   color={0,0,127}));
      connect(Weather.SolarRadiation_OrientedSurfaces[6], PVsystem.IcTotalRad)
        annotation (Line(points={{-87,48.15},{-54.5,48.15},{-54.5,55.5},{-43.8,55.5}},
            color={255,128,0}));
     annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end PVElektrolyseur;
  end Examples;

  record H2Verbrauch
    parameter Modelica.SIunits.MolarFlowRate n_H2[:, 2] = [0,9;
  3600,9;
  7200,40;
  10800,30;
  14400,12;
  18000,280;
  21600,80;
  25200,281;
  28800,283;
  32400,80;
  36000,9;
  39600,9;
  43200,12;
  46800,12;
  50400,9;
  54000,9;
  57600,9;
  61200,9;
  64800,9;
  68400,9;
  72000,9;
  75600,12;
  79200,283;
  82800,282;
  86400,297]
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));

  end H2Verbrauch;

  function DataFromTable
    // Linear interpolation of data
    // diese Tabelle benoetigen Sie erst in Labor 6 und 7
    input Real u "input value (first column of table)";
    input Real table[:, :] "table to be interpolated";
    output Real y "interpolated input value (icol column of table)";
  protected
    Integer i;
    Integer n "number of rows of table";
    Real u1;
    Real u2;
    Real y1;
    Real y2;
    Real alpha1;
    Real alpha2;
  algorithm
    n := size(table, 1);
    if u <= table[1, 1] then
      y := table[1, 2];
    elseif u < table[n, 1] then
      i := 2;
      while i < n and u >= table[i, 1] loop
        i := i + 1;
      end while;
      i := i - 1;
      u1 := table[i, 1];
      y1 := table[i, 2];
      u2 := table[i + 1, 1];
      y2 := table[i + 1, 2];
      alpha1 := (u2 - u) / (u2 - u1);
      alpha2 := (u - u1) / (u2 - u1);
      y := alpha1 * y1 + alpha2 * y2; // einfach u=u1
    else
      y := table[n, 2];
    end if;
  end DataFromTable;
  annotation (uses(Modelica(version="3.2.2")));
end HydrogenComponents;

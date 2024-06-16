TITLE '3PX3 Technical Analysis'    																																									! The problem identification

COORDINATES cartesian2  																																												! Coordinate system, 1D,2D,3D, etc

{
Model Assumptions:
- No convection or radiation from air
- Snow is between freshly fallen & lightly packed
- No effects on heating from windspeed
- 2D model is extrapolatable to 3D
- 1 sqft of pavement under avg snowfall & weather conditions for January in Hamilton
- No radiative heat transfer back to the concrete from snow
}

VARIABLES
Temp(threshold=1e-6)																																														! Temperature (C)

SELECT
ngrid = 19     																																																			! Method controls

DEFINITIONS
! Dynamic
k																																																								! Thermal Conductivity (W/cm*K^-1)
e																																																								! Emissivity (unitless)
rho																																																							! Density (g/cm^3)
C_p																																																							! Heat capacity at constant pressure (J/g*K^-1)
init_temp																																																				! Inital temperature (C)
qvol																																																							! Volumetric heat generation (W/cm^3)

! Static
sigma = 5.67e-5																																																	! Stefan-Boltzmann Constant (W/cm^2*K^-1)
Q_m  = 334																																																			! Specific Heat of Fusion of Ice (J/g)
T_inf = -4																																																				! Temperature of outdoor air (C)
L_pavement = 30.48																																															! Length of pavement (cm)
t_pavement = 7.62																																																! Thickness of pavement (cm)
t_snow = 8.64																																																		! Thickness of snow block (cm)

! Sweep
P_system = 37																																																	! System power (W)

! Equations
qdot = -k*grad(Temp)																																														! Heat flux for conduction (W/cm^2)
E_direct  = (0 - T_inf)*area_integral(C_p*rho*L_pavement, "Snow")																									! Energy required for raising snow to 0 C (J)
E_latent = Q_m*area_integral(rho*L_pavement, "Snow")																														! Energy required for phase change from ice to water (J)
A_snow = 2*(L_pavement^2 + 2*t_snow*L_pavement)																															! Surface are of snow block (cm^2)

INITIAL VALUES
Temp = init_temp																																																! Temperature (C)

EQUATIONS        																																																	! PDE's, one for each variable
rho*C_p*dt(Temp) = qvol - div(qdot)																																								! Heat equation

BOUNDARIES       																																																! The domain definition
	REGION "Pavement"
	k = 0.0225																																																			! Thermal Conductivity (W/cm*K^-1)
	e = 0																																																					! Emissivity (unitless)	
	rho = 2.3																																																			! Density (g/cm^3)
	C_p = 0.88																																																			! Heat capacity at constant pressure (J/g*K^-1)
	init_temp = T_inf																																																! Inital temperature (C)

	qvol = P_system/area_integral(L_pavement, "Pavement")																												! Volumetric heat generation due to hydronic system (W/cm^3)
	
	START (-15.24, 0)
		load(Temp) = 0																																															! Assume insulated sides apart from contact surface due to poor conductivity of air & no convective effects
	LINE TO (15.24, 0)
	LINE TO (15.24, -7.62)
	LINE TO (-15.24, -7.62)
	LINE TO CLOSE
	
	REGION "Snow"
	k = 0.00045																																																		! Thermal Conductivity (W/cm*K^-1)
	e = 0.98																																																				! Emissivity (unitless)
	rho = 0.1																																																			! Density (g/cm^3)
	C_p = 2.090																																																		! Heat capacity at constant pressure (J/g*K^-1)
	init_temp = T_inf																																																! Inital temperature (C)	
	
	qvol = sigma*e*A_snow*((Temp - 273.15)^4 - (eval(Temp, 0, -3.81) - 273.15)^4)/area_integral(L_pavement, "Snow")		! Volumetric heat generation due to radiant heat transfer (W/cm^3)

	START (-15.24, 0)
		load(Temp) = 0																																															! Assume insulated sides apart from contact surface due to poor conductivity of air & no convective effects
	LINE TO (15.24, 0)
	LINE TO (15.24, 8.64)
	LINE TO (-15.54, 8.64)
	LINE TO CLOSE

TIME 0 TO 3600 halt (val(Temp, 0, 8.64) > 0)																																				! Run is finished when snow is melted (i.e. all snow is above 0 C)
PLOTS          	  																																																		! Save result displays
for t = 0 by endtime/60 to endtime
	history(val(Temp, 0, 8.64))
	history(val(Temp, 0, -3.81))
	contour(Temp) painted
	vector(qdot) norm
SUMMARY																																																				! Report energy usage
	report(t*P_system) as "Energy Consumed without Phase Change  "
	report(E_direct) as "Energy Required without Phase Change  "
	report(E_direct/(t*P_system)) as "Efficiency  "
	report(E_latent) as "Energy Required for Phase Change  "
	report(t*P_system*(1 + 1/E_direct*E_latent)) as "Predicted Energy Consumed with Phase Change  "
	report(t) as "Time to Finish  "
END

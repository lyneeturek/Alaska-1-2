%[INPUT] Dimensions of structure (ft)
dims=[1.125,1.125,1.125];

%[INPUT] Thermal values
T_inf=23.8927; %Outside air temperature (C)
T_air=32.3847; %Internal air temperature (C)
th=.75; %Wall thickness (in)
K=0.036; %Conductivity of insulation (W/mK)
windspeed=6.7; %Wind speed (m/s)
Q=125; %Heat generation (W)

%Constant calculations
dims=dims*unitsratio('meters','feet');
th_metric=th*unitsratio('meters','in');
A=2*dims(1)*dims(2)+2*dims(1)*dims(3)+2*dims(2)*dims(3);
roughness=[8.23,4.0,-0.057]; %Surface roughness coefficients (clear pine)
h=roughness(1)+roughness(2)*windspeed+roughness(3)*windspeed^2;

%Thermal resistances
R_conv=1/(h.*A); %Convection resistance (K/W)
R_cond_exp=(T_air-T_inf)/Q - R_conv; %Experimental Conduction Resistance (K/W)

R_cond_mod=th_metric/(K*A); %Model Conduction Resistance

shape_factor=R_cond_exp/R_cond_mod;




%{
    Heat Loss Estimation

    Input: Dimensions of the structure, thermal values
    Output: Maximum heat loss of structure in steady state
    
    By Cole Thomson 3/7/19
    Edited by Daniel Chan 3/9/19
%}

%[INPUT] Dimensions of structure (ft)
dims=[20,10,10];

%[INPUT] Thermal values
T_inf=-20; %Outside air temperature (C)
T_air=18; %Optimal cabbage-growing temperature (C)
th=6.5; %Wall thickness (in)
    %SIP with R-value of 24
K=1/25.6005; %Conductivity of insulation (W/mK)
windspeed=4; %Wind speed (m/s) - cut in wind speed for turbines
F_s=0.7145; %Shape factor for conduction

%Constant calculations
dims=dims*unitsratio('meters','feet');
th_metric=th*unitsratio('meters','in');
A=2*dims(1)*dims(2)+2*dims(1)*dims(3)+2*dims(2)*dims(3);
roughness=[8.23,4.0,-0.057]; %Surface roughness coefficients (clear pine)
h=roughness(1)+roughness(2)*windspeed+roughness(3)*windspeed^2;

%Thermal resistances
R_cond=th_metric/(K.*A)*F_s;
R_conv=1/(h.*A);

q_ss = (T_air-T_inf)./(R_cond+R_conv)


%% Transient - Heating up a 10x10x20 container full of air from 0C to 18C in 30 minutes
vol = 56.63; %Volume of the space in m^3
density = 1.225; %Density of air in kg/m^3
m = vol * density; 
c = .718 %kJ/kg*K
T1 = 0;
T2 = 18;

E1 = m*c*T1; %in kJ - not correct in an absolute sense because T is not in Kelvin, but good for difference in E
E2 = m*c*T2; %in kJ

t = .5*60*60; %half an hour in seconds

q_trans = (E2-E1) / t * 1000; %Energy/time = power, *1000 to get units to Watts

%% Total
q_total = q_ss + q_trans

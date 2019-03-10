%{
    Input: Maximum heat loss of system, Minimum difference between temp
    of interior and temp of pipes

    Output: Q/degrees C needed, a metric of effectiveness of heat tranfer
    of a pump
%}

%% Steady State 

Q_out_max = 1200; %W
T_difference_min = 50-18; %The worst case temperature differential would be the water at the minimum temperature and the air inside at 18C.

Q_over_T_SS = Q_out_max / T_difference_min

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

W = (E2-E1) / t * 1000 %Energy/time = power, *1000 to get units to Watts

Q_over_T_Transient = W/T_difference_min %Additional load on the heater for the same worst case Temperature differential

%%
Q_over_T_Total = Q_over_T_SS + Q_over_T_Transient 
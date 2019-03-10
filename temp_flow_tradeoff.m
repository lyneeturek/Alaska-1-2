%{

Function: Shows the difference minimum tank temperatures and maximum flow rates
needed to achieve sufficient heating power
    
    Input: Worst case heat loss

    By Cole Thomson 3/7/19
    Edited by Daniel Chan 3/9/19

%}

%User input
q=1153; %worst-case heat loss (W)
t_end=48.5; %water temperature at end of heating loop (C)
    %ASSUMPTION: water equilibrates with greenhouse inner lower-bound temp
t_start=50; %water temperature at start of heating loop (C)


%Constants
cp=4181.5; %specific heat of water at avg temp of 50 C (J/kg K)
mdot=q/cp*1./(t_start-t_end); %mass flow rate (kg/s)
vol_flow_rate_imperial=mdot*15850.3231/1000; %mass flow rate in gpm
plot(t_start,mdot_imperial);
ylabel('Flow Rate (gpm)');
xlabel('Minimum Tank Temp (C)');
title('Tradeoff Between Flow Rate and Tank Temperature');
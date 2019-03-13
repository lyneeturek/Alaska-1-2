%{

    Function: Shows the difference minimum tank temperatures and maximum flow rates
    needed to achieve sufficient heating power
    
    Input: Pipe length, diameter, water input temperature, flow rate
    Output: Final temperature, final heat transfer

    By Lewis Buenrostro 3/1/19
    Edited by Daniel Chan 3/9/19

%}

clear; close all; clc;

%INPUTS
pipe_length = 8; % m - 40 feet
m_dot = 0.189; % kg / s - 3 gal/min


%CONSTANTS
D = 0.01905; % m - 0.75 in
rho = 1000; % kg / m^3
C_p = 4200; % 
h = 5; %Assumption for some fanning (THe h value calculated above is wayyy too high)
T_in = 90; % degrees Celsius
T_air = 16; % degrees Celsius - ~60 degrees Fahrenheit

T(1) = T_in;
Q_dot(1) = m_dot * C_p * T(1);
steps = 500;
step_length = pipe_length / steps; % m
A_s = pi * D * step_length; % m^2
Q_out_total = 0;

for i = 1:steps
    Q_out = A_s * h * (T(i) - T_air);
    Q_dot(i+1) = (Q_dot(i) - Q_out);
    T(i+1) = Q_dot(i+1) / m_dot / C_p;
    Q_out_total = Q_out_total + Q_out;
end

Q_out_total
T(end);
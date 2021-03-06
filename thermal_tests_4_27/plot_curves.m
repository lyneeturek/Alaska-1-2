clear all; close all; clc
load('temperature_data_4_27');

external=temps(:,1);
long_axis=temps(:,2);
short_axis=temps(:,3);
time=transpose(1:length(external));

adjusted_long = long_axis - external;
adjusted_short = short_axis - external;

hold on;
%plot(time,external);
plot(time,adjusted_short);
plot(time,adjusted_long);

legend('external','internal(short axis)','internal(long axis)');

lb=400;
hb=455./60;
time = time./60;
smoothed_diff = movmean(adjusted_short-adjusted_long, 10);
figure(2);
plot(time,smoothed_diff, 'LineWidth', 2);
ylim([-3 3])

ylabel 'Degrees C'
xlabel 'Time (Hours)'
title 'Temperature Difference of Internal Thermocouples vs Time'


T_inf = mean(external(400:455));
T_ss = (mean(long_axis(400:455))+mean(short_axis(400:455)))/2;

clear all

%water drum; https://www.uline.com/Product/Detail/S-14364/Drums/Steel-Drum-with-Lid-30-Gallon-Open-Top-Unlined?pricode=WB0299&gadtype=pla&id=S-14364&gclid=Cj0KCQiA7briBRD7ARIsABhX8aBXzahTKj4O7YuUj8KMNm21OF2ZJds-dQyWbdpo2vGf0dmV7Cs9M0kaAoHpEALw_wcB&gclsrc=aw.ds
%allowable drum temperature range: 233.15K to 435.928K (-40F to 325F)
t_d = 1.1E-3; %[m]
k_d = 385.0; %[W/mK]
k_ins = 0.04; %[W/mK] fiberglass http://hyperphysics.phy-astr.gsu.edu/hbase/Tables/thrcn.html
t_ins = .152; %[m]
ar_d = .633; %aspect ratio drum

%Water properties
rho_w = 963.33; %kg/m^3 https://water.usgs.gov/edu/density.html
%c_w = 3.8204*1000; % [kJ/kgK]*[J/kJ] = [J/kgK] Isochoric heat capacity of water at 363K https://www.engineeringtoolbox.com/specific-heat-capacity-water-d_660.html
c_w = 4.19 * 1000;
k_w = .606;

%flow rate
conv = 0.00378541; %1 gal = 0.00378541 m^3
%[m^3/s][60s/min][1gal/.0038m^3]
%Temperatures
T_outside = -23 + 273; %[K] coldest temperatur Kong could reach
T_w_0 = 95 + 273; %[K] intial water temperature
T_gh = 18 + 273; %[K] greenhouse temperature

%time assumptions
t_req = 4 * 24 * 3600; %seconds in a week
dt = 1; %1 second
t_storage = zeros(1, 100);

%Greenhouse assumptions
V_side1 = 3; %[m]
V_side2 = 3; %[m]
V_floor = 6; %[m]
V_house = V_side1 * V_side2 * V_floor; %[m^3 = 10ft^3]
%http://highperformanceinsulation.eu/wp-content/uploads/2016/08/Thermal_insulation_materials_made_of_rigid_polyurethane_foam.pdf
k_ins_wall = 0.025; %Rigid polyurethane foam (PUR/PIR) insulation; [(W/(mK)])
rho_ins = 30; %density [kg/m^3]
A_ins_2 = V_side1 * V_side2; %insulation Area [m^2]
A_ins_1 = V_side1 * V_floor;
t_ins_wall = 0.1524; %thickness; [m]
t_ins_floor = .0254;

Qdot_gh_1 = 2 * k_ins_wall * A_ins_2 * (T_gh - T_outside) / t_ins_wall; %Heat transfered per second from gh to environment; [W]
Qdot_gh_2 = 2 * k_ins_wall * A_ins_1 * (T_gh - T_outside) / t_ins_wall;
Qdot_gh_3 = 2* k_ins_wall * A_ins_1 * (T_gh - T_outside) / t_ins_floor;
m_in = (Qdot_gh_1 + Qdot_gh_2 + Qdot_gh_3) / (c_w * (T_w_0 - T_gh)); %mass required for a dt by the greenhouse
    
n_nodes = 10;
gal = linspace(1,1000,1000);
T_old = zeros(1,n_nodes);
[i1, i2] = size(gal);

for j = 1:i2
    
T_w = T_w_0;
 
t_s = 0;  

T_old(1:n_nodes) = T_w;
T = T_old;

%new amount of water (step up)
V_w = gal(j) * conv; %[m^3]; 
h_d = (4 * V_w / (pi * ar_d^2))^(1/3); % height drum [m]
r_d = ar_d * h_d / 2; %radius drum [m]
m_w = V_w * rho_w; %[kg] mass water

m_n = m_w/n_nodes; %[kg] mass of a node
h_n = h_d/n_nodes; %[m] height of a node

    while T(n_nodes) > 60 + 273
        
     t_s = t_s + 1; %time step up one
     
    %first node 
    T(1)  = ((m_n - m_in) * T_old(1) + m_in * T_gh) / m_n;

    %mass displacement into other nodes
        for i = 2:n_nodes
            T(i) = ((m_n - m_in) * T_old(i) + m_in * T_old(i-1))/ m_n;
        end
    
    %heat transfer out first node
    %heat out wall
    R_d = log((r_d + t_d)/r_d) / (2 * pi * k_d * h_n); %thermal resistance through wall
    R_ins = log((r_d + t_d + t_ins)/(r_d + t_d)) / (2 * pi * k_ins * h_n); %thermal resistance through drum insulation
    Qdot_wall = (T(1) - T_gh) / (R_d + R_ins); %[W]
    Q_wall = Qdot_wall * dt; %energy loss for the time period [J]

    %heat out bottom
    R_d = t_d / (k_d * pi * r_d^2); %bottom drum
    R_ins = t_ins / (k_ins * pi * r_d^2); %bottom insulation
    Qdot_bottom = (T(1) - T_gh) / (R_d + R_ins);
    Q_bottom = Qdot_bottom * dt;

    Q_out = Q_bottom + Q_wall;

    %heat in from top node
    Q_in = k_w * pi * r_d^2 * (T(2) - T(1)) / h_n;
    
    %new T
    Q_tot = dt * (-Q_out + Q_in);
    T(1) = Q_tot / (m_w * c_w) + T(1);

    for i = 2 : n_nodes -1

    %heat out wall
    R_d = log((r_d + t_d)/r_d) / (2 * pi * k_d * h_n); %thermal resistance through drum
    R_ins = log((r_d + t_d + t_ins)/(r_d + t_d)) / (2 * pi * k_ins * h_n); %thermal resistance through drum insulation
    Qdot_wall = (T(i) - T_gh) / (R_d + R_ins); %[W]
    Q_wall = Qdot_wall * dt; %energy loss for the time period [J]

    %heat out to bottom node
    Qdot_bottom = k_w * pi * r_d^2 * (T(i) - T(i-1));
    Q_bottom = Qdot_bottom * dt;

    Q_out = Q_bottom + Q_wall;

    %heat in from top node
    Q_in = k_w * pi * r_d^2 * (T(i+1) - T(i)) / h_n;
    
    Q_tot = dt*(Q_in - Q_out);
    T(i) = Q_tot / (m_w * c_w) + T(i);

    end

    %heat transfer out last node
    %heat out wall
    R_d = log((r_d + t_d)/r_d) / (2 * pi * k_d * h_n); %thermal resistance through drum
    R_ins = log((r_d + t_d + t_ins)/(r_d + t_d)) / (2 * pi * k_ins * h_n); %thermal resistance through drum insulation
    Qdot_wall = (T(n_nodes) - T_gh) / (R_d + R_ins); %[W]
    Q_wall = Qdot_wall * dt; %energy loss for the time period [J]

    %heat out top
    R_d = t_d / (k_d * pi * r_d^2); %bottom drum
   R_ins = t_ins / (k_ins * pi * r_d^2); %bottom insulation
    Qdot_top = (T(n_nodes) - T_gh) / (R_d + R_ins);
    Q_top = Qdot_top * dt;

    %heat out from below
    Q_down = k_w * pi * r_d^2 * (T(n_nodes) - T(n_nodes - 1)) / h_n;

    Q_out = Q_down + Q_wall + Q_top;

    Q_tot = dt*(-Q_out);
    T(n_nodes) = Q_tot / (m_w * c_w) + T(n_nodes);
T_old = T;
    end
    
    t_storage(j) = t_s;
end

plot(gal, t_storage/24/3600 , 'LineWidth',2)
xlabel('Volume [gallons]')
ylabel('Time Elapsed [days]')
title('Heating periods of varying water tank sizes')



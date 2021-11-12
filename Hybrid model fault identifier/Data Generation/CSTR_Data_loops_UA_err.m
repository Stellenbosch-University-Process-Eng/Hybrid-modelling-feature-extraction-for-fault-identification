%% Loading and generating training and test datasets for SVM

% This code is used to generate the SVM datasets used by the multiclass SVM
% fault detection method

clear
clc
load CSTR_In % Loading CSTR parameters
file = 'CSTR_gen_UA_err';
%load k0Data
%load CaData
t_samp = 1;
dur = 4000;
%dur = 3000; 
%I = 100
%% Setting parameters

t_Tj = 0;
t_fin = 0;
time.fin = 10;
time.serr = 50;
time.test = 80;
time.k0 = 15;
time.test2 = 40;
time.Ca_i_s = 15;

Tj_i = 290;
fin = 290;
serr_drift = 0.005;
serr = 0;

% Setting error magnitudes 
Ca_i_mag = -0.0005;
k_0_mag = -0.008;

%% Normal data trainig set

Tj_o = Tj_i;
t_serr = 0;


k0Switch = 0;
CaSwitch = 0;
UASwitch = 0;

sim (file);

% Storing Simulink results
Ca_sens_nrm = Ca_sensor.Data;
Ca_out_nrm = Ca_out.Data;
F_in_nrm = F_in.Data;
T_in_nrm = Ti_in.Data;
Tj_sens_nrm = Tj_sensor.Data;
T_out_nrm = T_out.Data;
F_crt_nrm = F_crt.Data;
Ca_in_nrm = Ca_in.Data;
k_nrm = k_rate.Data;

%% k0 error training set

t_serr = 0;
k_0_s = k_0_mag; %1e9 : 2e8 : 9e9
k0Switch = 1;
n = 1;
for mag =  4.5e10        %5e9 : 1e9 : 4.5e10
    [NOC_k0(n,:), t_k0] = k0Data_fxn(mag,dur);
    sim (file);

    Ca_sens_fin(:,n) = Ca_sensor.Data;
    Ca_out_fin(:,n) = Ca_out.Data;
    F_in_fin(:,n) = F_in.Data;
    T_in_fin(:,n) = Ti_in.Data;
    Tj_sens_fin(:,n) = Tj_sensor.Data;
    T_out_fin(:,n) = T_out.Data;
    F_crt_fin(:,n) = F_crt.Data;
    Ca_in_fin(:,n) = Ca_in.Data;
    k_fin(:,n) = k_rate.Data;
    n = n+1;
end
%% Ca_i error training set

k_0_s = 0;
Ca_i_s = Ca_i_mag;
k0Switch = 0;
%CaSwitch = 1;
UASwitch = 1;
n= 1;
%26 22:1:26 50:1:55
for mag = 75:1:80        %24:1:26
    [NOC_UA(n,:), t_Ca] = UAData_fxn(mag,dur);
    sim (file);

    Ca_sens_UA(:,n) = Ca_sensor.Data;
    Ca_out_UA(:,n) = Ca_out.Data;
    F_in_UA(:,n) = F_in.Data;
    T_in_UA(:,n) = Ti_in.Data;
    Tj_sens_UA(:,n) = Tj_sensor.Data;
    T_out_UA(:,n) = T_out.Data;
    F_crt_UA(:,n) = F_crt.Data;
    Ca_in_UA(:,n) = Ca_in.Data;
    k_UA(:,n) = k_rate.Data;
    UAp_UA(:,n) = UA_rate.Data;
    n = n+1;
end

save ('CSTR_gen_data_UA_2')
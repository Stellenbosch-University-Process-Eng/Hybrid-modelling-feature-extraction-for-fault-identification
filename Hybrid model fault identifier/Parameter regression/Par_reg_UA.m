%% Online parameter regression
%  This code will perform online regression to estimate the parameters 
%  k, ??? and UA associated with a jacketed CSTR
%  The unknown parameters k, dH and UA are contained in a parameter vector "p"
%  The remaining fixed parameters are contained in a parameter structure "pf"
%  Any inputs varying with time are contained in the function structure "f"
%% Initialize
clc 
clear 
clf

%  Load data from file
load 'collection 10-08'
load 'CSTR_gen_data_UA_2'
tspan = (0:1:dur)'; % Total data timespan  ???  Why is this not in the saved file ???

%  Guess initial parameter values
%  p(1) = k   (1/s),    reaction rate constant
%  p(2) = ???
%  p(3) = UA  (W/K),    overall heat transfer coefficient x heat transfer area

%pg = [0.3, 4000, 3540];
%pg = [0.05,49500/(237*100)];
% 49500/(237)
pg = [0.9,0.9];

% Define structure for fixed parameters
pf.V = V; 
pf.Ca_i = Ca_i; 
pf.HD = HD;
pf.trans = [0.03333, 2.32444 ]; %232.45

%% Perform online regression
options = optimoptions('lsqnonlin','display','none');   % Suppress regression output

% each column represents a diffrent magnitude of the fault, for loop goes
% through the fault magnitudes sequentially
for l = 4:1:6
    pm_UA = zeros(length(tspan),3);

    % Define interpolants for input functions for the given dataset
    f.Fest =   @(t) interp1(tspan, F_crt_UA(:,l),  t);   % (m3/s), inlet flowrate
    f.Ti_est = @(t) interp1(tspan, T_in_UA(:,l),   t);   % (K), inlet temperature
    f.Tj_est = @(t) interp1(tspan, Tj_sens_UA(:,l),t);   % (K), jacket temperature
    
    % Perform regression over a sliding window of data
    WindowLength = 16;
    for i = 1:1:(length(tspan)-WindowLength)
        Index   = (0:WindowLength) + i;   % Index of data corresponding to the current window
        Window.ts = tspan(Index);           % Time points for current window
        Window.Y  = [Ca_sens_UA(Index, l), T_out_UA(Index, l)]; % Measured data for current window
        
        % Regress parameters and store in array
        p = lsqnonlin(@(p) LSQ_int(p, Window, pf, f), pg, [0,0],[1,1], options);
        pm_UA(i,1) = p(1)*pf.trans(1);
        pm_UA(i,2) = p(2)*pf.trans(2);
        %pm_UA(i,3) = p(3)*pf.trans(3);
        pg = p; % Update initial guess
        
        % Track progress
        fprintf(' l = %d / 3, i = %d / %d \n', l, i, length(tspan)-20);
        
        % Simulate using current parameters and plot prediction vs data
        %The following code is not necessary but useful for tracking progress
        Ypred = funcSimulate(p.*pf.trans, Window, pf, f);
        subplot(2,1,1)
        plot(Window.ts, Ypred(:,1), Window.ts, Window.Y(:,1),'x:');
        subplot(2,1,2)
        plot(Window.ts, Ypred(:,2), Window.ts, Window.Y(:,2),'x:');
        drawnow
               
        
    end
    pm_UA_col(:,:,l) = pm_UA;
end

%load ('collection 10-08') 
pg = [0.9 0.9 ];

for l = 1
    pm_fin = zeros(length(tspan),3);

    % Define interpolants for input functions for the given dataset
    f.Fest =   @(t) interp1(tspan, F_crt_fin(:,l),  t);   % (m3/s), inlet flowrate
    f.Ti_est = @(t) interp1(tspan, T_in_fin(:,l),   t);   % (K), inlet temperature
    f.Tj_est = @(t) interp1(tspan, Tj_sens_fin(:,l),t);   % (K), jacket temperature
    
    % Perform regression over a sliding window of data
    WindowLength = 16;
    for i = 1:1:(length(tspan)-WindowLength)
        Index   = (0:WindowLength) + i;   % Index of data corresponding to the current window
        Window.ts = tspan(Index);           % Time points for current window
        Window.Y  = [Ca_sens_fin(Index, l), T_out_fin(Index, l)]; % Measured data for current window
        
        % Regress parameters and store in array
        p = lsqnonlin(@(p) LSQ_int(p, Window, pf, f), pg, [0,0,0],[1,1,1], options);
        pm_fin(i,1) = p(1)*pf.trans(1);
        pm_fin(i,2) = p(2)*pf.trans(2);
        %pm_fin(i,3)  = p(3)*pf.trans(3);
        pg = p; % Update initial guess
        
        % Track progress
        fprintf(' l = %d / 41, i = %d / %d \n', l, i, length(tspan)-20);
        
        % Simulate using current parameters and plot prediction vs data
        %The following code is not necessary but useful for tracking progress
        Ypred = funcSimulate(p.*pf.trans, Window, pf, f);
        subplot(2,1,1)
        plot(Window.ts, Ypred(:,1), Window.ts, Window.Y(:,1),'x:');
        subplot(2,1,2)
        plot(Window.ts, Ypred(:,2), Window.ts, Window.Y(:,2),'x:');
        drawnow
               
        
    end
    pm_fin_col_2(:,:,41) = pm_fin;
end   

for l = 41
    pm_serr = zeros(length(tspan),3);

    % Define interpolants for input functions for the given dataset
    f.Fest =   @(t) interp1(tspan, F_crt_serr(1:4001,l),  t);   % (m3/s), inlet flowrate
    f.Ti_est = @(t) interp1(tspan, T_in_serr(1:4001,l),   t);   % (K), inlet temperature
    f.Tj_est = @(t) interp1(tspan, Tj_sens_serr(1:4001,l),t);   % (K), jacket temperature
    
    % Perform regression over a sliding window of data
    WindowLength = 16;
    for i = 1:1:(length(tspan)-WindowLength)
        Index   = (0:WindowLength) + i;   % Index of data corresponding to the current window
        Window.ts = tspan(Index);           % Time points for current window
        Window.Y  = [Ca_sens_serr(Index, l), T_out_serr(Index, l)]; % Measured data for current window
        
        % Regress parameters and store in array
        p = lsqnonlin(@(p) LSQ_int(p, Window, pf, f), pg, [0,0],[1,1], options);
        pm_serr(i,1) = p(1)*pf.trans(1);
        pm_serr(i,2) = p(2)*pf.trans(2);
        %pm_fin(i,3)  = p(3)*pf.trans(3);
        pg = p; % Update initial guess
        
        % Track progress
        fprintf(' l = %d / 41, i = %d / %d \n', l, i, length(tspan)-20);
        
        % Simulate using current parameters and plot prediction vs data
        %The following code is not necessary but useful for tracking progress
        Ypred = funcSimulate(p.*pf.trans, Window, pf, f);
        subplot(2,1,1)
        plot(Window.ts, Ypred(:,1), Window.ts, Window.Y(:,1),'x:');
        subplot(2,1,2)
        plot(Window.ts, Ypred(:,2), Window.ts, Window.Y(:,2),'x:');
        drawnow
               
        
    end
    pm_serr_col_2(:,:,41) = pm_serr;
end   
save ('UA_regressed_par_5')

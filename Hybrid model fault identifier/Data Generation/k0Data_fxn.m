function [NOC_vec, time] = k0Data_fxn (mag,dur)

% Use this script to generate a mat-file that can be read by Simulink
% mag - fault magnitude
% dur - fault average duration

% Normal operating conditions (NOC)
NOC_k0 = 7.2e10;         % Value of k0 under NOC, default 1
NOC_T = 600;       % Average duration of NOC (i.e. time between faults) default 1000 for T = 3000 used 600, 60
NOC_T_dev = 150;    % Standard deviation of NOC duration (i.e. deviation of time between faults) default 200 used 150, 15

% Faulty conditions)
%fault_slope = 8e-2; % Rate at which k0 degrades under fault condition (per second) default 1e-2
fault_slope = mag;   % setting fault ramp as fault magnitude
fault_T = 80;       % Average fault duration default 50 , for T = 3000, used 80, 11
fault_T_dev = 25;   % Standard deviation of fault duration default 20 used 25, 2.5

% Initialize values
T = dur;                               % Simulation time
time = linspace(0,T,1e4);               % Time vector
dt = time(2) - time(1);                 % Time-step
k0 = zeros(1, length(time));            % Allocate memory for a k0 vector
NOC_vec = zeros(1, length(time));       % Allocate memory for a vector to track NOC/faulty conditions
NOC = 1;                                % NOC = 1 if NOC, NOC = 0 if faulty
time_since_switch = 0;                  % Duration of current condition

% Loop over all timesteps
for i = 1 : length(time)
    time_since_switch = time_since_switch + dt;     % Increase duration of current condition
    r = rand;                                       % Generate random number to determine if switch occurs
    if NOC
        k0(i) = NOC_k0;  
        NOC_vec(i) = 1;
        if r < cdf('Normal',time_since_switch, NOC_T, NOC_T_dev)
            % Switch if the random number between 0 and 1 is less than the
            % normal cumulative distribution function with mean "NOC_T" and
            % standard deviation "NOC_T_dev", evaluated at "time_since_switch"
            NOC = 0;
            time_since_switch = 0;
        end
    else
        k0(i) = k0(i-1) - fault_slope*dt;
        NOC_vec(i) = 0;
        if r < cdf('Normal',time_since_switch, fault_T, fault_T_dev)
            % As above
            NOC = 1;
            time_since_switch = 0;
        end
    end
end
InputData = [time; k0];
% storing results so that it can be read by the simulink model
save('k0Data', 'InputData');
end
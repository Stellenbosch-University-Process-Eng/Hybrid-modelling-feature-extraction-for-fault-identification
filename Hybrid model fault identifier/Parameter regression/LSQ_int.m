function Err = LSQ_int(p, Window, pf, f)
% Function to estimate the parameters given below through regression
%  p(1) = k   (1/s),   reaction rate constant
%  p(2) = ???
%  p(3) = UA  (W/K),   overall heat transfer coefficient x heat transfer area
%  V          (m3),     reactor volume
%  Ca_i       (mol/m3), inlet concentration
%  HD         ????
%  f.Fest     (m3/s),   inlet flowrate
%  f.Ti_est   (K),      inlet temperature
%  f.Tj       (K),      jacket temperature

% converting intially parameters guess to the actual parameter values
p_trns = p;
p_trns(1) = p_trns(1)*pf.trans(1);
p_trns(2)= p_trns(2)*pf.trans(2);
%p_trns(3)= p_trns(3)*pf.trans(3);
% Simulate the system using guessed values of "p"
Ypred = funcSimulate(p_trns, Window, pf, f);

% Compare predicted values (given below) to measured values. Errors are
% scaled by the mean value of each measurement
%  Y(:,1) = Ca (mol/m3), outlet concentration
%  Y(:,2) = T  (K),      reactor temperature
Ca_err = ( Ypred(:,1) - Window.Y(:, 1) )   ./ mean( Window.Y(:,1) );
T_err  = ( Ypred(:,2) - Window.Y(:, 2) )   ./ mean( Window.Y(:,2) );
%Err = Ca_err + T_err;
Err=[Ca_err;T_err];
end
        
 
 

   

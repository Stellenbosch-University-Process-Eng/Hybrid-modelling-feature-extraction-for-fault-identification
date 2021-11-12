function Ypred = funcSimulate(p, Window, pf, f)
% Function to simulate the CSTR from initial values "y0" over timespan "ts"
%
%  p(1) = k   (1/s),    reaction rate constant
%  p(2) = ???
%  p(3) = UA  (W/K),    overall heat transfer coefficient x heat transfer area
%  V          (m3),     reactor volume
%  Ca_i       (mol/m3), inlet concentration
%  HD         ????
%  f.Fest     (m3/s),   inlet flowrate
%  f.Ti_est   (K),      inlet temperature
%  f.Tj       (K),      jacket temperature


model = @(t,y) [(f.Fest(t)/pf.V) * ( pf.Ca_i        - y(1) )...
                      - p(1)*y(1); ...
                (f.Fest(t)/pf.V) * ( f.Ti_est(t) - y(2) ) ...
                      + (50000/(pf.HD))*p(1)*y(1)...    %(50000/(pf.HD))*p(1)*y(1)
                      - p(2)*(y(2) - f.Tj_est(t)) ];

[~, Ypred] = ode45(model, Window.ts, Window.Y(1,:));        

end
        
 
 

   
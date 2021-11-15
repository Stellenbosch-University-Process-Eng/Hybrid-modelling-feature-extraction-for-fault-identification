function [res, Ypred,par_est] = res_est (X,Y,par,bsize,Nblock,par_guess,ncomp,const,time,lambda)

F_crt = X(:,2);
Tj_sens = X(:,1);
T_in = X(:,3);
mu = 0;

Fest = @(t) interp1(time,F_crt,t);
Ti_est = @(t) interp1(time,T_in,t);
Tj_est = @(t) interp1(time,Tj_sens,t);

Ca_i = const(1);
EaR = const(2);
V = const (3);
HA = const (4);
HD = const(5);
H = const(6);

%par_est(1:50,1) = par_guess
%par_est(1:bsize,1) = par_guess;
par_est(1:bsize,1) = par(1:bsize);
% First set of residual estimates 

 for k = 1:1:bsize
        
        model = @(t,y) [(Fest(t)/V)*(Ca_i-y(1))-par_est(k)*exp(-EaR/(y(2)))*y(1);
       (Fest(t)/(V))*(Ti_est(t)-y(2))+H/(HD)*par_est(k)*exp(-EaR/(y(2)))*y(1)-HA/(HD*V)*(y(2)-Tj_est(t))];
        
    [~, yp] = ode45(model,[time(k) time(k+1)],[Y(k,1) Y(k,2)]);
   
    Ypred(k,:) = yp(end,:);    
        
 end

res(1:bsize,:) = Ypred(1:bsize,:) - [Y(1:bsize,1) Y(1:bsize,1)];
 
for i = 1:1:Nblock
    %[~,~,b,~,P,Q,W] = npls([X(1+bsize*(i-1):bsize*i,:) res(1+bsize*(i-1):bsize*i,:)],par_est(1+bsize*(i-1):bsize*i,1), ncomp);
    [~,~,b,~,P,Q,W] = npls([X(1+bsize*(i-1):bsize*i,:) res(1+bsize*(i-1):bsize*i,:)],par(1+bsize*(i-1):bsize*i,1), ncomp);
    
    B = diag(b);
    
    
    if i>1
        
    X0 = [lambda*Pold';P'];
    Y0 = [lambda*Bold*Qold';B*Q'];
    
    
    [~,~,b,~,P,Q,W] = npls(X0,Y0,ncomp);
    
    B = diag(b);
    end
    
    Bold = B;
    Pold = P;
    Qold = Q;
    
    if i < Nblock
    
        strt = 1+bsize*i;
        fin = bsize*(i+1);
        rs = bsize*(i-1)+1;
        rf = i*bsize;
        
    elseif i == Nblock
        
        strt = 1+bsize*i;
        fin = size(par,1);
        rs = bsize*(i-1)+1;
        rf = bsize*(i-1)+1+fin-strt;
    end
    
    [tscale,mu] = xscale([X(strt:fin,:) res(rs:rf,:)],W,P,ncomp,mu);
    	
    par_hold = tscale.*b.*Q;
    
    par_est(strt:fin,:)= sum(par_hold,2);
    par_est(strt:fin,:)= par_est(strt:fin,:)+mean(par);
    %par_est(strt:fin,:)= par_est(strt:fin,:)+mean(par_est(1+bsize*(i-1):bsize*i,1));
    
    for k = strt:1:fin
        
        model = @(t,y) [(Fest(t)/V)*(Ca_i-y(1))-par_est(k)*exp(-EaR/(y(2)))*y(1);
       (Fest(t)/(V))*(Ti_est(t)-y(2))+H/(HD)*par_est(k)*exp(-EaR/(y(2)))*y(1)-HA/(HD*V)*(y(2)-Tj_est(t))];
        
    [~, yp] = ode45(model,[time(k-1) time(k)],[Y(k-1,1) Y(k-1,1)]);
   
    Ypred(k,:) = yp(end,:);    
        
    end
    
    res(strt:fin,:) = Ypred(strt:fin,:) - Y(strt:fin,:); 
    
end

%% Scales X and determines Xscores for testing/ prediction input set as defined in the PLS algorithm provided by S.Wold
function [Tscale,mu] =  xscale (X,W,P,ncomp,mean_prev)

if size(X,1) == 1
    X1 = X - mean_prev;
    mu = 0;
else
X1 = X - mean(X);
mu = mean(X);
end

for j = 1:1:ncomp
    
    Tscale(:,j) = X1*W(:,j);
    X1 = X1-Tscale(:,j)*P(j);
end

end
end
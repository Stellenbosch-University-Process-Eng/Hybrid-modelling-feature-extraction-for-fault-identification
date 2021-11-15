function [T,U,B,mu,P,Q,W,Beta] = npls (X,Y, ncomp)

% PLS regression using NIPALS algorithm
% C represents Q from the Wold (Wold et al, 1987) description of the 
% algorithm

% X - Training input used for regression
% Y - Training outputs which PLS model should predict
% ncomp - number of PLS components
%% Preparing input and output data


X0 = X - mean(X); % centering matrices
Y0 = Y - mean(Y);

mu = std(Y0);

X0 = X0./ repmat(std(X0), size(X,1), 1); % scaling to unit variance
Y0 = Y0./ repmat(std(Y0), size(Y,1), 1);

% Centering and scaling matrices handled by zscores function

%X0 = zscore(X);
%Y0 = zscore(Y);

Xpls = X0;
Ypls = Y0;
u = Ypls(:,1); % setting initial estimate of u

%% NIPALS algorithm

% for loop run to generate parameters for desired number of components
for i=1:ncomp 
    
    maxi = 1; % initialising max iteration counter
    conv(1) = 10; %intialising convergence tracker
    %iteration criteria, loop will run until either 100 iterations have
    %occured or convergence is less than 0.01
    while (maxi <= 100 && conv > 0.01) 
        w = Xpls'*u/(u'*u);
        w = w/norm(w);
        uold = u;
        t = Xpls*w/(w'*w);
        c = Ypls'*t/(t'*t);
        c = c/norm(c);
        u = Ypls*c/(c'*c);
        conv = norm(u-uold)/norm(u);
        p = Xpls'*t/(t'*t);
        maxi = maxi+1;
    end
    q = Ypls'*u/(u'*u);
    
    % Additional conversions required by Wold NIPALS algorithm to ensure
    % predictions are correct
    t = t*norm(p');
    wp = w'*norm(p');
    w = wp';
    p = p/norm(p);
    
    % Storing PLS parameters
    b = (u'*t)/(t'*t);
    Q(:,i) = q;
    W(:,i) = w;
    T(:,i) = t;
    P(:,i) = p;
    C(:,i) = c;
    U(:,i) = u;
    B(i)= b;

    % Residuals calculation which will be used as X and Y for next PLS
    % component
    Xpls = Xpls - t*p';
    Ypls = Ypls - b*t*c';
end

Beta = W*(P'*W)^(-1)*C';

% Feeds parameter predictions to SVM model
%clear
%load 'SVM_break'

%% Obtaining Parameter estimates


Ca_in_pg = [Ca_in_nrm;Ca_in_fin;Ca_in_serr];
Tj_sens_pg = [Tj_sens_nrm;Tj_sens_fin;Tj_sens_serr];
Ca_sens_pg = [Ca_sens_nrm;Ca_sens_fin;Ca_sens_serr];
F_crt_pg = [F_crt_nrm;F_crt_fin;F_crt_serr];
T_in_pg = [T_in_nrm;T_in_fin;T_in_serr];
T_out_pg = [T_out_nrm; T_out_fin; T_out_serr];
In = [Tj_sens_pg,F_crt_pg, T_in_pg];
Out = [Ca_sens_pg, T_out_pg];

pg = [100, 4000, 3540];

[par] = svm_hybrid_v2(In,Out,pg,V,Ca_i,HD);

%% Finding optimum forgetting factor (lambda)

Dataset = [In,Out];

CVdata = cvpartition(par(:,1),'HoldOut');
X = Dataset(CVdata.training,:);
Y = par(CVdata.training,1);
Xtest = Dataset(CVdata.test,:);
Ytest = par(CVdata.test,1);
ncomp = 2;
%lambda = 0.95;
bsize = 50;
mu = 0;
x0 = 0.8;

[lambda, fval, exitflag, output] = fminunc(@(x) obj_recurs(X, Y, Xtest, Ytest, ncomp, x, bsize),x0);

%% Recursive PLS and predictions 

[b,Q,P,W,T,U,bcat,Qcat,Yp,Wcat,Pcat] = recursive_pls(X,Y,ncomp,lambda,bsize);

Xtest = Dataset(CVdata.test,:);
Ytest = par(CVdata.test,1);
ts = xscale(Xtest,W,P,ncomp);

for i=1:1:size(ts,1)
Yts(i,:) = ts(i,:).*b.*Q;
end

Yt = sum(Yts,2);
Yt = Yt+mean(Y);
N = floor(size(X,1)/bsize);

% for i=1:1:N
%     
%     if i>N
%         t = xscale(X(1+bsize*(i):bsize*i+1,:),W,P,ncomp);
%     end
% end

%Determining the Y predictions for the next block of data in the recursive
%process, i.e. the y values for the second data block is estimated using
%PLS parameters (B and Q) determined from previous data block
for i=1:1:(N-1)

    [t,~] = xscale(X(1+bsize*(i):bsize*(i+1),:),Wcat(:,:,i),Pcat(:,:,i),ncomp,mu);
    
    for j=1:1:size(t,1)
        prep(j,:) = t(j,:).*bcat(i,:).*Qcat(i,:);
    end 
    
    Ylog(1+bsize*(i-1):bsize*i,:) = sum(prep,2);
    Ylog(1+bsize*(i-1):bsize*i,:) = Ylog(1+bsize*(i-1):bsize*i,:)+mean(Y);
    
end

%acc = abs(Ytest-Yt)./abs(Ytest)*100;

%% Determining residuals
% The residuals are determined by feeding the parameter estimates to the
% CSTR model. The residuals are calculated from the difference between the
% model predictions and the actual values. As inputs are discontinuos
% interpolation is used to estimate smooth linear functions between the
% data points. The residuals are then fed to recursive PLS as an additional
% feature as well as to the SVM classifier (possibly).

% set1 = par(1:701)';
% time = (0:0.1:70)';
% %m_F = (F_crt_nrm(2:end)-F_crt_nrm(1:end-1))./(time(2:end)-time(1:end-1));
% %c_F = -1*m_F.*time(2:end)+F_crt_nrm(2:end);
% k_0 = Ylog(1:701,1);
% Ypred(1,:) = [Ca_sens_nrm(1) T_out_nrm(1)];
% 
% for i = 1:1:size(k_0,1)-1
%    
%     Fest = @(t) interp1(time,F_crt_nrm,t);
%     Ti_est = @(t) interp1(time,T_in_nrm,t);
%     
%     model = @(t,y) [(Fest(t)/V)*(Ca_i-y(1))-k_0(i)*exp(-EaR/(y(2)))*y(1);
%         (Fest(t)/(V))*(Ti_est(t)-y(2))+H/(HD)*k_0(i)*exp(-EaR/(y(2)))*y(1)-HA/(HD*V)*(y(2)-interp1(time,Tj_sens_nrm,t))];
%     
%    [ts, ys] = ode45(model,[time(i) time(i+1)],[Ca_sens_nrm(i) T_out_nrm(i)]);
%    
%     Ypred(i+1,:) = ys(end,:);
% end

%% Hybrid model with active frist principle estimates and residual estimates

par_reg = par(1:701,1);
par_est(1:50,1) = 10;
%par_est(1:50,1) = par_reg(1:50,1);

M = floor (size(par_reg,1)/bsize);
time = (0:0.1:70)';
Fest = @(t) interp1(time,F_crt_nrm,t);
Ti_est = @(t) interp1(time,T_in_nrm,t);

% First set of residual estimates 

 for k = 1:1:bsize
        
        model = @(t,y) [(Fest(t)/V)*(Ca_i-y(1))-par_est(k)*exp(-EaR/(y(2)))*y(1);
       (Fest(t)/(V))*(Ti_est(t)-y(2))+H/(HD)*par_est(k)*exp(-EaR/(y(2)))*y(1)-HA/(HD*V)*(y(2)-interp1(time,Tj_sens_nrm,t))];
        
    [tp, yp] = ode45(model,[time(k) time(k+1)],[Ca_sens_nrm(k) T_out_nrm(k)]);
   
    Ypred(k,:) = yp(end,:);    
        
 end

res(1:50,:) = Ypred(1:50,:) - [Ca_sens_nrm(1:50,:) T_out_nrm(1:50,:)];
mu = mean([Dataset(1:50,:) res(1:50,:)]);
% remaining residuals

for i = 1:1:M
    %[~,~,b,~,P,Q,W] = npls([Dataset(1+bsize*(i-1):bsize*i,:) res(1+bsize*(i-1):bsize*i,:)],par_est(1+bsize*(i-1):bsize*i,1), ncomp);
    [~,~,b,~,P,Q,W] = npls([Dataset(1+bsize*(i-1):bsize*i,:) res(1+bsize*(i-1):bsize*i,:)],par_reg(1+bsize*(i-1):bsize*i,1), ncomp);
    
    % this is what the code originally was where the parameter estimates
    % where the PLS model is trained on was obtained from non-lin LSQ
    % estimation, it was changed because it made more sense to only use
    % non-lin LSQ to obtain initial estimates and then train PLS model on
    % subsequent parameter estimates obtained from PLS prediction on
    % previous data block.
    
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
    
    if i < M
    
        strt = 1+bsize*i;
        fin = bsize*(i+1);
        rs = bsize*(i-1)+1;
        rf = i*bsize;
        
    elseif i == M
        
        strt = 1+bsize*i;
        fin = size(par_reg,1);
        rs = bsize*(i-1)+1;
        rf = bsize*(i-1)+1+fin-strt;
    end
    
    [tscale,mu] = xscale([Dataset(strt:fin,:) res(rs:rf,:)],W,P,ncomp,mu);
    	
    par_hold = tscale.*b.*Q;
    
    par_est(strt:fin,:)= sum(par_hold,2);
    par_est(strt:fin,:)= par_est(strt:fin,:)+mean(par(:,1));
    
    
    for k = strt:1:fin
        
        model = @(t,y) [(Fest(t)/V)*(Ca_i-y(1))-par_est(k)*exp(-EaR/(y(2)))*y(1);
       (Fest(t)/(V))*(Ti_est(t)-y(2))+H/(HD)*par_est(k)*exp(-EaR/(y(2)))*y(1)-HA/(HD*V)*(y(2)-interp1(time,Tj_sens_nrm,t))];
        
    [tp, yp] = ode45(model,[time(k-1) time(k)],[Ca_sens_nrm(k-1) T_out_nrm(k-1)]);
   
    Ypred(k,:) = yp(end,:);    
        
    end
    
    res(strt:fin,:) = Ypred(strt:fin,:) - [Ca_sens_nrm(strt:fin,:) T_out_nrm(strt:fin,:)]; 
    
end

P_est_nrm = par_est;

%% k_0 error residuals estimation

par_k0 = par(702:1402,1);
Xk0 = [Tj_sens_fin,F_crt_fin,T_in_fin,Ca_sens_fin,T_out_fin];
Yk0 = [Ca_sens_fin,T_out_fin];
cons = [Ca_i EaR V HA HD H];

[res_k0, Ypred_k0,P_est_k0] = res_est(Xk0,Yk0,par_k0,bsize,M,10,ncomp,cons,time,lambda);

%% Inlet error residuals estimation

par_Ca = par(1403:2103,1);
XCa = [Tj_sens_serr,F_crt_serr,T_in_serr,Ca_sens_serr,T_out_serr];
YCa = [Ca_sens_serr,T_out_serr];


[res_Ca, Ypred_Ca,P_est_Ca] = res_est(XCa,YCa,par_Ca,bsize,M,10,ncomp,cons,time,lambda);

%% Prep and normal SVM classifier training 

t = Ca_in.time; 
input_nrm = [Ca_sens_nrm, F_crt_nrm, T_in_nrm, T_out_nrm, Tj_sens_nrm, P_est_nrm, res];

input_fin = [Ca_sens_fin, F_crt_fin, T_in_fin, T_out_fin, Tj_sens_fin, P_est_k0, res_k0];

input_serr= [Ca_sens_serr, F_crt_serr, T_in_serr, T_out_serr, Tj_sens_serr, P_est_Ca, res_Ca];

nrm_SVM = fitcsvm(input_nrm,ones(size(input_nrm,1),1),'KernelFunction','rbf');

% Preparing labels for SVM models 
Label_norm = strings (701,1);
k0Labels = strings (701,1);
Ca_iLabels = strings (701,1);

Label_norm(1:701) = 'Normal';
%% Identifying input error and labelling - k0 error

[~, scores] = predict(nrm_SVM, input_fin); 

% Setting up error identification
[k, j] = size(input_fin);
ei = 1;
ni = 1;
k0Labels = strings(k,1);

% Identifying input error

for i = 1:1:k
    
    if scores (i) < 0
        errorT (ei,:) = input_fin(i,:);
        te (ei) = t(i);
        ei = ei+1;
        k0Labels(i) = "k0 Error";
    else 
        nonT (ni,:) = input_fin(i,:);
        tn (ni) = t(i);
        ni = ni+1;
        k0Labels(i) = "Normal";
    end
    
end

k0Labels(1:151) = 'Normal';

%% Identifying input error and labelling - Ca_i error

[~, scores] = predict(nrm_SVM, input_serr); 

% Setting up error identification
[k, j] = size(input_fin);
ei = 1;
ni = 1;
Ca_iLabels = strings(k,1);

% Identifying input error

for i = 1:1:k
    
    if scores (i) < 0
        errorT (ei,:) = input_fin(i,:);
        te (ei) = t(i);
        ei = ei+1;
        Ca_iLabels(i) = "Ca_i Error";
    else 
        nonT (ni,:) = input_fin(i,:);
        tn (ni) = t(i);
        ni = ni+1;
        Ca_iLabels(i) = "Normal";
    end
    
end

Ca_iLabels(1:151) = 'Normal';

%% Training two class SVM's - Normal & k0 Error

Data = [input_nrm;input_fin];
Labels = [Label_norm;k0Labels];

Norm_fin_SVM = fitcsvm(Data,Labels,'KernelFunction','rbf');

%% Training two class SVM's - Normal & Ca Error

Data = [input_nrm;input_serr];
Labels = [Label_norm;Ca_iLabels];

Norm_Serr_SVM = fitcsvm(Data,Labels,'KernelFunction','rbf');

%% Ca_i and k0 error SVM
i = find(Ca_iLabels == 'Ca_i Error');
j = find(k0Labels == 'k0 Error');
TrainSet_Serr_fin = [input_serr(i,:); input_fin(j,:)];
Labels = [Ca_iLabels(i,:); k0Labels(j,:)];
Serr_fin_SVM = fitcsvm(TrainSet_Serr_fin,Labels,'KernelFunction','rbf');

%% Multiclass SVM testing - k0 training set

%test_all = [Ca_sens_test, F_crt_test, T_in_test, T_out_test, Tj_sens_test];
test_all = input_fin;
[m, ~] = size(test_all);  
All_Label = strings(m,1);

% Running all SVM classifiers on test set 

[Serr_fin, score_1] = predict(Serr_fin_SVM,test_all);
[Norm_fin, score_2] = predict(Norm_fin_SVM,test_all);
[Norm_Serr, score_3] = predict(Norm_Serr_SVM,test_all);

k = 0;
l = 0;
s = 0;

for i = 1:1:m

% Initialise counters
ns = 0;
nf = 0;
nr = 0;
    
    if strncmp(Serr_fin(i), 'Ca_i Error',5)
        ns = ns+1;
    else
        nf = nf+1;
    end
    
    if strncmp(Norm_fin(i), 'k0 Error',5)
        nf = nf+1;
    else
        nr = nr+1;
    end
    
    if strncmp(Norm_Serr(i), 'Ca_i Error',5)
        ns = ns+1;
    else
        nr = nr+1;
    end
       
    if nf == 2
        All_Label(i) = 'k0 Error';
        k = k+1;
        in_err(k,:) = test_all(i,:);
        t_err(k) = t(i);
    
    elseif ns == 2
        All_Label(i) = 'Ca_i Error';
        l = l+1;
        sens_err(l,:) = test_all(i,:);
        t_sens(l) = t(i);
    
    else
        All_Label(i) = 'Normal';
        s = s+1;
        no_err(s,:) = test_all(i,:);
        t_no(s) = t(i);
    end
    
end

[mis1, sens1, spec1] = ClassMes(All_Label,15,'k0 Error');

%% Multiclass SVM testing - Inlet error training case

%test_all = [Ca_sens_test_2, F_crt_test_2, T_in_test_2, T_out_test_2, Tj_sens_test_2];
test_all = input_serr;
[m, ~] = size(test_all);  
All_Label = strings(m,1);

% Running all SVM classifiers on test set 

[Serr_fin, score_1] = predict(Serr_fin_SVM,test_all);
[Norm_fin, score_2] = predict(Norm_fin_SVM,test_all);
[Norm_Serr, score_3] = predict(Norm_Serr_SVM,test_all);

k = 0;
l = 0;
s = 0;
t_sens = 0;
t_err = 0;
%t_no = 0;
%in_err = zeros(1,6);
for i = 1:1:m
    
% Initialise counters
ns = 0;
nf = 0;
nr = 0;
    
    if strncmp(Serr_fin(i), 'Ca_i Error',5)
        ns = ns+1;
    else
        nf = nf+1;
    end
    
    if strncmp(Norm_fin(i), 'k0 Error',5)
        nf = nf+1;
    else
        nr = nr+1;
    end
    
    if strncmp(Norm_Serr(i), 'Ca_i Error',5)
        ns = ns+1;
    else
        nr = nr+1;
    end
       
    if nf == 2
        All_Label(i) = 'k0 Error';
        k = k+1;
        in_err(k,:) = test_all(i,:);
        t_err(k) = t(i);
    
    elseif ns == 2
        All_Label(i) = 'Ca_i Error';
        l = l+1;
        sens_err(l,:) = test_all(i,:);
        t_sens(l) = t(i);
    
    else
        All_Label(i) = 'Normal';
        s = s+1;
        no_err(s,:) = test_all(i,:);
        t_no(s) = t(i);
    end
    
end

[mis2, sens2, spec2] = ClassMes(All_Label,15,'Ca_i Error');
%% Multiclass SVM testing - k0 test set

testX = [Ca_sens_test, F_crt_test, T_in_test];

testY = [T_out_test, Tj_sens_test];

[par_test] = svm_hybrid_v2(testX,testY,pg,V,Ca_i,HD);

[res_test, Ypred_test,P_est_test] = res_est(testX,testY,par_test(:,1),bsize,M,10,ncomp,cons,time,lambda);

test_all = [testX,testY,P_est_test,res_test];

[m, ~] = size(test_all);  
All_Label = strings(m,1);

% Running all SVM classifiers on test set 

[Serr_fin, score_1] = predict(Serr_fin_SVM,test_all);
[Norm_fin, score_2] = predict(Norm_fin_SVM,test_all);
[Norm_Serr, score_3] = predict(Norm_Serr_SVM,test_all);

k = 0;
l = 0;
s = 0;

for i = 1:1:m

% Initialise counters
ns = 0;
nf = 0;
nr = 0;
    
    if strncmp(Serr_fin(i), 'Ca_i Error',5)
        ns = ns+1;
    else
        nf = nf+1;
    end
    
    if strncmp(Norm_fin(i), 'k0 Error',5)
        nf = nf+1;
    else
        nr = nr+1;
    end
    
    if strncmp(Norm_Serr(i), 'Ca_i Error',5)
        ns = ns+1;
    else
        nr = nr+1;
    end
       
    if nf == 2
        All_Label(i) = 'k0 Error';
        k = k+1;
        in_err(k,:) = test_all(i,:);
        t_err(k) = t(i);
    
    elseif ns == 2
        All_Label(i) = 'Ca_i Error';
        l = l+1;
        sens_err(l,:) = test_all(i,:);
        t_sens(l) = t(i);
    
    else
        All_Label(i) = 'Normal';
        s = s+1;
        no_err(s,:) = test_all(i,:);
        t_no(s) = t(i);
    end
    
end

[mis3, sens3, spec3] = ClassMes(All_Label,30,'k0 Error');

%% Plot k0 for different errors

figure;
hold on

p1 (1) = plot (t,P_est_nrm);
p1 (2) = plot (t,P_est_k0);
p1 (3) = plot (t,P_est_Ca);
%p1 (4) = plot (t,P_test);
legend (p1,'Normal Operation', 'Catalyst deactivation', 'Inlet reagent decline'); 
xlabel ('Time (min)')
ylabel ('Reaction rate constant estimate (min^-1)')
hold off



%% Scales X and determines Xscores for testing/ prediction input set as defined in the PLS algorithm provided by S.Wold
function [Tscale,mu] =  xscale (X,W,P,ncomp,mean_prev)

if size(X,1) == 1
    X0 = X - mean_prev;
    mu = 0;
else
X0 = X - mean(X);
mu = mean(X);
end

for i = 1:1:ncomp
    
    Tscale(:,i) = X0*W(:,i);
    X0 = X0-Tscale(:,i)*P(i);
end

end

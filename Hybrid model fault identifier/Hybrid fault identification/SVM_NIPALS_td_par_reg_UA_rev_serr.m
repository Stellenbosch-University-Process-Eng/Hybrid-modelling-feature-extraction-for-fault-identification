%% NIPALS SVM loop for sensitivity and specificity data

% The purpose of this MATLAB script is to run the NIPALS-SVM model in a
% loop in order to generate sensitivity and specificity curves continously
% at different error magnitudes.

clear
%load 'CSTR_gen_data_UA 20-08'
%load 'par reg UA 25-08'
%load 'CSTR_gen_data_UA'

%load 'par reg UA 05-09'
%load 'CSTR_gen_data - 08-07'


load 'UA_regressed_par_5'
%load 'collection 10-08'

rng(12)

par(1) = EaR;       %Activation Energy
par(2) = H;         %Heat of reaction
par(3) = HA;        %Overall heat transfer coefficient
par(4) = HD;        %Grouped density and heat capacity
par(5) = k_0;       %Pre-exponential factor
%par(6) = p;
par(6) = cp;        %Heat capacity
par(7) = V;         %Volume
par(8) = Ca_i;      %Inlet concentration

r=0;                % setting window size for dynamic PLS

for k = 41          % loop runs for each fault magnitude considered

%% Obtaining Parameter estimates

ncomp = 2; % def 2
bsize = 100;
mu = 0;

% Ca_in_pg = [Ca_in_nrm; Ca_in_fin(:,41);Ca_in_UA(:,k)];
% Tj_sens_pg = [Tj_sens_nrm; Tj_sens_fin(:,41);Tj_sens_UA(:,k)];
% Ca_sens_pg = [Ca_sens_nrm; Ca_sens_fin(:,41);Ca_sens_UA(:,k)];
% F_crt_pg = [F_crt_nrm; F_crt_fin(:,41);F_crt_UA(:,k)];
% T_in_pg = [T_in_nrm; T_in_fin(:,41);T_in_UA(:,k)];
% T_out_pg = [T_out_nrm; T_out_fin(:,41); T_out_UA(:,k)];
% In = [Tj_sens_pg,F_crt_pg, T_in_pg];
% Out = [Ca_sens_pg, T_out_pg];




%% Preparing datasets
for z = 1:10            % for loop runs through each CV case

%Allocating training and test sets for normal operating conditions      
% nrm = [Ca_sens_nrm, F_crt_nrm, T_in_nrm, T_out_nrm, Tj_sens_nrm];
% nrm = nrm(5:14986,:);
% par = pm_nrm(5:14986,:);
% CVdata_n = cvpartition(nrm(:,1),'KFold',10);
% Xtr_n = nrm(CVdata_n.training(z),:);
% Ytr_n = par(CVdata_n.training(z),:);
% Xtst_n = nrm(CVdata_n.test(z),:);
% Ytst_n = par(CVdata_n.test(z),:);

%Allocating training and test sets for catalyst deactivation fault
k_err = [Ca_sens_fin(:,1), F_crt_fin(:,1), T_in_fin(:,1), T_out_fin(:,1), Tj_sens_fin(:,1)];
k_err_all = k_err;
k_err = k_err(r+1:dur-8,:);
k_err_td = t_delay(k_err,k_err_all,r); % generating time delayed data
pk = pm_fin(r+1:dur-8,1:2);
CVdata_k = cvpartition(pk(:,1),'KFold',10);
Xtr_k = k_err_td(CVdata_k.training(z),:);
Ytr_k = pk(CVdata_k.training(z),:);
Xtst_k = k_err_td(CVdata_k.test(z),:);
Ytst_k = pk(CVdata_k.test(z),:);

%Allocating training and test sets for heat transfer coefficient fault
UA_err = [Ca_sens_serr(1:4001,k), F_crt_serr(1:4001,k), T_in_serr(1:4001,k), T_out_serr(1:4001,k), Tj_sens_serr(1:4001,k)];
UA_all = UA_err;
UA_err = UA_err(r+1:dur-8,:);
UA_err_td = t_delay(UA_err,UA_all,r); % generating time delayed data
CVdata_UA = cvpartition(UA_err(:,1),'KFold',10);
pUA = pm_serr_col_2(r+1:dur-8,1:2,41);
Xtr_UA= UA_err(CVdata_UA.training(z),:);
Ytr_UA = pUA(CVdata_UA.training(z),:);
Xtst_UA = UA_err(CVdata_UA.test(z),:);
Ytst_UA = pUA(CVdata_UA.test(z),:);

Data = [k_err;UA_err];
Data_td = [k_err_td;UA_err_td];
pT = [pk;pUA];
CVdata = cvpartition(Data(:,1),'KFold',10);
Xtr = Data_td(CVdata.training(z),:);
Ytr = pT(CVdata.training(z),:);
Xtst = Data_td(CVdata.test(z),:);
Ytst = pT(CVdata.test(z),:);
%% NIPALS and predictions

% Pre-allocating space for matrices storing the estimated parameters
P_est_k0 = zeros(size(Ytr_k));
P_est_UA = zeros(size(Ytr_UA));
P_est_t = zeros(size(Ytst_k));

% Training PLS model
[T,U,B,~,P,Q,W] = npls (Xtr,Ytr, ncomp);

Ym = mean(pT(CVdata.training(z),:));
k_var = std(pT(CVdata.training(z),:)-Ym);

for w =1:1:2

% tscale = xscale(Xtr_n,W,P,ncomp,mu);
% 
% Ytemp = tscale.*B.*Q(w,:);
% P_est_nrm(:,w) = sum(Ytemp,2);
% P_est_nrm(:,w) = P_est_nrm(:,w).*repmat(k_var(w), size(P_est_nrm,1), 1);
% P_est_nrm(:,w) = P_est_nrm(:,w) + Ym(w);

%Predicting CSTR parameters for training set using trained PLS models
tscale = xscale(k_err_td(CVdata_k.training(z),:),W,P,ncomp);

Ytemp = tscale.*B.*Q(w,:);
P_est_k0(:,w) = sum(Ytemp,2);
P_est_k0(:,w) = P_est_k0(:,w).*repmat(k_var(w), size(P_est_k0,1), 1);
P_est_k0(:,w) = P_est_k0(:,w) + Ym(w);

tscale = xscale(UA_err_td(CVdata_UA.training(z),:),W,P,ncomp);

Ytemp = tscale.*B.*Q(w,:);
P_est_UA(:,w) = sum(Ytemp,2);
P_est_UA(:,w) = P_est_UA(:,w).*repmat(k_var(w), size(P_est_UA,1), 1);
P_est_UA(:,w) = P_est_UA(:,w) + Ym(w);

end

%input_nrm = [Xtr_n, P_est_nrm];
input_fin = [Xtr_k, P_est_k0];
input_serr = [Xtr_UA, P_est_UA];

%% Data Labelling
% Preparing labels for SVM models
index = find (CVdata_k.training(z));

%Pre-allocating varaibles to store data labels
Label_norm = strings (size(index,1),1);
Label_k = strings(15001,1);
Label_C = strings(15001,1);

Label_k (1:end) = 'Normal';
Label_C (1:end) = 'Normal';
Label_norm (1:end) = 'Normal';
t = Tj_sensor.time;


% k0 err Labels 
ind = find (NOC_k0(1,:)==0); % identify observations where fault is active
err = t_k0 (ind);% identify time points corresponding to fault observations
for i = 1:length(k_err)
    
    for j = 1:size(ind,2)-1 
        %identifying observation which are fault based on fault condition
        %allocation done during data generation
        if t(i) >= err(j) && t(i) <= err(j+1) && (err(j+1) - err(j)) <= 1.6
            Label_k(i) = 'k0 Error'; 
        end
    end
    
end
Label_k = Label_k(r+1:dur-8);
k0Labels = Label_k(CVdata_k.training(z)); 

% Ca err Labels 
ind = find (NOC_Ca(41,:)==0);
err = t_k0 (ind);
for i = 1:length(UA_err)
    
    for j = 1:size(ind,2)-1
        if t(i) >= err(j) && t(i) <= err(j+1) && (err(j+1) - err(j)) <= 1.6
            Label_C(i) = 'Ca_i Error';
        end
    end
    
end
Label_C = Label_C(r+1:dur-8);
Ca_iLabels = Label_C(CVdata_k.training(z));

%% Training two class SVM's - Normal & k0 Error

Data = [input_fin];
Labels = [k0Labels];

Norm_fin_SVM = fitcsvm(Data,Labels,'KernelFunction','rbf');

%% Training two class SVM's - Normal & Ca Error

Data = [input_serr];
Labels = [Ca_iLabels];

Norm_Serr_SVM = fitcsvm(Data,Labels,'KernelFunction','rbf');

%% Ca_i and k0 error SVM
i = find(Ca_iLabels == 'Ca_i Error');
j = find(k0Labels == 'k0 Error');
TrainSet_Serr_fin = [input_serr(i,:); input_fin(j,:)];
Labels = [Ca_iLabels(i,:); k0Labels(j,:)];
Serr_fin_SVM = fitcsvm(TrainSet_Serr_fin,Labels,'KernelFunction','rbf');

%% Multiclass SVM testing - k0 training set CV test

M = floor(size(Ytst_k,1)/bsize);

tscale = xscale(k_err_td(CVdata_k.test(z),:),W,P,ncomp);

%k_var = std(pT(CVdata.training(z),1)-Ym);
for w = 1:1:2
Ytemp = tscale.*B.*Q(w,:);
P_est_t(:,w) = sum(Ytemp,2);
P_est_t(:,w) = P_est_t(:,w).*repmat(k_var(w), size(P_est_t,1), 1);
P_est_t(:,w) = P_est_t(:,w) + Ym(w);
end

%Pk(:,:,k)=P_est_t;
Pk = P_est_t;
SqErr = (P_est_t-Ytst_k).^2;
PerErr = abs(Ytst_k-P_est_t)./abs(Ytst_k);
MSE(z,:) = mean(SqErr);
MPE(z,:) = mean(PerErr);

test_all = [Xtst_k,P_est_t];
[m, ~] = size(test_all);  
All_Label = strings(m,1);

% Running all SVM classifiers on test set 

[Serr_fin, score_1] = predict(Serr_fin_SVM,test_all);
[Norm_fin, score_2] = predict(Norm_fin_SVM,test_all);
[Norm_Serr, score_3] = predict(Norm_Serr_SVM,test_all);

w = 0;
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
        w = w+1;
        in_err(w,:) = test_all(i,:);
        t_err(w) = t(i);
    
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

[~, sens1(k), spec1(k)] = sens_test(All_Label,Label_k(CVdata_k.test(z)),'k0 Error');
Sens_k (z,k) = sens1(k);

 %% Multiclass SVM testing - Ca_i training set CV test

M = floor(size(Ytst_UA,1)/bsize);

tscale = xscale(UA_err_td(CVdata_UA.test(z),:),W,P,ncomp);

for w = 1:1:2
Ytemp = tscale.*B.*Q(w,:);
P_est_t(:,w) = sum(Ytemp,2);
P_est_t(:,w) = P_est_t(:,w).*repmat(k_var(w), size(P_est_t,1), 1);
P_est_t(:,w) = P_est_t(:,w) + Ym(w);
end

PUA=P_est_t;
SqErr = (P_est_t-Ytst_k).^2;
PerErr = abs(Ytst_k-P_est_t)./abs(Ytst_k);
MSE_UA(z,:) = mean(SqErr);
MPE_UA(z,:) = mean(PerErr);

test_all = [Xtst_UA, P_est_t];
[m, ~] = size(test_all);  
All_Label = strings(m,1);

% Running all SVM classifiers on test set 

[Serr_fin, score_1] = predict(Serr_fin_SVM,test_all);
[Norm_fin, score_2] = predict(Norm_fin_SVM,test_all);
[Norm_Serr, score_3] = predict(Norm_Serr_SVM,test_all);

w = 0;
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
        w = w+1;
        in_err(w,:) = test_all(i,:);
        t_err(w) = t(i);
    
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

[~, sens2(k), spec2(k)] = sens_test(All_Label,Label_C(CVdata_UA.test(z)),'Ca_i Error');
Sens_UA (z,k) = sens2(k);


end
end

%% Scales X and determines Xscores for testing/ prediction input set as defined in the PLS algorithm provided by S.Wold
function [Tscale,mu] =  xscale (X,W,P,ncomp,mean_prev)

if size(X,1) == 1
    X0 = X - mean_prev;
    X0 = X0./ repmat(std(X0), size(X,1), 1);

    mu = 0;
else
X0 = X - mean(X);
X0 = X0./ repmat(std(X0), size(X,1), 1);

mu = mean(X);
end

for z = 1:1:ncomp
    
    Tscale(:,z) = X0*W(:,z);
    X0 = X0-Tscale(:,z)*P(z);
end

end

%% Time delay function
function [cnd_td] = t_delay(input,all_mat,delay)
cnd_td = input;

for i = 1:1:delay

    td = all_mat ((delay+1)-i:3992-i,:);
    cnd_td=[cnd_td,td];
end

end
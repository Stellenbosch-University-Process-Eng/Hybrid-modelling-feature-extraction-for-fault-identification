%% NIPALS SVM loop for sensitivity and specificity data

% The purpose of this MATLAB script is to run the NIPALS-SVM model in a
% loop in order to generate sensitivity and specificity curves continously
% at different error magnitudes.

clear
load 'UA_regressed_par_5'
load 'UA_regressed_par_6'

rng(12)

par(1) = EaR;
par(2) = H;
par(3) = HA;
par(4) = HD;
par(5) = k_0;
%par(6) = p;
par(6) = cp;
par(7) = V;
par(8) = Ca_i;

for k = 1

%% Obtaining Parameter estimates

ncomp = 2; % def 2
bsize = 100;
mu = 0;
r = 8;





%% Preparing datasets
for z = 1:10
    
 nrm = [Ca_sens_nrm, F_crt_nrm, T_in_nrm, T_out_nrm, Tj_sens_nrm];
% nrm_all = nrm;
% nrm = nrm(1:3985,:);
% %nrm_td = t_delay(nrm,nrm_all,8);
% par = pm_nrm(1:3985,:);
% CVdata_n = cvpartition(nrm(:,1),'KFold',10);
% Xtr_n = nrm(CVdata_n.training(z),:);
% Ytr_n = par(CVdata_n.training(z),1:2);
% Xtst_n = nrm(CVdata_n.test(z),:);
% Ytst_n = par(CVdata_n.test(z),1:2);

k_err = [Ca_sens_fin(:,k), F_crt_fin(:,k), T_in_fin(:,k), T_out_fin(:,k), Tj_sens_fin(:,k)];
k_err_all = k_err;
k_err = k_err(r+1:3985,:);
k_err_td = t_delay(k_err,k_err_all,r);
pk = pm_fin(r+1:3985,1:2);
CVdata_k = cvpartition(pk(:,1),'KFold',10);
Xtr_k = k_err(CVdata_k.training(z),:);
Ytr_k = pk(CVdata_k.training(z),1:2);
Xtst_k = k_err(CVdata_k.test(z),:);
Ytst_k = pk(CVdata_k.test(z),1:2);

Ca_err = [Ca_sens_serr(:,1), F_crt_serr(:,1), T_in_serr(:,1), T_out_serr(:,1), Tj_sens_serr(:,1)];
Ca_all = Ca_err;
Ca_err = Ca_err(r+1:3985,:);
Ca_err_td = t_delay(Ca_err,Ca_all,r);
CVdata_Ca = cvpartition(Ca_err(:,1),'KFold',10);
pC = pm_serr(r+1:3985,1:2);
Xtr_Ca = Ca_err(CVdata_Ca.training(z),:);
Ytr_Ca = pC(CVdata_Ca.training(z),1:2);
Xtst_Ca = Ca_err(CVdata_Ca.test(z),:);
Ytst_Ca = pC(CVdata_Ca.test(z),1:2);

Data = [k_err;Ca_err];
Data_td = [k_err_td;Ca_err_td];
pT = [pk;pC];
CVdata = cvpartition(pT(:,1),'KFold',10);
Xtr = Data_td(CVdata.training(z),:);
Ytr = pT(CVdata.training(z),1:2);
Xtst = Data_td(CVdata.test(z),:);
Ytst = pT(CVdata.test(z),1:2);

%% NIPALS and predictions

[T,U,B,~,P,Q,W] = npls (Xtr,Ytr, ncomp);

Ym = mean(Ytr);
k_var = std(Ytr-Ym);

P_est_k0 = zeros(size(Ytr_k));
P_est_Ca = zeros(size(Ytr_Ca));

for w =1:1:2

%Predicting CSTR parameters for training set using trained PLS models
tscale = xscale(k_err_td(CVdata_k.training(z),:),W,P,ncomp);

Ytemp = tscale.*B.*Q(w,:);
P_est_k0(:,w) = sum(Ytemp,2);
P_est_k0(:,w) = P_est_k0(:,w).*repmat(k_var(w), size(P_est_k0,1), 1);
P_est_k0(:,w) = P_est_k0(:,w) + Ym(w);

tscale = xscale(Ca_err_td(CVdata_Ca.training(z),:),W,P,ncomp);

Ytemp = tscale.*B.*Q(w,:);
P_est_Ca(:,w) = sum(Ytemp,2);
P_est_Ca(:,w) = P_est_Ca(:,w).*repmat(k_var(w), size(P_est_Ca,1), 1);
P_est_Ca(:,w) = P_est_Ca(:,w) + Ym(w);

end

%input_nrm = [Xtr_n, P_est_nrm];
input_fin = [Xtr_k];
input_serr = [Xtr_Ca];

% tscale = xscale(k_err_td,W,P,ncomp);
% 
% Ytemp = tscale.*B.*Q;
% P_est_k = sum(Ytemp,2);
% P_est_k = P_est_k.*repmat(k_var, size(P_est_k,1), 1);
% P_est_k = P_est_k + Ym;

%% Data Labelling
% Preparing labels for SVM models
index = find (CVdata_k.training(z));
Label_norm = strings (size(index,1),1);
Label_k = strings(15001,1);
Label_C = strings(15001,1);

Label_k (1:end) = 'Normal';
Label_C (1:end) = 'Normal';
Label_norm (1:end) = 'Normal';
t = Tj_sensor.time;

Label_SVM = fitcsvm(nrm,ones(4001,1));

% k0 err Labels 
ind = find (NOC_k0(k,:)==0);
err = t_k0 (ind);
for i = 1:length(k_err)
    
    for j = 1:size(ind,2)-1
        if t(i) >= err(j) && t(i) <= err(j+1) && (err(j+1) - err(j)) <= 1.6
            Label_k(i) = 'k0 Error';
        end
    end
    
end

% [K,sK] = predict(Label_SVM,k_err);
% 
% for i = 1:length(k_err)
%     
%     if sK(i) < 0
%         Label_k(i) = 'k0 Error';
%     end
% end

Label_k = Label_k(r+1:3985);
k0Labels = Label_k(CVdata_k.training(z)); 

% Ca err Labels 
ind = find (NOC_Ca(k,:)==0);
err = t_k0 (ind);
for i = 1:length(Ca_err)
    
    for j = 1:size(ind,2)-1
        if t(i) >= err(j) && t(i) <= err(j+1) && (err(j+1) - err(j)) <= 1.6
            Label_C(i) = 'Ca_i Error';
        end
    end
    
end

% SVM labelling

% [C,sC] = predict(Label_SVM,Ca_err);
% 
% for i = 1:length(Ca_err)
%     
%     if sC(i) < 0
%         Label_C(i) = 'Ca_i Error';
%     end
% end

Label_C = Label_C(r+1:3985);
Ca_iLabels = Label_C(CVdata_Ca.training(z));


%% Training two class SVM's - Normal & k0 Error

Data = [input_fin];
Labels = [k0Labels];

Norm_fin_SVM = fitcsvm(Data,Labels);

%% Training two class SVM's - Normal & Ca Error

Data = [input_serr];
Labels = [Ca_iLabels];

Norm_Serr_SVM = fitcsvm(Data,Labels);

%% Ca_i and k0 error SVM
i = find(Ca_iLabels == 'Ca_i Error');
j = find(k0Labels == 'k0 Error');
TrainSet_Serr_fin = [input_serr(i,:); input_fin(j,:)];
Labels = [Ca_iLabels(i,:); k0Labels(j,:)];
Serr_fin_SVM = fitcsvm(TrainSet_Serr_fin,Labels);

%% Multiclass SVM testing - k0 training set CV test

M = floor(size(Ytst_k,1)/bsize);

tscale = xscale(k_err_td(CVdata_k.test(z),:),W,P,ncomp);
P_est_t = zeros(size(Ytst_k));

%k_var = std(pT(CVdata.training(z),1)-Ym);

for w = 1:1:2
Ytemp = tscale.*B.*Q(w,:);
P_est_t(:,w) = sum(Ytemp,2);
P_est_t(:,w) = P_est_t(:,w).*repmat(k_var(w), size(P_est_t,1), 1);
P_est_t(:,w) = P_est_t(:,w) + Ym(w);
end

SqErr = (P_est_t-Ytst_k).^2;
PerErr = abs(Ytst_k-P_est_t)./abs(Ytst_k);

test_all = [Xtst_k];
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
Spec_k (z,k) = spec1(k);

 %% Multiclass SVM testing - Ca_i training set CV test

M = floor(size(Ytst_Ca,1)/bsize);

tscale = xscale(Ca_err_td(CVdata_Ca.test(z),:),W,P,ncomp);
P_est_t = zeros(size(Ytst_Ca));

for w = 1:1:2
Ytemp = tscale.*B.*Q(w,:);
P_est_t(:,w) = sum(Ytemp,2);
P_est_t(:,w) = P_est_t(:,w).*repmat(k_var(w), size(P_est_t,1), 1);
P_est_t(:,w) = P_est_t(:,w) + Ym(w);
end


SqErr = (P_est_t-Ytst_k).^2;
PerErr = abs(Ytst_k-P_est_t)./abs(Ytst_k);


test_all = [Xtst_Ca];
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

[~, sens2(k), spec2(k)] = sens_test(All_Label,Label_C(CVdata_Ca.test(z)),'Ca_i Error');
Sens_UA (z,k) = sens2(k);
Spec_UA (z,k) = spec2(k);

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

    td = all_mat ((delay+1)-i:3985-i,:);
    cnd_td=[cnd_td,td];
end

end
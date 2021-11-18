This following folder contains the code used in the hybrid fault identification 
The hybrid fault identification utilising PLS and dPLS hybrid models is implemented in the script _SVM_NIPALS_td_par_reg._ 
The script loads the generated datasets as well the the parameters regressed during the non-linear regression performed by the parameter regression
scripts.
The generated datasets and regressed parameter estimates are used to train PLS models. The window length for the dPLS model is set using the parameter r, where if r=0 standard PLS is performed.
The PLS model training is done using the npls function. The training parameter estimates for the training set are predicted using the the relevant PLS model.
The data is labelled based the fault and NOC durations set during data generation
The training set and training parameter estimates are used to train the binary classifiers. The performance of the fault identification is done by using the multiclass svm fault identification approach where the fault is based on the majority of labels obtained for the observation. 
The sensitivity and specificity are calculated using the _sens_test_ function based on the predicited data labels and the true data labels.

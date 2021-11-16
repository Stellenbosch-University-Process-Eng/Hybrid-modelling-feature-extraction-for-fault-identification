The following folder contains the scripts used during the parameter regression. The script Par_reg_UA is repsonsible for performing the regression. 
The non-linear regression is implemented using a sliding window as is described in the methodolgy chapter of the thesis. The data generated is first loaded 
for use by the script. THe window width is set and the sliding window is implememented using a for loop which ensure that the sliding window moves procedurely through the dataset.
The CSTR model is described using the model equations and the regession is performed for each observation in the window using the lsqnonlin function. The regression minimises 
an error measure by adjusting the unknown model parameters (_k_ and _UA_) with the error measured calulated using LSQ_int function. LSQ_int function is responsible for calculating the error between the true model outputs and those estimated using the parameters adjusted during the regression. The funcSimulate function is responsible for calculating the 
model outputs using the parameters adjusted during the regression. The funcSimulate function is used within the LSQ_int function. The resutls of the non-linear regression are then 
saved to MATLAB data file labelled _UA_regressed_par.mat_.

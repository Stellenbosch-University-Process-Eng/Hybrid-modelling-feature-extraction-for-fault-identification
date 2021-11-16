The scripts included in the following folder are responsible for generating the data used by the fault identifier. 
The CSTR model is implemented in the Simulink model **_CSTR_gen_UA_err.slx_**. The model parameters used by the Simulink 
model are stored in the data file **_CSTR_In.mat_**. The parameters can be loaded when needed for the Simulink model.
The script **_CSTR_Data_Generation_** is used to generate the model datasets which can be used by the fault identification models. 
The cript loads the model parameters and runs the CSTR model to generate the datasets for the normal, catalyst deactivation, inlet 
concentraion and heat transfer fault cases. The faults are implemented using function files where for each function the fault magnitude,
mean duration and duration variance is specified.
The function file which corresponds to each fault is listed as follows:

* k0Data_fxn - Catalyst deactivation fault
* CaData_fxn - Inlet concentration fault
* UAData_fxn - Heat transfer fault

The time points where a fault is active and when the process is under normal operating conditions 
is tracked and the information stored so that each obseravation can be lablled according to the operating conditions at that particular time.  
The generated data is then saved as a MATLAB data file labelled _CSTR_gen_data_ which can then be loaded as needed by the parameter
regression and fault identification scripts.

The variable naming conventions followed by the generated data is described in the file Read_Me_CSTR_Data_gen.txt

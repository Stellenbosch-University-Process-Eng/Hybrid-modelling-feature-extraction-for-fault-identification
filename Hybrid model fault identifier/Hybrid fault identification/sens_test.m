function [mis, sens, spec] = sens_test (PredLabels, TestLabels, name)

% This script is used to calculate the sensitivity and specificity
% measures used to evaluate the performance of the fault identifcation
% models 

[k,~] = size (PredLabels);
FP = 0; %False positive counter
FN = 0; %False negative counter
TP = 0; %True positive counter 
TN = 0; %True negative counter

for i = 1:1:k
    %for loop where the predicted labels are compared to the true labels
    %and the appropriate counter is increased
    if strcmp(PredLabels(i),name) == 0 && strcmp(TestLabels(i),name) == 1
        FN = FN+1;
    end
    
    if strcmp(PredLabels(i),name) == 1 && strcmp(TestLabels(i),name) == 0
        FP = FP+1;
    end

    if strcmp(PredLabels(i),name) == 0 && strcmp(TestLabels(i),name) == 0
        TN = TN+1;
    end
    
    if strcmp(PredLabels(i),name) == 1 && strcmp(TestLabels(i),name) == 1
        TP = TP+1;
    end
    
end

mis = (FP+FN)/(k);
sens = TP/(TP+FN); %sensitivity is calculated
spec = TN/(TN+FP); %specificity is calculated

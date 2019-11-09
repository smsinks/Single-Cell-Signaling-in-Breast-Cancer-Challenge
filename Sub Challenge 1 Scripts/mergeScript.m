% ==================== Merge script ================= 

% This script merges the prediction to the template data so that they are
% in the same order before submitting to synampe for grading

% load the predictions and the template 
predictions = readtable('finalPredictions100_5_2.csv');
template = readtable('subchallenge_1_template_data.csv');
predictions = movevars(predictions,'cell_line','Before','treatment');

% merge the tables 
subC1Predictions = innerjoin( template(:,1:6), predictions, 'Key',...
    {'cell_line','treatment','time','cellID','fileID'} ) ;
subC1Predictions = sortrows(subC1Predictions,'glob_cellID','ascend') ;

% check that the values are alligned
for ii = 1:6
    % for numeric columns
    if isnumeric(template.(ii))
        assert(all(template.(ii) == subC1Predictions.(ii) ) )
    else
        assert(all(strcmp(template.(ii),subC1Predictions.(ii))))
    end   
    fprintf('\n Assert passed for column %d \n',ii)
end

% load the training data and change the prediction that are below the
% minimin in train to the minimun in train 
fprintf('\n Changing the negative values \n')
train = readtable('trainingData_Time_1.csv');
proteins = {'p_ERK','p_PLCg2','p_HER2','p_S6','p_Akt_Ser473_'};

for ii = 1:length(proteins)
    numLow = subC1Predictions.(proteins{ii}) < ...
        min( train.(proteins{ii}) ) ;
    if any(numLow)
        subC1Predictions.(proteins{ii})(numLow) = ...
            min( train.(proteins{ii}) ) ;
    end
end

% replaces the 
subC1Predictions.Properties.VariableNames(7:end) = ...
        ["p.Akt.Ser473.", "p.ERK",  "p.HER2", "p.PLCg2","p.S6"] ;

% save to excel 
fprintf('\n saving to excel \n')
writetable(subC1Predictions,'subC1_XGBoostPred.csv')

fprintf('\n Done \n')
clc; close all; clear

% cd(['/Users/sinkala/Documents/MATLAB/Dream Challenge Breast',...
%     ' Cancer Network'])

fprintf('\n Loading the data \n')
% load the template
templateFile = readtable('subchallenge_2_template_data.csv');
templateFile.treatment = categorical(templateFile.treatment) ;
templateFile.cell_line = categorical( templateFile.cell_line ) ;

% load the data for xgboost 
firstColumnsXGboost = readtable('subC2FirstColumns.csv') ;
ogXGboost = readtable('subC2PredictionsNonMedianXGBoostOG.csv') ;

% clean up the rmse table by moving the protein to the front
ogXGboost = movevars(ogXGboost, {'treatment','time'},'Before',...
    'cleavedCas');

% combine the two table into one 
finalResults = [firstColumnsXGboost(:, 1:5), ogXGboost(:,3:end) ] ;

clear ogXGboost firstColumnsXGboost
%% Clean up the table 
fprintf('\n Cleaining up the prediction data \n')

% remove cell line whose ID is great than that present in template file
finalResults = finalResults(...
    finalResults.cellID <= max(templateFile.cellID), : ) ;
finalResults.fileID = [] ;

% throw in an assertion 
assert(max(finalResults.cellID) == max(templateFile.cellID))

% convert the treatment and cell_line into categorical 
finalResults.cell_line = categorical(finalResults.cell_line);
finalResults.treatment = categorical(finalResults.treatment);

% sort the column of the final results tables
[~, theseArrange] = ismember(templateFile.Properties.VariableNames, ...
    finalResults.Properties.VariableNames ) ;
finalResults = finalResults(:, theseArrange) ;

fprintf('\n removing values that do not exist in the data \n') 
% delete the rows that have 16 as a time but are not for
toGo = finalResults.treatment == 'iEGFR' & ...
    finalResults.cell_line ~= 'MDAMB231' & finalResults.time == 16 ; 
finalResults(toGo, : ) = [] ;

% remove the 13 time point for the  'MCF12A' cell ine because it does not
% exist in the data
toGo = finalResults.time == 13 & finalResults.cell_line == 'MCF12A' & ...
    finalResults.treatment == 'iEGFR' ;
finalResults(toGo, :) = [] ;

% replace time point 17 for time = 16 for the 'MDAMB231' cell line treated
% with iPI3K
toChange = finalResults.treatment == 'iPI3K' & ...
    finalResults.cell_line == 'MDAMB231' & finalResults.time == 17 ;
finalResults.time(toChange) = 16 ;

% find the group mean of the final results because some cell line have more
% than one replication 
finalResults = grpstats(finalResults,...
        {'cell_line','treatment','time','cellID'}, 'mean',...
        'DataVars',finalResults.Properties.VariableNames(5:end) );
    
% remove the row names and the groupcount
finalResults(:, {'GroupCount'}) = [];
finalResults.Row = [] ;
finalResults.Properties.VariableNames(5:end) = extractAfter(...
        finalResults.Properties.VariableNames(5:end),'mean_' ) ;

clear toGo toChange theseArrange testCellLines theseUnique
%% Return only cell lines that are present in the templateFile 

fprintf('\n Laoding data for the prevous predictions \n')
% add the 17 minutes time point from the previous predictions
previousPred = readtable('subC2_predictions_All_NCA_Lasso.csv');
previousPred.cell_line = categorical(previousPred.cell_line) ;
previousPred.treatment = categorical(previousPred.treatment) ;

% merge with the current predictions
results =  [ finalResults; previousPred];

% get the unique rows of the dataset 
[~ , theseUnique ] = unique(results(:,1:4) ,'first') ;
results = results(theseUnique, :) ;

results = innerjoin(templateFile(:,1:4), results);

% assert that the results in the template match those in the predictions
assert(height(results) == height(templateFile) )

results.Properties.VariableNames(4:end) = ...
    strrep( results.Properties.VariableNames(4:end) ,'_','.') ;

% write the table to excel 
writetable(results, 'subC2_XGBoost_plut_Tree_First.csv')

clear theseUnique them


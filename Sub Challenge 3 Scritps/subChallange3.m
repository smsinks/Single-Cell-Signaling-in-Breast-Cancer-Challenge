%% Sub Challange 3 Prediction

% Subchallenge III: Predict how single cells respond to a novel kinase
% inhibitor.

% The goal of model building is often to make predictions for conditions
% that the model was not trained on. To test this capability of the
% methods, we do not provide imTOR perturbation data for any of the cell
% lines, but we ask the participants to predict the effects of inhibition
% of mTOR for all cell lines that have otherwise complete single cell data,
% see Table in Figure 3.

% Participants have to computationally simulate the effect of inhibition
% for the given time points for each reporter. We ask the participants to
% derive a model that can predict the value of all the single cell markers
% in the missing imTOR condition. Please note that in each condition and
% time point we ask for 10’000 representative cells. For the exact time
% points and reporters please consult the “Templates for submitting the
% predictions” section.

%% Load the data
clear; close all force; clc

% cd(['/Users/sinkala/Documents/MATLAB/Dream Challenge Breast',...
%     ' Cancer Network'])

% laod the template and make the variables compatible
fprintf('\n Loading the template data \n')
template = readtable('subchallenge_3_template_data.csv') ;
template = movevars(template,'treatment','Before','cell_line') ;

% specifiy the cell line that are presents & specify the sampling times
templateTimes = unique(template.time);
cellLines = unique( template.cell_line) ;

% specificy the an array in which the best correlations and drugs will be
% put for each time point
bestDrugAll = cell(1,1) ;

for ii = 1:length(cellLines)
    
    % load the data for the particular time point Load the training data
    fprintf("\n Loading the training data for cell line %s\n", ...
        cellLines{ii})
    
    data = readtable(sprintf("%s.csv", cellLines{ii} )) ;
    
    % fill the missing values sometime linear fill does not work. This is
    % because some of the data dont have protein measurements
    try
        fprintf('\n Filling the missing values \n')
        data = fillmissing(data,'linear' ,...
            'DataVariables',data.Properties.VariableNames(6:end)) ;
    catch
        fprintf('\n Turns out the missing data are more complex \n')
        % load the prediction from subchallenge 1 and use the prediction to
        % fill in the missing protein levels
        subC1preds = readtable('subC1_XGBoostPredictions100_0.68.csv');
        
        % get the current cell line from the data
        myCellLinePred = subC1preds(ismember(...
            subC1preds.cell_line, data.cell_line) , :) ;
        
        % check if the cell line has prediction made then just use a fill
        % missing to replace the missing values
        if isempty(myCellLinePred)
            for kk = 5:width(data)
                if iscell(data.(kk))
                    data.(kk) = str2double(data.(kk)) ;
                end
            end
            % now replace the data 
            data = fillmissing(data,'linear' ,...
            'DataVariables',data.Properties.VariableNames(6:end)) ;
        else
            fprintf(['\n It turns out prediction are available for',...
                'this dataset \n'])
            % if the data has prediction then use those to replace that
            % data that are missing
            data = fillmissing(data,'linear' ,...
            'DataVariables',data.Properties.VariableNames(6:end)) ;       
        end   
    end
    
    % convert to categorical
    data.treatment = categorical(data.treatment) ;
    
    for jj = 1:length(templateTimes)
        % count the cell lines that are available for making prediction get
        % the samples for the that have prediction from the data
        curTemplate = template(ismember(...
            template.cell_line,cellLines{ii}),:);
        curTemplate = curTemplate( ...
            curTemplate.time == templateTimes(jj) , : );
        
        % get the data that matched the template
        curData = data(ismember(data.time ,curTemplate.time), :) ;
        curData = curData(ismember(curData.cellID, curTemplate.cellID),:);
        curData = curData(:,ismember(curData.Properties.VariableNames,...
            curTemplate.Properties.VariableNames)) ;
        
        % get the data for PI3K and MAPK inhibitors for each cell line
        PI3K = curData{curData.treatment == 'iPI3K', 5:end} ;
        
        % find the drug that is most correlated with iPI3K response for
        % this particular cell line at this particular time points
        
        % specify other treatment that aside the iPI3K
        otherDrugs = unique(curData.treatment(curData.treatment ~='iPI3K'));
        bestCor = 0;
        for yy = 1:length(otherDrugs)
            otherDrugData = curData{ ...
                curData.treatment == otherDrugs(yy), 5:end};
            % make the two datasets of the same size so they could be
            % compared using corr2
            if size(otherDrugData, 1) > size(PI3K,1)
                otherDrugData = otherDrugData(1:size(PI3K,1), :) ;
            elseif size(otherDrugData, 1) < size(PI3K,1)
                PI3K = PI3K(1:size(otherDrugData, 1), :) ;
            end
                
            % now get the correlation coefficent with other drugs
            corWithiPI3K = corr2(PI3K , otherDrugData);
            if corWithiPI3K > bestCor
                bestCor = corWithiPI3K ; 
                bestDrug =  otherDrugs(yy);
            end         
        end      
        
        % save a copy of the best drug and the correction coefficent
        bestDrugAll(jj*ii,1) = cellLines(ii);
        bestDrugAll(jj*ii,2) = cellstr(bestDrug);
        bestDrugAll(jj*ii,3) = num2cell( templateTimes(jj) );
        bestDrugAll(jj*ii,4) = num2cell( bestCor);
        
        % get the data for PI3K and MAPK inhibitors for each cell line
        PI3K = curData{curData.treatment == 'iPI3K', 5:end} ;
        theOtherDrug = curData{curData.treatment == bestDrug, 5:end} ;
        
        % reduce the data if they are large or increase if they are small
        if size(PI3K,1) > 10000
            PI3K = PI3K(1:10000,:);
        elseif ~isempty(PI3K)
            PI3K = repmat(PI3K,100,1);
            PI3K = PI3K(1:10000,:);
        end
        
        if size(theOtherDrug,1) > 10000
            theOtherDrug = theOtherDrug(1:10000, :) ;
        elseif ~isempty(theOtherDrug)
            theOtherDrug = repmat(theOtherDrug,100,1);
            theOtherDrug = theOtherDrug(1:10000,:);
        end
        
        % calculate the weight mean of the PI3K and mTOR inhibitors
        % sometimes the MEK data is not complate
        if size(PI3K,1) == 10000 && size(theOtherDrug,1) == 10000
            mTOR = PI3K ;
        elseif size(theOtherDrug,1) < 10000
            mTOR = PI3K ;
        elseif size(PI3K,1) < 10000
            mTOR = theOtherDrug ;
        end
        
        % check that mTOR has 10000 rows
        loopResults = [ curTemplate(:,1:4), ...
            array2table(mTOR,'VariableNames',...
            template.Properties.VariableNames(5:end)) ] ;
        
        fprintf(['\n Adding the predictions at time %0.2f for %s to',...
            ' the template\n'],templateTimes(jj), cellLines{ii} )
        % add the predictions to the template
        if jj == 1 && ii == 1
            finalResults = loopResults ;
        else
            finalResults = [finalResults; loopResults ] ;
        end
        
    end
end

% throw in an asserttion
assert(all( strcmp( finalResults.cell_line ,template.cell_line) ) )

finalResults.Properties.VariableNames(5:end) = ...
   strrep(finalResults.Properties.VariableNames(5:end),'_','.') ; 

fprintf('\n Saving the data to excel \n')
% save to excel
writetable(finalResults,'subC3_Predictions.csv')

% create a table showing drugs correclated with Pi3K
drugCorrelations = array2table(bestDrugAll,'VariableNames', ...
    {'cell_line','Treatment','Time','Correlation'});
drugCorrelations(cellfun(@isempty ,drugCorrelations.Treatment), :) = [];
drugCorrelations.Treatment = categorical(drugCorrelations.Treatment);
summary(drugCorrelations.Treatment)

clear ii jj MEK mTOR PI3K samplingTime curCellLine cellLines ...
    curTemplate curData mTOR myCellLinePred subC1preds templateTimes ...  
    kk loopResults bestCor otherDrugs theOtherDrug corWithiPI3K yy ...
    otherDrugData bestDrugAll bestDrug bestCorrAll

fprintf('\n Done \n')

%% Subchallenge 2: 

% process data to be used for making predictions in subchallange 2

% check if the allPhospho data exist if yes then just load the data

% change directory to the folder where the data are
clc; close; clear
% 
% cd(['/Users/sinkala/Documents/MATLAB/Dream Challenge Breast ',...
%     'Cancer Network/subChallenge2 Data'])

% laod the template
templateFile = readtable('subchallenge_2_template_data.csv') ;

% =========== return data for only only time point ============
samplingTime = unique(templateFile.time) ; % 16 Is missing

for kk = 1 :length(samplingTime)
    fprintf('\n Getting data for time point %d \n',kk)
    
    % read in the phosphoprotein data from the csv files in a loop from
    % the folder
    
    Files = dir('*.*');
    for ii = 3:length(Files) % since the actual files start from four
        fprintf('\n Loading data for cellLines %s: Number %d \n',...
            Files(ii).name , ii-2)
        fprintf('\n The current time point is %d : %d mins \n',...
            kk, samplingTime(kk) )
        
        curCellLines = readtable( Files(ii).name ) ;
        
        % change the missing proteins to double
        for jj = 5:width(curCellLines)
            % check if the values are of cell type
            if iscell( curCellLines.(jj) )
                curCellLines.(jj) = ones(height(curCellLines),1)* NaN ;
            end
        end
        
        % get data for the current time point
        curCellLines = curCellLines( ...
            curCellLines.time == samplingTime(kk), :) ;
        
        if ii == 3
            trainingData = curCellLines ;
        else
            trainingData = vertcat(trainingData, curCellLines) ;
        end
    end
    
    % go to the folder where I am saving the processed data
%     cd(['/Users/sinkala/Documents/MATLAB/Dream Challenge Breast',...
%         ' Cancer Network/subC2 Processed Data'])
    
    % fill the missing data 
    fprintf('\n Filling the missing values \n')
    trainingData = fillmissing(trainingData,'linear' ,...
        'DataVariables',trainingData.Properties.VariableNames(6:end)) ;

    % save to excel
    writetable( trainingData,strcat( strcat(...
        'trainingData_C2_Time_',num2str(kk)),'.csv') )

    clear trainData testData
    
    % return back to the folder where all the excel file for the
    % subchallenge are hosted
%    cd(['/Users/sinkala/Documents/MATLAB/Dream Challenge Breast ',...
%     'Cancer Network/subChallenge2 Data'])
end


% get back to the current directory
% cd(['/Users/sinkala/Documents/MATLAB/Dream Challenge Breast ',...
%     'Cancer Network'])

fprintf('\n Done \n')
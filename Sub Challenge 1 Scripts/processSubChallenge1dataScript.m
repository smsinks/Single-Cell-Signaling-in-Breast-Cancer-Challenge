%% Single Cell Signaling in Breast Cancer Challenge

% This script runs data for single cell signalling in breast cancer for the
% dream challenge with the same name

%  =========================== Overview ==================================
% Signaling underlines nearly every cellular event. Individual cells, even
% if genetically identical, respond to perturbation in different ways. This
% underscores the relevance of cellular heterogeneity, in particular in how
% cells respond to drugs. This is of high relevance since the fact that a
% subset of cells do not respond (or only weakly) to drugs can render this
% drug an ineffective treatment. In spite of its relevance to many
% diseases, comprehensive studies on the heterogeneous signaling in single
% cells are still lacking. We have generated the, to our knowledge,
% currently largest single cell signaling dataset on a panel of 67
% well-characterized breast cancer cell lines by mass cytometry (3’015
% conditions, ~80 mio single cells, 38 markers; Bandura et al. 2009;
% Bendall et al., 2011; Bodenmiller et al., 2012; Lun et al., 2017; Lun et
% al., 2019). These cell lines are, among others, also characterized at the
% genomic, transcriptomic, and proteomic level (Marcotte et al., 2016). We
% ask the community to use these measurements to predict the signaling
% responses of the cell lines to drug treatment at the single-cell level.
% Firstly, we ask to predict certain markers in specific single cells.
% Secondly, we ask to predict the time-dependent response of single cells
% upon treatment with drugs. Finally, we aim to predict the time- dependent
% response of cell lines for which only static, unperturbed data is given.
% The methods developed to answer these questions will allow us to better
% understand the determinants that control single-cell signaling, the
% heterogeneity in drug response of cancer cells, and to push the limits of
% signaling modeling.

% =========================== Challenge questions ========================
% Within the scope of the Challenge, we invite participants to provide
% single-cell predictions for four questions with gradually increasing
% difficulty (see also Figure 4): Subchallenge I: the Challenge starts by
% predicting missing markers of kinase activity in single cells.
% Subchallenge I: the Challenge starts by predicting missing markers of
% kinase activity in single cells. Subchallenge II: we ask the participants
% to predict time-dependent single-cell response upon kinase inhibitor
% treatments for a known perturbation. Subchallenge III: the participants
% predict single-cell response to an mTOR-inhibiting perturbation, for
% which no data is available. Subchallenge IV: Finally, we ask the
% participants to predict perturbation response (on population level) for
% cell-lines for which only unperturbed data is available. Here, the
% participants can rely on the unperturbed data of the cell-lines and the
% perturbation data from other cell lines.

% clear variables and change the directory
clc ; clear ; close all force

% change to the directory where the downloaded subchallenge 1 datasets are
% located
% cd(['/Users/sinkala/Documents/MATLAB/Dream Challenge Breast ',...
% 'Cancer Network'])

%% Subchallenge I: Predict missing markers at the single cell level.

% Predicting missing markers is experimentally a very relevant question as
% missing markers might result from experimental errors (for example in
% this dataset p-PLCg2 and p-HER2 are not available for all cell lines) or
% due to limitations in the panel size (such as in fluorescence
% measurements). In a nutshell, we ask participants in subchallenge I. to
% predict the missing values based on other reporters of the measured
% cells. For six cell lines we do not provide the measurements for the
% following markers: p-ERK, p-PLCg2, p-HER2, p-S6, p-AKT_S473 in any
% treatment conditions - as shown in Figure 1.a and 1.b - and we ask the
% participants to predict their levels on single cell level, in specific
% time points after perturbation. Note, that the data on the other cell
% lines are also available and might be used to improve the prediction.

% =========== return data for only only time point ============
samplingTime = [0, 5.5, 7, 9, 13, 17, 23, 30, 40, 60] ;

for kk = 1:length(samplingTime)

    fprintf('\n Getting data for time point %d \n',kk)

    % read in the phosphoprotein data from the csv files in a loop from
    % the folder

    Files = dir('*.*');
    for ii = 4:length(Files) % since the actual files start from four
        fprintf('\n Loading data for cellLines %s: Number %d \n',...
            Files(ii).name , ii)
        fprintf('\n The current time point is %d : %d \n',...
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

        if ii == 4
            curTimeData = curCellLines ;
        else
            curTimeData = vertcat(curTimeData, curCellLines) ;
        end
    end

    % convert the data to double because xgboost does not seem to work
    % and scikitlearn is bullshit
    %         curTimeData.treatment =
    %         double(categorical(curTimeData.treatment));
    %         curTimeData.cell_line =
    %         double(categorical(curTimeData.cell_line));

    % now we have data for one time point now let get data for each
    % phophoproteins Slice up the data into different time points from
    % 0 to the last get the training data for each protein
    proteins = {'p_ERK','p_PLCg2','p_HER2','p_S6','p_Akt_Ser473_'};

%     % go to the folder where I am saving the processed data
%     cd(['/Users/sinkala/Documents/MATLAB/Dream Challenge Breast',...
%         ' Cancer Network/subC2 Processed Data'])

    for yy = 1:length(proteins)

        % print the protein data being processed
        fprintf('\n Getting data for protein %d : %s \n',yy, ...
            proteins{yy} )

        % get some temporary data
        tempTimeData = curTimeData;

        % get the training data for the current proteins
        if ~exist('trainData','var')
            % only save on training set for eachh time point
            trainData = tempTimeData(...
                ~isnan(tempTimeData.(proteins{yy})),:);

            % save to excel
            writetable( trainData,strcat( strcat(...
                'trainingData_Time_' ,num2str(kk)),'.csv') )
        end

        % get the other proteins that are to be removed from the data
        otherProteins = proteins(...
            ~ismember(proteins, proteins(yy))) ;

        % get the test data for the current proteins
        testData = tempTimeData(isnan(tempTimeData.(proteins{yy})),:);

        % remove the other proteins from the data
        testData = removevars(testData, [otherProteins,proteins{yy}]);

        fprintf('\n Time the saving to excel \n')
        % save the data to excel

        writetable( testData, strcat( strcat( strcat( strcat( ...
            'testData_',proteins{yy}) ,'_Time_') ,num2str(kk)),'.csv'))

    end

    clear trainData testData

    % return back to the folder where all the excel file for the
    % subchallenge are hosted
%     cd(['/Users/sinkala/Documents/MATLAB/Dream Challenge Breast',...
%         ' Cancer Network/subChallenge1 cellLine data'])

end


%% get the training and test data

% split the datasets into training the test sets for each protein

% get the training data for each protein
proteins = {'p_ERK','p_PLCg2','p_HER2','p_S6','p_Akt_Ser473_'};
samplingTime = [0, 5.5, 7, 9, 13, 17, 23, 30, 40, 60] ;

for ii = 1:length(proteins)

    % print the protein data being processed
    fprintf('\n Getting data for protein %d : %s \n',ii, ...
        proteins{ii})

    % get the training data for the current proteins
    trainData = curTimeData( ~isnan(curTimeData.(proteins{ii})), :);

    % get the other proteins that are to be removed from the data
    otherProteins = proteins(~ismember(...
        proteins, proteins(ii))) ;

    % remove the other proteins from the data that will not be used for
    % training
    trainData = removevars(trainData , otherProteins );

    % get the test data for the current proteins
    testData = curTimeData( isnan(curTimeData.(proteins{ii})), :) ;

    % remove the other proteins from the data
    testData = removevars(testData, [otherProteins,proteins{ii}] );

    fprintf('\n Time the saving to excel \n')

    for jj = 1:length(samplingTime)

        % get data for the current time point
        curTrain = trainData(trainData.time == samplingTime(jj), :) ;
        curTest = testData(testData.time == samplingTime(jj), :);

        save( strcat( strcat(proteins{ii},'trainTime') ,jj) , ...
            'curTrain','curTest','-v7.3')

    end

    % save the training and test data for that particular protein into
    % a file so that I dont run out of memory and also delete the train
    % data
    %         write( strcat(testProteins{ii}, 'train'), trainData)
    %         write( strcat(testProteins{ii}, 'test'), testData)
    clear testData trainData curTrain curTest

end


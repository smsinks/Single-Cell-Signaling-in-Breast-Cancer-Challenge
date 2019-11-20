%% Subchallenge 2 Prediction Script

% Cell-lines tend to respond to the same perturbations similarly, however
% some cell-lines show distinct patterns. For example, some markers show an
% activated then relaxed response upon EGF perturbation, but in some a
% sustained activity is observed. This subchallenge aims at inferring the
% signaling dynamics in the presence and absence of drugs in different cell
% lines. We ask the participants to predict the distinct response of the
% cell-lines from the unperturbed data. For this, we selected five cell
% lines (see Figure 7) for which we provide measurements in unperturbed
% states only. This data includes proteomics, transcriptomic, genomics and
% base-line single cell (and median) phosphoproteomics measurements.
% Participants can use the perturbation data of other cell lines to build
% their models.

% Here, we are interested in which dataset and which measures in it are
% important in predicting how a cell line signals. In the longer term,
% understanding the genotype-phenotype relationship on the signaling
% network level will be a significant step toward the clinical use of the
% patient’s genomic information to predict drug sensitivity in personalized
% medicine. Please consult the “Templates for submitting the predictions”
% subsection to see the expected format of the submissions for this
% subchallenge.

%% Load the data
clc; close all; clear

% cd(['/Users/sinkala/Documents/MATLAB/Dream Challenge Breast',...
%     ' Cancer Network'])

fprintf('\n Loading the data \n')
% load the template
templateFile = readtable('subchallenge_2_template_data.csv');
templateFile.treatment = categorical(templateFile.treatment) ;

% read the test data
testCellLines = unique(templateFile.cell_line);

% load the test data
xTest = [] ;
for ii = 1:length(testCellLines)
    fprintf('\n Loading the test data for cell line %d : %s\n', ...
        ii, testCellLines{ii})
    % load the test data
    curTest = readtable(sprintf('%s.csv',testCellLines{ii} ) );
    
    % add to the test data but sometimes this fails
    for jj = 6:width(curTest)
        if iscell(curTest.(jj))
            % convert to the cell array to double
            curTest.(jj) = str2double(curTest.(jj)) ;
        end
    end
    fprintf('\n Filling the missing values \n')
    curTest = fillmissing(curTest,'linear' ,...
        'DataVariables',curTest.Properties.VariableNames(6:end)) ;
    
    % remove cell line that have more cellID than those present in the
    % template file and return only the unique rows in the final
    % results table
    % remove cell line that have more cellID than those present in the
    % template file and return only the unique rows in the final
    % results table
    curTest = curTest(curTest.cellID <= max(templateFile.cellID), :);
    curTest.cellID = categorical(curTest.cellID);
    curTest = grpstats(curTest,...
        {'treatment','cell_line','time','cellID'}, 'mean',...
        'DataVars',curTest.Properties.VariableNames(6:end) );
    
    % clean up the varariables in the table
    curTest.GroupCount = [];
    curTest.Row = [] ;
    curTest.Properties.VariableNames(5:end) = extractAfter(...
        curTest.Properties.VariableNames(5:end),'mean_' ) ;
    
    % clean up the data and make sure it only has 10000 cell lines
    curTest.treatment = categorical(curTest.treatment);
    curTest = curTest(curTest.treatment == 'full', :) ;
    
    if height(curTest) < 10000
        % if the samples are less than 5000
        if height(curTest) < 5000
           curTestMore = curTest;
           curTestMore.cellID = categorical( ...
               (height(curTest)+1:1:height(curTest)*2)' );
           % finally add the current test
           curTest = [curTest; curTestMore];
            
        end
        % get the number of cell 
        numOfCell = max(double(curTest.cellID));
        
        % get the values that take us to 10000 
        curTestMore = curTest(1:10000-numOfCell, :);
        curTestMore.cellID = categorical( (numOfCell+1:10000)' );
        
        % add these to the test data
        curTest = [curTest; curTestMore] ;  
        
    elseif height(curTest) > 10000
        curTest = curTest(1:10000, :);
    end
    
    % finally add to the table
    xTest = [xTest ;curTest ]  ;
end

% clean up the data in xTest and only return the full treatment
xTest = movevars(xTest, 'cell_line','Before','treatment');

% as agree remove pHER2 and pPLCg2
xTest = removevars(xTest,{'p_HER2','p_PLCg2'});

% check which cell line have fewer data than the template files and build
% up these in the xTest file

fprintf('\n Loading the predictions median phospho  data \n')
% load time point 2 data that i will put together to create a train sets
phosphoMedian = readtable('median_phospho_data.csv');
for ii = 4:width(phosphoMedian)
    if iscell(phosphoMedian.(ii))
        phosphoMedian.(ii) = str2double( phosphoMedian.(ii) ) ;
    end
end

% =============== Load the Process the rnaSeq Data ================

% get the gene expression data and add the it to the phospho data
essentialGenesDown = readtable('RNA essential Genes.xlsx', 'Sheet',1) ;
essentialGenesUp = readtable('RNA essential Genes.xlsx', 'Sheet',2) ;
essentialGenesUp = essentialGenesUp(:,[1:3]);

% clean up the data
essentialGenesDown.Properties.VariableNames = ...
    essentialGenesDown{2,:}  ;
essentialGenesUp.Properties.VariableNames = ...
    essentialGenesUp{3,:};
essentialGenesDown([1,2],:) = [];
essentialGenesUp([1,2],:) = [];

% get the essential genes
essentialGenes = [essentialGenesDown.Gene; essentialGenesUp.Gene] ;

% load the rnaSeq data
fprintf('\n Loading the rnaSeq data \n')
rnaSeq = readtable('RNAseqExonOnly_Marcotte.csv');

% clean up the rna seq data and only return the essetioal genes
rnaSeq = rnaSeq(:,3:end);
rnaSeq = movevars(rnaSeq,'cell_line','Before',1);
toKeep = [true, ismember(rnaSeq.Properties.VariableNames,essentialGenes)];
rnaSeq = rnaSeq(:,toKeep) ;

% combine the rnaSeq data with phosphoprotein data and also add the rnaSeq
% data to the subchallenge 1 data
data = innerjoin(phosphoMedian, rnaSeq) ;
xTest = innerjoin(xTest, rnaSeq);

% ============== Load the Process the Proteomics Data ==========

% load the protein data and process it by cleaning it
fprintf('\n Loading and processing proteomics data \n')
proteomics = readtable('Proteomics_log2raw.csv');
proteomics.Var1 = regexprep(proteomics.Var1,'\;+\w*','');
for ii = 2:width(proteomics)
    proteomics.(ii) = str2double(proteomics.(ii)) ;
end

% delete the row with more than 95% NaN values obtain expression
% measurements
ProtExpression = proteomics{:,2:end};
samples = proteomics.Properties.VariableNames(2:end) ;
genes = proteomics.Var1; % obtain the genes

% remove nan values
nanIndices = any(isnan(ProtExpression),2);
ProtExpression(nanIndices,:) = [];
genes(nanIndices) = [];
numel(genes)

% Gene profiling experiments typically include genes that exhibit little
% variation in their profile and are generally not of interest. These genes
% are commonly removed from the data.

for times = 1:4
    mask = genevarfilter(ProtExpression);
    
    ProtExpression = ProtExpression(mask,:);
    genes = genes(mask);
    numel(genes)
    
    % filter out genes below a certain fold change threshold
    [~,ProtExpression,genes] = ...
        genelowvalfilter(ProtExpression,genes,'absval',log2(2));
    numel(genes)
    
    % filter genes below a certain percentile: 
    [~,ProtExpression,genes] = ...
        geneentropyfilter(ProtExpression,genes,'prctile',40);
    numel(genes)
end

% put together all the data
processedProteomics = [genes, array2table(ProtExpression,...
    'VariableNames',samples)] ;
processedProteomics = rows2vars(processedProteomics ,...
    'VariableNamesSource','Var1') ;
processedProteomics.Properties.VariableNames(1) = "cell_line" ;

% clean the cell line names
processedProteomics.cell_line = regexprep( regexprep( ...
    processedProteomics.cell_line ,'\_+\w*','') , '\x','');

% get the mean expression levels across replicates of cell lines
processedProteomics.cell_line = categorical(processedProteomics.cell_line);
processedProteomics = grpstats(processedProteomics,'cell_line') ;
processedProteomics( :,{'GroupCount'} ) = [];
processedProteomics.Properties.VariableNames(2:end) = ...
    extractAfter( processedProteomics.Properties.VariableNames(2:end),...
    'mean_') ;
processedProteomics.Row = []; % remove the row naems
processedProteomics.cell_line = cellstr(processedProteomics.cell_line) ;

% ============== Merge the proteomics with the other data ==========
data = innerjoin(data,processedProteomics);
xTest = innerjoin(xTest, processedProteomics);

% convert treatment to a categorical variable and also fill in the missing
% values
data.treatment = categorical(data.treatment);
data = fillmissing(data,'linear',...
    'DataVariables',data.Properties.VariableNames(4:end)) ;

% remove the variable that are not in the template
data = removevars(data,{'p_HER2','p_PLCg2'});

% get the training data
xTrainAll = data(data.treatment== 'full' & ...
    ~ismember(data.cell_line, testCellLines), :) ;

yTrainAll = data(data.treatment~= 'full' & data.time ~= 0, :) ;

% save the training data and ytrain and xTest for use with XGBoost in Julia
writetable(xTest,'xTestSubC2.csv')
writetable(yTrainAll, 'yTrainSubC2.csv');
writetable(xTrainAll,'xTrainSubC2.csv');

% clear some of the VAriables
clear ii jj curTest SC2cellLines curTrain templateSC2 templateSC1 ...
    SC1cellLines rnaSeq proteomics essentialGenesDown essentialGenesUp ...
    essentialGenes toKeep genes nanIndices ProtExpression times ...
    mask processedProteomics phosphoMedian curTest numOfCell curTestMore

%% Start Training Models for Each Cell Line

% save the training data and ytrain and xTest for use with XGBoost in Julia
xTest = readtable('xTestSubC2.csv');
yTrainAll = readtable('yTrainSubC2.csv');
xTrainAll = readtable('xTrainSubC2.csv');

% define the validation cell line names
validationCellLines = {'184B5','BT474','CAL148','CAMA1','UACC893'};

% the cell line that will be used to make predictions here are
SubC2cellLines = unique(templateFile(:,1:2)) ;

% run this in a loop for each sampling time
treatment = categories(templateFile.treatment) ;
samplingTime = unique(templateFile.time);

% create a table to put the rmse of the model prediction to be used for
% diagnosing the model find where stat5 of located the last proteins
pSTAT4pos = find(ismember(xTrainAll.Properties.VariableNames,'p_STAT5')) ;
rmseTable = table('Size',[1, width(xTrainAll(:,4:pSTAT4pos)) ],...
    'VariableType', repmat("double",1,width(xTrainAll(:,4:pSTAT4pos))),...
    'VariableNames', xTrainAll.Properties.VariableNames(4:pSTAT4pos) ) ;
rmseTable = addvars(rmseTable, 0 ,'NewVariableNames',{'time'} ) ;
rmseTable = addvars(rmseTable, {'k'} ,'NewVariableNames',{'treatment'}) ;

% now run the loop
for ii = 1:length(treatment)
    % make predictions for each treatment given in the training data
    
    % train the model for each time point
    for kk = 1:length(samplingTime)
        
        % load the yTrain data for the current time point and add it to
        % the median cell line data that have the yTrain with median
        % values. This is why my RMSE was low but my dream score very high
        % because the distribution of the median data is not the same as
        % that of the current time point
        
        % here their is no EGF treatment at time point 0 and there are not
        % drug treatmetn measurements taken at time point 5.5 therefore, I
        % have to skip these loops if strcmp(treatment(ii), 'EGF') &&
        % samplingTime(kk) == 0
        if ~strcmp(treatment(ii), 'EGF') && ...
                any(samplingTime(kk) == [5.5, 23,30])
            continue
        end
        
        % only iEGFR has a time point of 16 dont know where they got it
        if samplingTime(kk) == 16 && ~strcmp(treatment(ii),'iEGFR')
            % time point 16 does not appear in other treatments
            continue
        end
            
        % make the training set the validation set when sampling time is
        % equal to zero
        if samplingTime(kk) == 0
            curYtrain = xTrainAll ;
            % 16 DOES NOT EXIST IN THE DATA
        elseif samplingTime(kk) == 16 && strcmp(treatment(ii),'iEGFR')
            % get data for a particular drug from the yTrains and add an
            % identifier to differential the PI3K from the train and y
            % values
            curYtrain = yTrainAll(yTrainAll.treatment == treatment(ii) &...
                yTrainAll.time == 17, :) ;
        else
            % get data for a particular drug from the yTrains and add an
            % identifier to differential the PI3K from the train and y
            % values
            curYtrain = yTrainAll(yTrainAll.treatment == treatment(ii) &...
                yTrainAll.time == samplingTime(kk), :) ;
        end
        
        % inner join with the xTrain with and the yTrainDrug data but
        % remove the treatment type and the time from the training data
        % because we already know that these are both full and 0
        
        % first make sure the xTrain is the same size with variable arrange
        % the same way as in cur_yTrain
        loopCellLines = unique(curYtrain.cell_line) ;
        curXtrain = [];
        for pp = 1:length(loopCellLines)
            % get the variable of xTrainAll that are of the current cell
            % lines
            smallTrain = xTrainAll( ismember(xTrainAll.cell_line ,...
                loopCellLines(pp)), :) ;
            % reduce the size of small train to only that contained in the
            % cur_yTrain data
            try
                smallTrain = smallTrain(1: sum(ismember( ...
                    curYtrain.cell_line,loopCellLines(pp)) ), :) ;
            catch
                fprintf('\n smallTrain has fewer variable than in train\n')
                % just delete that cell line from that yData
                curYtrain(ismember( ...
                    curYtrain.cell_line,loopCellLines(pp)) ,:) = [] ;
                continue
            end
            
            % merge the current small with in another table that will keep
            % all the values
            curXtrain = [curXtrain; smallTrain] ;
        end
        
        % throw in an assertion
        assert(all(strcmp(curYtrain.cell_line,curXtrain.cell_line)))
        
        % finally create matched dataset
        curYtrain.Properties.VariableNames(2:end) = ...
            strcat(curYtrain.Properties.VariableNames(2:end),...
            '_curYtrain');
        curTrainData = horzcat(curXtrain(:,[2,4:end]), curYtrain);
        
        % now the data are matched: let separate them into a train and test
        % set
        trainingData = curTrainData(~ismember(curTrainData.cell_line, ...
            validationCellLines) ,:) ;
        validationData = curTrainData(ismember(curTrainData.cell_line, ...
            validationCellLines) , :) ;
        
        % get the location of b_catenin and p_STAT5 the last point
        cateninStartPos = find(contains(...
            trainingData.Properties.VariableNames ,'b_CATENIN')) ;
        pSTAT4pos = find(ismember(trainingData.Properties.VariableNames,...
            'p_STAT5')) ;
        
        % now train the models in a loop for jj =
        % cateninStartPos(2):width(trainingData)
        for jj = cateninStartPos(1):pSTAT4pos
            
            % find the location of the first of the training set
            trainEndPos = find(contains(...
                trainingData.Properties.VariableNames,'cell_line'))-1;
            
            % here its is the yTrain values that will be changes all the
            % time make sure to exclude the variable names for the time in
            % the test set and treatment
            yTrain = trainingData.(jj);
            yVal = validationData.(jj) ;
            xTrain = trainingData{:,cateninStartPos(1): trainEndPos} ;
            xVal = validationData{:,cateninStartPos(2): end} ;
            
            % create a validation set if the xVal is empty
            if isempty(xVal)
                cvp = cvpartition(trainingData.b_CATENIN,'HoldOut',0.1);
                % get the training set and the validation set
                xTrain = trainingData{cvp.training, ...
                    cateninStartPos(1): trainEndPos} ;
                xVal = trainingData{cvp.test, ...
                    cateninStartPos(2): end} ;
                yTrain = trainingData.(jj)(cvp.training);
                yVal = trainingData.(jj)(cvp.test) ;
            end
            
            
            
            % train a tree because it is fast
            fprintf('\n Training an ensemble regression model \n')
            rng(0,'twister'); % For reproducibility
            template = templateTree('MinLeafSize', 8);
            
            % train multiple models and pick the best performing model to
            % make the final train and predictions
            fprintf('\n Selecting the best model \n')
            modelType = ["BoostedTree","BaggedTree","FineTrue",...
                "CourseTree","nca","ncaBoostedTree","ncaBaggedTree",...
                "ncaFineTree","ncaCourseTree", "ncaGPRrelationalQaud",...
                "ncaGPRexponentGaussian","ncaSVMqaudratic", ...
                "ncaSVMcubic","Lasso","lassoBoostedTree",...
                "lassoBaggedTree","lassoFineTree","lassoCourseTree",...
                "lassoGPRrelationalQaud","lassoGPRexponentGaussian",...
                "lassoSVMqaudratic","lassoSVMcubic"] ;
             bestRMSE = ones(length(modelType),1);
            
            % ================== train the models =====================
            for gg = 1:length(bestRMSE)
                fprintf('\n Training model %d',gg)
                if gg == 1
                    % train an ensemble model with LSboost
                    trainedModel = fitrensemble( xTrain, yTrain, ...
                        'Method', 'LSBoost', 'NumLearningCycles', 100, ...
                        'Learners', template, 'LearnRate', 0.01);
                elseif gg == 2
                    % train a bagged ensemble model
                    trainedModel = fitrensemble(xTrain, yTrain, ...
                        'Method', 'Bag', ...
                        'NumLearningCycles', 100, 'Learners', template);
                elseif gg == 3
                    % train an ensemble model with fine tree
                    trainedModel = fitrtree(xTrain, yTrain, ...
                        'MinLeafSize', 4, 'Surrogate', 'off');
                elseif gg == 4
                    % train the course decision true
                    trainedModel  = fitrtree(xTrain, yTrain, ...
                        'MinLeafSize', 36, 'Surrogate', 'off') ;
                elseif gg == 5
                    % train an NCA model ===== select features using NCA
                    % and lasso ======
                    [toGetIdx, bestlambda] = ...
                        ncaFeaturesSelection(xTrain,yTrain) ;
                    nca = fsrnca( xTrain(:,toGetIdx), yTrain,...
                        'FitMethod','exact','Solver','sgd','Lambda', ...
                        bestlambda,'LossFunction',...
                        'epsiloninsensitive','Epsilon',0.8);
                elseif gg == 6
                    % train an ensemble model with LSboost on NCA data
                    trainedModel = fitrensemble(xTrain(:,toGetIdx),...
                        yTrain,'Method', 'LSBoost','NumLearningCycles',...
                        100,'Learners', template, 'LearnRate', 0.01);
                elseif gg == 7
                    % train a bagged ensemble model on NCA data
                    trainedModel = fitrensemble(...
                        xTrain(:,toGetIdx), yTrain, 'Method', 'Bag', ...
                        'NumLearningCycles', 100, 'Learners', template);
                elseif gg == 8
                    % train a  fine tree on NCA data
                    trainedModel = fitrtree(xTrain(:,toGetIdx),...
                        yTrain, 'MinLeafSize', 4, 'Surrogate', 'off');
                elseif gg == 9
                    % train a course tree on NCA data
                    trainedModel  = fitrtree( xTrain(:,toGetIdx), ...
                        yTrain,'MinLeafSize', 36, 'Surrogate', 'off') ;
                elseif gg == 10
                    % train a guassian processes regression model
                    trainedModel =  fitrgp( xTrain(:,toGetIdx), ...
                        yTrain, 'BasisFunction', 'constant', ...
                        'KernelFunction', 'rationalquadratic', ...
                        'Standardize', true);
                elseif gg == 11
                    % train a guassian processes regression model
                    trainedModel =  fitrgp( xTrain(:,toGetIdx), ...
                        yTrain, 'BasisFunction', 'constant', ...
                        'KernelFunction', 'squaredexponential', ...
                        'Standardize', true);
                elseif gg == 12
                    % This code specifies all the model options and trains
                    % the model.
                    responseScale = iqr(yTrain);
                    if ~isfinite(responseScale) || responseScale == 0.0
                        responseScale = 1.0;
                    end
                    boxConstraint = responseScale/1.349;
                    epsilon = responseScale/13.49;
                    % train a support vector machines regression model
                    trainedModel =  fitrsvm( xTrain(:,toGetIdx), ...
                        yTrain, 'KernelFunction', 'polynomial', ...
                        'PolynomialOrder', 2,'KernelScale', 'auto', ...
                        'BoxConstraint', boxConstraint, 'Epsilon', ...
                        epsilon,'Standardize', true);
                elseif gg == 13
                    % This code specifies all the model options and trains
                    % the model.
                    responseScale = iqr(yTrain);
                    if ~isfinite(responseScale) || responseScale == 0.0
                        responseScale = 1.0;
                    end
                    boxConstraint = responseScale/1.349;
                    epsilon = responseScale/13.49;
                    % train a support vector machines regression model
                    trainedModel =  fitrsvm( xTrain(:,toGetIdx), yTrain,...
                        'KernelFunction', 'polynomial', ...
                        'PolynomialOrder', 3,'KernelScale', 'auto', ...
                        'BoxConstraint', boxConstraint, ...
                        'Epsilon', epsilon, 'Standardize', true);
                elseif gg == 14
                    % trained lasso model
                    [toSelectIdxLasso, B, FitInfo] = ...
                        lassoFeatureSelection(xTrain,yTrain) ;
                    % train a lasso model Use the largest Lambda value such
                    % that the mean squared error (MSE) is within one
                    % standard error of the minimum MSE.
                    idxLambda1SE = FitInfo.Index1SE;
                    coef = B(:,idxLambda1SE);
                    coef0 = FitInfo.Intercept(idxLambda1SE);
                elseif gg == 15
                    % train an ensemble model with LSboost on Lasso data
                    trainedModel = fitrensemble( ...
                        xTrain(:,toSelectIdxLasso), yTrain, ...
                        'Method', 'LSBoost', 'NumLearningCycles', 100, ...
                        'Learners', template, 'LearnRate', 0.01);
                elseif gg == 16
                    % train a bagged ensemble model on Lasso data
                    trainedModel = fitrensemble(...
                        xTrain(:,toSelectIdxLasso), yTrain, 'Method',...
                        'Bag','NumLearningCycles',100,'Learners',template);
                elseif gg == 17
                    % train a  fine tree on Lasso data train a  fine tree
                    % on NCA data
                    trainedModel = fitrtree(xTrain(:,toSelectIdxLasso),...
                        yTrain, 'MinLeafSize', 4, 'Surrogate', 'off');
                elseif gg == 18
                    % train a  course tree on Lasso data
                    trainedModel  = fitrtree( ...
                        xTrain(:,toSelectIdxLasso),yTrain, ...
                        'MinLeafSize', 36, 'Surrogate', 'off') ;
                elseif gg == 19
                    % train a guassian processes regression model
                    trainedModel =  fitrgp( xTrain(:,toSelectIdxLasso),...
                        yTrain, 'BasisFunction', 'constant', ...
                        'KernelFunction', 'rationalquadratic', ...
                        'Standardize', true);
                elseif gg == 20
                    % train a guassian processes regression model
                    trainedModel =  fitrgp( xTrain(:,toSelectIdxLasso), ...
                        yTrain, 'BasisFunction', 'constant', ...
                        'KernelFunction', 'squaredexponential', ...
                        'Standardize', true);
                elseif gg == 21
                    % This code specifies all the model options and trains
                    % the model.
                    responseScale = iqr(yTrain);
                    if ~isfinite(responseScale) || responseScale == 0.0
                        responseScale = 1.0;
                    end
                    boxConstraint = responseScale/1.349;
                    epsilon = responseScale/13.49;
                    % train a support vector machines regression model
                    trainedModel =  fitrsvm(xTrain(:,toSelectIdxLasso), ...
                        yTrain, 'KernelFunction', 'polynomial', ...
                        'PolynomialOrder', 2,'KernelScale', 'auto', ...
                        'BoxConstraint', boxConstraint, 'Epsilon', ...
                        epsilon,'Standardize', true);
                elseif gg == 22
                    % This code specifies all the model options and trains
                    % the model.
                    responseScale = iqr(yTrain);
                    if ~isfinite(responseScale) || responseScale == 0.0
                        responseScale = 1.0;
                    end
                    boxConstraint = responseScale/1.349;
                    epsilon = responseScale/13.49;
                    % train a support vector machines regression model
                    trainedModel =  fitrsvm(xTrain(:,toSelectIdxLasso),...
                        yTrain,'KernelFunction', 'polynomial', ...
                        'PolynomialOrder', 3, 'KernelScale', 'auto', ...
                        'BoxConstraint', boxConstraint, ...
                        'Epsilon', epsilon,'Standardize', true);
                end
                
                % make the prediction on the validation set and find the
                % crossvalidation loss
                if gg >= 5 && gg < 14
                    % make prediction using NCA
                    if gg == 5
                        % the actual nca model
                        yValFit = predict(nca,xVal(:,toGetIdx) );
                    else
                        % the other regression models
                        yValFit = predict(trainedModel , xVal(:,toGetIdx));
                    end
                elseif gg >= 14
                    % the lasso model
                    if gg == 14
                        yValFit = xVal*coef + coef0;
                    else
                        % the other regression models
                        yValFit = predict( trainedModel , ...
                            xVal(:,toSelectIdxLasso));
                    end
                else
                    % predictions for other trained machine learning model
                    yValFit = predict(trainedModel , xVal);
                end
                
                % find the validatation rmse error of the model
                predictionError = yVal  - yValFit ;
                rmse = sqrt(mean(predictionError.^2)) ;
                
                fprintf('\n the rmse for model %s is %0.4f \n',...
                    modelType(gg), rmse)
                
                % save the keep the copy of the RMSE
                bestRMSE(gg) = rmse;
            end
            
            % error cheching and finding the saving the variables to a
            % table
            
            if isnan(rmse)
                error('\nthe RMSE is NaN\n')
            end
            
            % get the best model
            bestModel = modelType( find( bestRMSE == min(bestRMSE),1 ));
            
            % print the current best rmse
            fprintf(['\n best validation (%s) RMSE is: %0.4f for %s ',...
                'after treatment with %s, at %0.1f minute\n' ],...
                bestModel, min(bestRMSE), ...
                strrep( trainingData.Properties.VariableNames{jj}, ...
                '_curYtrain',''), string(treatment(ii)) ,samplingTime(kk))
            
            % add the rmse errors to the table and add the treatment to the
            % table and add the treatment to the table
            rmseTable(kk,jj-1) = num2cell(min(bestRMSE)) ;
            rmseTable(kk,end-1) = num2cell(samplingTime(kk) ) ;
            rmseTable(kk,end) = treatment(ii) ;
            
            
            % put all the data together
            xTrain = [xTrain; xVal] ;
            yTrain = [yTrain; yVal] ;
            
            fprintf('\n Training the final model \n')
            rng(0,'twister'); % For reproducibility
            
            % ======= now make prediction using the best model =========
            if find( bestRMSE == min(bestRMSE),1 ) == 1
                % train an ensemble model with LSboost
                trainedModel = fitrensemble( xTrain, yTrain, ...
                    'Method', 'LSBoost', 'NumLearningCycles', 100, ...
                    'Learners', template, 'LearnRate', 0.01);
            elseif find( bestRMSE == min(bestRMSE),1 ) == 2
                % train a bagged ensemble model
                trainedModel = fitrensemble(...
                    xTrain, yTrain, 'Method', 'Bag', ...
                    'NumLearningCycles', 100, 'Learners', template);
            elseif find( bestRMSE == min(bestRMSE),1 ) == 3
                % train an ensemble model with LSboost
                trainedModel = fitrtree(xTrain, yTrain, ...
                    'MinLeafSize', 4, 'Surrogate', 'off');
            elseif find( bestRMSE == min(bestRMSE),1 ) == 4
                % train the course decision true
                trainedModel  = fitrtree(xTrain, yTrain, ...
                    'MinLeafSize', 36, 'Surrogate', 'off') ;
            elseif find( bestRMSE == min(bestRMSE),1 ) == 5
                % train an nca model
                nca = fsrnca( xTrain(:,toGetIdx), yTrain,...
                    'FitMethod','exact','Solver','sgd','Lambda', ...
                    bestlambda,'LossFunction',...
                    'epsiloninsensitive','Epsilon',0.8);
            elseif find( bestRMSE == min(bestRMSE),1 ) == 6
                % train an ensemble model with LSboost on NCA data
                trainedModel = fitrensemble( ...
                    xTrain(:,toGetIdx), yTrain, ...
                    'Method', 'LSBoost', 'NumLearningCycles', 100, ...
                    'Learners', template, 'LearnRate', 0.01);
            elseif find( bestRMSE == min(bestRMSE),1 ) == 7
                % train a bagged ensemble model on NCA data
                trainedModel = fitrensemble(...
                    xTrain(:,toGetIdx), yTrain, 'Method', 'Bag', ...
                    'NumLearningCycles', 100, 'Learners', template);
            elseif find( bestRMSE == min(bestRMSE),1 ) == 8
                % train a  fine tree on NCA data
                trainedModel = fitrtree(xTrain(:,toGetIdx),...
                    yTrain, 'MinLeafSize', 4, 'Surrogate', 'off');
            elseif find( bestRMSE == min(bestRMSE),1 ) == 9
                % train a course tree on NCA data
                trainedModel  = fitrtree( xTrain(:,toGetIdx), ...
                    yTrain,'MinLeafSize', 36, 'Surrogate', 'off') ;
            elseif find( bestRMSE == min(bestRMSE),1 ) == 10
                % train a guassian processes regression model
                trainedModel =  fitrgp( xTrain(:,toGetIdx), ...
                    yTrain, 'BasisFunction', 'constant', ...
                    'KernelFunction', 'rationalquadratic', ...
                    'Standardize', true);
            elseif find( bestRMSE == min(bestRMSE),1 ) == 11
                % train a guassian processes regression model
                trainedModel =  fitrgp( xTrain(:,toGetIdx), ...
                    yTrain, 'BasisFunction', 'constant', ...
                    'KernelFunction', 'squaredexponential', ...
                    'Standardize', true);
            elseif find( bestRMSE == min(bestRMSE),1 ) == 12
                % This code specifies all the model options and trains the
                % model.
                responseScale = iqr(yTrain);
                if ~isfinite(responseScale) || responseScale == 0.0
                    responseScale = 1.0;
                end
                boxConstraint = responseScale/1.349;
                epsilon = responseScale/13.49;
                % train a support vector machines regression model
                trainedModel =  fitrsvm( xTrain(:,toGetIdx), ...
                    yTrain, 'KernelFunction', 'polynomial', ...
                    'PolynomialOrder', 2,'KernelScale', 'auto', ...
                    'BoxConstraint', boxConstraint, 'Epsilon', ...
                    epsilon,'Standardize', true);
            elseif find( bestRMSE == min(bestRMSE),1 ) == 13
                % This code specifies all the model options and trains the
                % model.
                responseScale = iqr(yTrain);
                if ~isfinite(responseScale) || responseScale == 0.0
                    responseScale = 1.0;
                end
                boxConstraint = responseScale/1.349;
                epsilon = responseScale/13.49;
                % train a support vector machines regression model
                trainedModel =  fitrsvm( xTrain(:,toGetIdx), yTrain,...
                    'KernelFunction', 'polynomial', ...
                    'PolynomialOrder', 3,'KernelScale', 'auto', ...
                    'BoxConstraint', boxConstraint, ...
                    'Epsilon', epsilon, 'Standardize', true);
            elseif find( bestRMSE == min(bestRMSE),1 ) == 14
                % trained lasso model
                [~, B, FitInfo] = ...
                    lassoFeatureSelection(xTrain,yTrain) ;
                % train a lasso model Use the largest Lambda value such
                % that the mean squared error (MSE) is within one standard
                % error of the minimum MSE.
                idxLambda1SE = FitInfo.Index1SE;
                coef = B(:,idxLambda1SE);
                coef0 = FitInfo.Intercept(idxLambda1SE);
            elseif find( bestRMSE == min(bestRMSE),1 ) == 15
                % train an ensemble model with LSboost on Lasso data
                trainedModel = fitrensemble( ...
                    xTrain(:,toSelectIdxLasso), yTrain, ...
                    'Method', 'LSBoost', 'NumLearningCycles', 100, ...
                    'Learners', template, 'LearnRate', 0.01);
            elseif find( bestRMSE == min(bestRMSE),1 ) == 16
                % train a bagged ensemble model on Lasso data
                trainedModel = fitrensemble(...
                    xTrain(:,toSelectIdxLasso), yTrain, 'Method',...
                    'Bag','NumLearningCycles',100,'Learners',template);
            elseif find( bestRMSE == min(bestRMSE),1 ) == 17
                % train a  fine tree on Lasso data train a  fine tree on
                % NCA data
                trainedModel = fitrtree(xTrain(:,toSelectIdxLasso),...
                    yTrain, 'MinLeafSize', 4, 'Surrogate', 'off');
            elseif find( bestRMSE == min(bestRMSE),1 ) == 18
                % train a  course tree on Lasso data
                trainedModel  = fitrtree( ...
                    xTrain(:,toSelectIdxLasso),yTrain, ...
                    'MinLeafSize', 36, 'Surrogate', 'off') ;
            elseif find( bestRMSE == min(bestRMSE),1 ) == 19
                % train a guassian processes regression model
                trainedModel =  fitrgp( xTrain(:,toSelectIdxLasso), ...
                    yTrain, 'BasisFunction', 'constant', ...
                    'KernelFunction', 'rationalquadratic', ...
                    'Standardize', true);
            elseif find( bestRMSE == min(bestRMSE),1 ) == 20
                % train a guassian processes regression model
                trainedModel =  fitrgp( xTrain(:,toSelectIdxLasso), ...
                    yTrain, 'BasisFunction', 'constant', ...
                    'KernelFunction', 'squaredexponential', ...
                    'Standardize', true);
            elseif find( bestRMSE == min(bestRMSE),1 ) == 21
                % This code specifies all the model options and trains the
                % model.
                responseScale = iqr(yTrain);
                if ~isfinite(responseScale) || responseScale == 0.0
                    responseScale = 1.0;
                end
                boxConstraint = responseScale/1.349;
                epsilon = responseScale/13.49;
                % train a support vector machines regression model
                trainedModel =  fitrsvm(xTrain(:,toSelectIdxLasso), ...
                    yTrain, 'KernelFunction', 'polynomial', ...
                    'PolynomialOrder', 2,'KernelScale', 'auto', ...
                    'BoxConstraint', boxConstraint, 'Epsilon', ...
                    epsilon,'Standardize', true);
            elseif find( bestRMSE == min(bestRMSE),1 ) == 22
                % This code specifies all the model options and trains the
                % model.
                responseScale = iqr(yTrain);
                if ~isfinite(responseScale) || responseScale == 0.0
                    responseScale = 1.0;
                end
                boxConstraint = responseScale/1.349;
                epsilon = responseScale/13.49;
                % train a support vector machines regression model
                trainedModel =  fitrsvm( xTrain(:,toSelectIdxLasso),...
                    yTrain,'KernelFunction', 'polynomial', ...
                    'PolynomialOrder', 3, 'KernelScale', 'auto', ...
                    'BoxConstraint', boxConstraint, ...
                    'Epsilon', epsilon,'Standardize', true);
                
            end
            
            %  ========= make predictions on the test data =============
             % get the test cell line for each drug that will be used to
            % make predictions
            fprintf('\n Getting the test data \n')
            curCellLines = SubC2cellLines.cell_line( ...
                ismember( SubC2cellLines.treatment, treatment(ii) ) ) ;
            curXtest = xTest(ismember(xTest.cell_line,curCellLines), :) ;
            
            % make this into a matrix
            curXtest = curXtest{:,5:end} ;
            
            if find( bestRMSE == min(bestRMSE),1 ) >= 5 && ...
                    find( bestRMSE == min(bestRMSE),1 ) < 14
                if find( bestRMSE == min(bestRMSE),1 ) == 5
                    % make prediction using NCA
                    yFit = predict(nca, curXtest(:,toGetIdx) );
                else
                    % make predictions for other machine learning models
                    yFit = predict(trainedModel, curXtest(:,toGetIdx)  );
                end
            elseif find( bestRMSE == min(bestRMSE),1 ) == 14
                % predictions for Lass Predict exam scores for the test
                % data. Compare the predicted values to the actual exam
                % grades using a reference line.
                yFit = curXtest*coef + coef0;
            elseif find( bestRMSE == min(bestRMSE),1 ) > 14
                % make prediction of the lasso based models
                yFit = predict(trainedModel, curXtest(:,toSelectIdxLasso));
            else
                % make predictions for model
                yFit = predict(trainedModel, curXtest  );
            end
            
          % make predictions and add the predictions to the table
            if jj == cateninStartPos(1)
                % get the first there columns of the test data and change
                % the times point and treatement of that of the current
                % loop
                trainResults = xTest(ismember(...
                    xTest.cell_line,curCellLines),1:4) ;
                trainResults.treatment = repmat(treatment(ii), ...
                    size(curXtest,1),1) ;
                trainResults.time = repmat(samplingTime(kk), ...
                    size(curXtest,1),1) ;
                
                % add the predictions to the table
                trainResults = addvars( trainResults, ...
                    yFit, 'NewVariableNames', ...
                    strrep( trainingData.Properties.VariableNames(jj),...
                    '_curYtrain','') );
                
            else
                trainResults = addvars( trainResults, ...
                    yFit, 'NewVariableNames', ...
                    strrep( trainingData.Properties.VariableNames(jj),...
                    '_curYtrain','') );
            end
        end
        
        % Combine the tables for the diffferent time points
        if ~exist('finalResults','var')
            finalResults = trainResults ;
        else
            finalResults = [ finalResults ; trainResults] ;
        end
        
        if ii == 1 && kk == 1
            finalRMSE = rmseTable ;
        else
            finalRMSE = [ finalRMSE ; rmseTable ];
        end     
    end
    
end

% return back to the folder where the results will be saved
% cd(['/Users/sinkala/Documents/MATLAB/Dream Challenge Breast',...
%     ' Cancer Network'])

rmseTable = finalRMSE ;

% clean up the rmse table by moving the protein to the front
rmseTable = movevars(rmseTable, {'treatment','time'},'Before',...
    'b_CATENIN');

% tree LSBoost boasted Dream Score RMSE = 0.3063 
% Your submission has been scored! Score: 84.81431

% rmse for all model + NCA + lasso: Mean RMSE = 0.1857
% Your submission has been scored! Score: 75.868852

disp(rmseTable)
fprintf('\n Mean RMSE = %0.4f \n', mean(mean(rmseTable{:,3:end})) )

% set up the data for comparison
templateFile.cell_line = categorical(templateFile.cell_line);
finalResults.cell_line = categorical(finalResults.cell_line);
finalResults.treatment = categorical(finalResults.treatment);

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

% assert that the results in the template match those in the predictions
assert(height(finalResults) == height(templateFile) )

finalResults = sortrows(finalResults,'cell_line','ascend');
finalResults.Properties.VariableNames(3:end) =...
    strrep(finalResults.Properties.VariableNames(3:end),'_','.') ;

% save to excel
writetable(finalResults,'subC2_predictions.csv')

clear ii jj cur_yTrain curTrainData testCellLines treatment ...
    validationData validationCellLines predictionError regressionGP...
    rmse  xTrain xVal yFit yTrain yVal yValFit kk samplingTime ...
    curXTrain loopCellLines pp trainEnd varsStart statsResults ans...
    cateninStartPos curXtrain curYtrain phosphoMedian pSTAT4pos ...
    samples trainEndPos  trainingData smallTrain trainResults ...
    curCellLines toChange toGo curXtest curTestMore SubC2cellLines ...
    numOfCell template regressionEnsemble finalRMSE 



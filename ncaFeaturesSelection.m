%% Tune Regularization Parameter in NCA for Regression
% Load the sample data.

% This function is taken from the MATLAB code for Neighborhood Component
% Analysis (NCA) feature selection for regression

function [toGetIdxNCA, bestlambda] = ncaFeaturesSelection(xData,yData)

% INPUT:
% geneticData: a table containing genetic data with the drug on the end
% column.
% numberOfFeatures is the number of features to select from the NCA values

% OUTPUT:
% genetic_train: training data with specific features selected

rng(1); % For reproducibility
cvp = cvpartition(yData,'holdout',0.2)  ;% 56
Xtrain = xData(cvp.training,:);
ytrain = yData(cvp.training,:);
Xtest  = xData(cvp.test,:);
ytest  = yData(cvp.test,:);

% The |robotarm| (pumadyn32nm) dataset is created using a robot arm simulator
% with 7168 training observations and 1024 test observations with 32 features [1][2].
% This is a preprocessed version of the original data set. The data are
% preprocessed by subtracting off a linear regression fit, followed by
% normalization of all features to unit variance.
% Perform neighborhood component analysis (NCA) feature selection for
% regression with the default $\lambda$ (regularization parameter) value.

nca = fsrnca(Xtrain,ytrain,'FitMethod','exact','Solver','lbfgs');

% More than half of the feature weights are nonzero. Compute the loss using the test
% set as a measure of the performance by using the selected features.
fprintf('\n Running NCA for regression \n')
L = loss(nca,Xtest,ytest);

% Try improving the performance. Tune the regularization
% parameter $\lambda$ for feature selection using five-fold
% cross-validation. Tuning $\lambda$ means finding the $\lambda$ value that
% produces the minimum regression loss. To tune $\lambda$ using cross-validation:
%
% 1. Partition the data into five folds. For each fold, |cvpartition|
% assigns 4/5th of the data as a training set, and 1/5th of the data as
% a test set.
rng(1) % For reproducibility
n = length(ytrain);
cvp = cvpartition(length(ytrain),'kfold',5);
numvalidsets = cvp.NumTestSets;

% Assign the $\lambda$ values for the search. Create an array to store the
% loss values.
lambdavals = linspace(0,50,20)*std(ytrain)/n;
lossvals = zeros(length(lambdavals),numvalidsets);

% 2. Train the NCA model for each
% $\lambda$ value, using the training set in each fold.
%
% 3. Compute the regression loss for the corresponding test
% set in the fold using the NCA model. Record the loss value.

% ============ Use built-in robust loss function ================
% -insensitive loss function is more robust to outliers than mean squared
% error.
try
    % 4. Repeat this for each $\lambda$ value and each fold.
    for i = 1:length(lambdavals)
        for k = 1:numvalidsets
            X = Xtrain(cvp.training(k),:);
            y = ytrain(cvp.training(k),:);
            Xvalid = Xtrain(cvp.test(k),:);
            yvalid = ytrain(cvp.test(k),:);
            
            nca = fsrnca(X,y,'FitMethod','exact', 'Solver',...
                'minibatch-lbfgs','Lambda',lambdavals(i), ...
                'GradientTolerance',1e-6,'IterationLimit',60,...
                'LossFunction','epsiloninsensitive',...
                'Standardize',true);
            
            lossvals(i,k) = loss(nca,Xvalid,yvalid,'LossFunction','mse');
        end
    end
catch
    % 4. Repeat this for each $\lambda$ value and each fold.
    for i = 1:length(lambdavals)
        for k = 1:numvalidsets
            X = Xtrain(cvp.training(k),:);
            y = ytrain(cvp.training(k),:);
            Xvalid = Xtrain(cvp.test(k),:);
            yvalid = ytrain(cvp.test(k),:);
            
            nca = fsrnca(X,y,'FitMethod','exact', 'Solver',...
                'minibatch-lbfgs','Lambda',lambdavals(i), ...
                'GradientTolerance',1e-6,'IterationLimit',60,...
                'LossFunction','epsiloninsensitive',...
                'Standardize',true);
            
            lossvals(i,k) = loss(nca,Xvalid,yvalid,'LossFunction','mse');
        end
    end
end

% Compute the average loss obtained from the folds for each $\lambda$ value.
meanloss = mean(lossvals,2);

% Find the $\lambda$ value that gives the minimum loss value.
[~,idx] = min(meanloss) ;

bestlambda = lambdavals(idx);
bestloss = meanloss(idx);

% Fit neighborhood component analysis model using -insensitive loss
% function and best lambda value.

nca = fsrnca(X,y,'FitMethod','exact','Solver','sgd',...
    'Lambda',bestlambda,'LossFunction','epsiloninsensitive',...
    'Standardize',true);

% Plot the selected features.
% figure(101)
% plot(nca.FeatureWeights,'ro')
% xlabel('Feature Index')
% ylabel('Feature Weight')
% grid on

% Most of the feature weights are zero. |fsrnca| identifies the four most
% relevant features.

% Compute the loss for the test set.

L2 = loss(nca,Xtest,ytest);


%% Get the Feature than have weights above a threshold

% get the index of the features that have passed the selection step and the
% varaible names and make a new table
all_features = nca.FeatureWeights';
% to_get_idx = all_features > 1e-8 ;

% get the features top 30 features
toGetIdxNCA = find( maxk(nca.FeatureWeights, 35) , 35) ;



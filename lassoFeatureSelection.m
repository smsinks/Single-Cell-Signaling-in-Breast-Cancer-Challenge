%This function select feature for machine learning using lasso elasltic net

function [toSelectIdxLasso,B,FitInfo]=lassoFeatureSelection(XTrain,yTrain)

% The alpha value changes lasso to elastic net which can select more
% feature than the number of samples.
alphas = [0.1:0.1:0.8] ;

% Option to cross-validate in parallel and specify random streams
opts = statset('UseParallel',true, 'UseSubstreams',true);

% UseSubstreams — Set to true to compute in parallel in a reproducible
% fashion. For reproducibility, set Streams to a type allowing substreams:
% 'mlfg6331_64' or 'mrg32k3a'. The default is false.

% Create a single stream and designate it as the current global stream:
s = RandStream('mlfg6331_64','Seed',1);
RandStream.setGlobalStream(s);

% Find the coefficients of a regularized linear regression model using
% 10-fold cross-validation and the elastic net method with Alpha = 0.1 to 0.8.
% Use the largest Lambda value such that the mean squared error (MSE) is
% within one standard error of the minimum MSE.
for kk = 1:length(alphas)
%     fprintf('\n Extracting features for alpha = %0.1f\n',alphas(1,kk) )
    cur_alpha = alphas(1,kk) ;
    [B,FitInfo] = lasso(XTrain,yTrain,'Alpha',cur_alpha,'CV',10 ,...
        'Options',opts);
    MSE = FitInfo.MSE(FitInfo.Index1SE) ;
    if kk == 1
        idxLambda1SE = FitInfo.Index1SE; 
        mseBestSoFar = MSE;
        Best_beta = B ;
    else
        if mseBestSoFar > MSE 
            mseBestSoFar = MSE ;
            idxLambda1SE = FitInfo.Index1SE ; 
            Best_beta = B ;
        end
    end
end

% Find the nonzero model coefficients corresponding to the two identified
% points. Sometimes, there are no coefficents that are nonezero; in such
% cases the original input table is returned.
fprintf('\n Elastic Net Best Lambda1SE is at %d \n', idxLambda1SE);
toSelectIdxLasso =  Best_beta(:,idxLambda1SE)~=0;

end

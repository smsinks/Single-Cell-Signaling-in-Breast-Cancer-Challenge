# This code processed data for the dream challange using cloud computing
# available at Rescale

# # Installing the requared packages: Activate this line to install the packanges
# using Pkg
# Pkg.add("DataFrames")
# Pkg.add("Statistics")
# Pkg.add("Printf")
# Pkg.add("CSV")
# Pkg.add("DataFrames")
# Pkg.add("Plots")
# Pkg.add("StatsPlots")
# Pkg.add("XGBoost")

# # Load the requared packages
println("\nLoading packages\n")
using XGBoost:predict
using XGBoost
using DataFrames
using Plots, StatsPlots
using CSV, DataFrames
using Statistics
using Printf
using MLDataUtils

# Matlab-like plotting. Installed via Pkg.add("Winston")
# using Winston

# change the directory to were the pre-processed subchallenge 1 are
try
    cd("/Users/sinkala/Documents/MATLAB/Dream Challenge Breast Cancer Network/subC4 Processed Data")
catch
    println("\nThe folder only changes on computer \n")
end

 # load the machine learning data
 # readtable is deprecated, use CSV.read from the CSV package instead
println("\nLoading files \n")

# specify the sampling time
samplingTime = reshape(1:12, 1, 12)
proteins = ["p_Akt_Ser473_", "p_ERK", "p_HER2", "PLCg2", "p_S6"]
global rmseTable = ones(12,6) ;

# make the predictions in a loop
for ii = 1:length(samplingTime)
    # for the sampleing times
    for jj = 1:length(proteins)

        # Load the training data
        println("\n Loading the training data \n")
        toLoadTrain = @sprintf("trainingData_Time_%d.csv",ii)
        train = CSV.read(toLoadTrain)

        # delete the cell line name from the data
        train = deletecols!(train,[:cell_line ])

        # INSERT ANOTHER LOOP HERE FOR EACH CELL LINE

        # Load the test dataset
        # The test data will remain the same it the model that will be changing
        println("\n Getting test set \n")
        toLoadTest = @sprintf("testData_Time_%d.csv",ii)
        xTest = CSV.read(toLoadTest)

        # save a copy of the original dataset. this will make the data smaller
        toLoadTestOG = @sprintf("OGtestData_Time_%d.csv",ii)
        testOG = CSV.read(toLoadTestOG)
        testOG = testOG[:,1:5]

        # getindex(df::DataFrame, col_ind::ColumnIndex)` is deprecated, use `
        # df[!, col_ind] = v
        # delete the cell line names from the data
        xTest = deletecols!(xTest,[:cell_line ])

        # convert to array
        xTest = convert(Array, xTest)

        # Get the training and test data for the current protein
        curProt = proteins[jj]
        print("\n Running model for time: $ii, protein: $curProt \n" )

        # Specify the proteins using the symbol notation used in julia. I dont
        # understand why this is the case
        protSymbol = [:p_Akt_Ser473_,:p_ERK,:p_HER2, :p_PLCg2,:p_S6]
        curProtSymbol = protSymbol[jj]

        # Get the training data first by removing the column for the y test
        # column and then remove all other proteins from the table that are also
        # to be predicted later
        println("\n Getting trainnig and test set \n")
        y = convert(Array, train[:,curProtSymbol] )
        X = deletecols!(train,[:p_Akt_Ser473_,:p_ERK,:p_HER2, :p_PLCg2,
                :p_S6] )
        X = convert(Array, X)

        # # romove the rows with nan values from the training set
        # nanRows = any(isnan(X), 2)
        # X = X[!vec(nanRows), :]
        # y = y[!vec(nanRows), :]

        # split the data into a training and test set
        # shuffle the data so its not in order when we split it up
        Xs, Ys = shuffleobs((transpose(X), y))

        #now split the data into training sets and validation sets
        (xTrain, yTrain), (xTrainTest, yTrainTest) =splitobs((Xs, Ys); at =0.10)

        # need to convert the split data back into arrays
        xTrain = Array(transpose(xTrain))
        yTrain = Array(yTrain)
        xTrainTest = Array(transpose(xTrainTest))
        yTrainTest = Array(yTrainTest)

        # PROVIDE SETTING FOR THE XGBOOST algorithm
        # change booster to gblinear, so that we are fitting a linear model
        # alpha is the L1 regularizer
        # lambda is the L2 regularizer
        # you can also set lambda_bias which is L2 regularizer on the bias term
        # gblinear: linear models gbtree: tree-based models "eta" => 1,

        param = [
            "booster" => "gbtree",
            "silent" => 0,
            "objective" => "reg:linear",
            "alpha" => 0.001,
            "lambda" => 2,
            "eval_metric" => "rmse",
            "max_depth" => 5
        ]

        # normally, you do not need to set eta (step_size)
        # XGBoost uses a parallel coordinate descent algorithm (shotgun),
        # there could be affection on convergence with parallelization on certain
        # cases,setting eta to be smaller value, e.g 0.5 can make the
        # optimization more stable

        # Fit Model
        println("\nTraining Model\n")
        num_round = 100
        bst = xgboost(xTrain, label = yTrain, num_round, param = param)

        # Make prediction on the largeTest test and add the predictions to the table
        yFit = predict(bst, xTrainTest)
        predictionError = yTrainTest  - yFit
        rmse = sqrt(mean(predictionError.^2))
        rmseTable[ii,jj] = rmse
        rmseTable[ii,6] = ii

        print("\n the validation RMSE prediction is $rmse\n")

        # Predict
        println("\n Making predictions \n")
        global pred = predict(bst, xTest)

        if isnan(sum(pred))
            error("\nthe prediction is nan \n")
        end

        # add the predictions to the table and save to excel by load the data
        # which has all the names of the data
        if jj == 1
            global testOG[!,curProtSymbol] = pred
            global curPredictions = testOG
        else
            global curPredictions[!, curProtSymbol] = pred
        end
    end
    # Combine the tables for the diffferent time points
    if ii == 1
        global finalResults = curPredictions
    else
        global finalResults = [ finalResults ; curPredictions]
    end
end

# Convert the array to a table
RMSE = mean(rmseTable[:,1:5])
rmseTable = convert( DataFrame, rmseTable)
colNamesString = ["p_Akt_Ser473_","p_ERK","p_HER2","p_PLCg2","p_S6","Time"]
names!(rmseTable, Symbol.(colNamesString) )

# 100 give RMSE of 0.6216295405183834
# 100 round, 5 depth, lambda 2 , gives a 0.61528492196774
# write to a csv file
print("\n Saving final predictions to excel \n" )
CSV.write("finalPredictions.csv", finalResults)

println("\nDone Processing the data\n")

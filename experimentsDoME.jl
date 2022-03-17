
using FileIO
using JLD2

include("usefulFunctions.jl");
include("DoME.jl");

import Base.OneTo

datasetType=Float64;

Strategies = [StrategyExhaustive, StrategyExhaustiveWithConstantOptimization, StrategySelectiveWithConstantOptimization, StrategySelective];
strategy = 1;



configurations = Array{Int,3}(undef,length(MinimumReductionsMSE),length(MaxNumNodes),numFolds);
n=1; for i1=1:length(MinimumReductionsMSE), i2=1:length(MaxNumNodes), i3=1:numFolds global n; configurations[i1,i2,i3] = n; n += 1; end;
configuration(number::Int) = ( index=findfirst(configurations.==number); (index==nothing)&&(exit()); (MinimumReductionsMSE[index[1]],MaxNumNodes[index[2]],index[3]) );
configurationIndices(minimumReductionMSE, maxNumNodes, numFold) = (findfirst(MinimumReductionsMSE.==minimumReductionMSE), findfirst(MaxNumNodes.==maxNumNodes), numFold);

SourceCodeCompiled = false;

println("Number of executions: ", n-1);


arguments = parseARGS();
if isempty(arguments)
    datasetName = datasetNames()[1];
    configurationsRun = 1:(length(Strategies)*length(MinimumReductionsMSE)*length(MaxNumNodes)*numFolds);
else
    @assert(length(arguments)>=2);
    datasetName = arguments[1];
    strategy = arguments[2];
    if length(arguments)==3
        configurationsRun = arguments[3];
    else
        configurationsRun = 1:(length(Strategies)*length(MinimumReductionsMSE)*length(MaxNumNodes)*numFolds);
    end;
end;



# If the folder "results" does no exist, create it
(!isdir("results")) && mkpath("results");

function loadResults(datasetName, strategy, minimumReductionMSE, maxNumNodes, numFold)
    fileName = resultsFileName(datasetName, strategy);
    results = load(fileName, "results");
    (i1, i2, i3) = configurationIndices(minimumReductionMSE, maxNumNodes, numFold);

    if i1>size(results,1) || i2>size(results,2) || i3>size(results,3)
        return [];
    else
        executionResults = results[i1, i2, i3, :];
        return (any(isnan.(executionResults[[1,3,4,5]])) || any(executionResults[3:end].==0)) ? [] : executionResults;
    end;
end;

function saveResults(datasetName, strategy, minimumReductionMSE, maxNumNodes, numFold, executionResults)
    fileName = resultsFileName(datasetName, strategy);
    results = load(fileName, "results");
    (i1, i2, i3) = configurationIndices(minimumReductionMSE, maxNumNodes, numFold);
    if i1>size(results,1) || i2>size(results,2) || i3>size(results,3)
        newResults = Array{Float64,4}(undef, max(i1,size(results,1)), max(i2,size(results,2)), max(i3,size(results,3)), 5);
        newResults[:] .= NaN;
        if !isempty(results)
            newResults[1:size(results,1), 1:size(results,2), 1:size(results,3), :] .= results[:,:,:,:];
        end;
        results = newResults;
    end;
    results[i1, i2, i3, :] .= executionResults;
    save(fileName, "results", results);
end;


# The function to test
function functionTest(inputs, targets, indicesTest, minimumReductionMSE, maxNumNodes, strategy)
    GC.gc();
    @timed dome(inputs, targets;
        testIndices = indicesTest,
        minimumReductionMSE = minimumReductionMSE,
        maximumHeight = Inf ,
        maximumNodes = maxNumNodes ,
        strategy = Strategies[strategy] ,
        showText = false
    );
end;



if (!isfile(resultsFileName(datasetName, strategy)))
    save(resultsFileName(datasetName, strategy), "results", Array{Float64,4}(undef,0,0,0,0));
end;

(inputs, targets) = loadDataset(datasetName; datasetType=datasetType)

# Create the same cross-validation indices used in the rest of experiments
indicesKFold = crossvalidationIndices(length(targets), numFolds);

println("I'm going to run the following configurations:")
for currentConfiguration in configurationsRun

    (minimumReductionMSE, maxNumNodes, numFold) = configuration(currentConfiguration)
    println("   ", currentConfiguration, ": k-fold with strategy ", strategy, ", minimum reduction in MSE: ", minimumReductionMSE, ", maximum number of nodes: ", maxNumNodes, ", fold ", numFold);

    if isempty(loadResults(datasetName, strategy, minimumReductionMSE, maxNumNodes, numFold))

            global SourceCodeCompiled;
            if (!SourceCodeCompiled)
                # Compile with a "fast" execution
                functionTest(inputs, targets, collect(1:Int(trunc(length(targets)/10))), 1e-4, 50, strategy);
                SourceCodeCompiled = true;
            end;

            ((trainingMSE, validationMSE, testMSE, bestTree), elapsedTime) = functionTest(inputs, targets, findall(indicesKFold.==numFold), minimumReductionMSE, maxNumNodes, strategy);

            executionResults = [trainingMSE, testMSE, elapsedTime, Float64(height(bestTree)), Float64(numNodes(bestTree))];
            saveResults(datasetName, strategy, minimumReductionMSE, maxNumNodes, numFold, executionResults);

            println("      Finished training fold ", numFold, "/", numFolds, " MSE training/test: ", trainingMSE, "/", testMSE, " in ", elapsedTime, " seconds");

    end;

end;

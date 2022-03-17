
using FileIO
using DelimitedFiles
using Random: seed!

# Experiment values
MinimumReductionsMSE = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7];
MaxNumNodes = 5:5:200;
numFolds = 10;

function loadDataset(datasetName::String; datasetType::DataType=Float64)
    dataset = DelimitedFiles.readdlm(string("datasets/", datasetName, ".tsv"));
    @assert(dataset[1,end]=="target");
    @assert(all(dataset[1,1:end-1].!="target"));
    dataset = datasetType.(dataset[2:end,:]);
    inputs = dataset[:,1:end-1];
    targets = dataset[:, end];
    return (inputs, targets);
end;

using Random:shuffle!

function crossvalidationIndices(N, numFolds)
    # Create the same cross-validation indices used in the rest of experiments
    indicesKFold = repeat(1:numFolds, Int(trunc(N/numFolds)+1))[1:N];
    seed!(1);
    shuffle!(indicesKFold);
    return indicesKFold;
end;

function readFilesInFolder(folder::String; extension::String="")
    fileNames = String[];
    for fileName in readdir(folder)
        if isempty(extension) || ((length(fileName)>length(extension)+1) && (uppercase(fileName[end-length(extension):end]) == uppercase(string(".", extension))))
            push!(fileNames, fileName);
        end;
    end;
    return fileNames;
end;


function datasetNames()
    names = readFilesInFolder("datasets"; extension="tsv")
    return [name[1:end-4] for name in names];
end;

function datasetNameFromFile(name::String)
    i = findfirst("_", name); i = i[1];
    @assert(tryparse(Int, name[1:i-1])!=nothing);
    name = replace(name, "_"=>" ");
    name = name[i+1:end];
    return name;
end;

resultsFileName(datasetName, strategy) = string("results/results_DoME_", datasetName, "_strategy_", strategy, ".jld2");

function parseARGS()
    args = Array{Any,1}(undef, length(ARGS));
    args[:] .= ARGS[:];
    integers = tryparse.(Int,ARGS);
    args[integers.!=nothing] .= integers[integers.!=nothing];
    return args;
end;

function experimentsDoMEFinished(results::Array{<:Real,4})
    if size(results,1)<length(MinimumReductionsMSE) || size(results,2)<length(MaxNumNodes) || (size(results,3)<numFolds)
        return false
    else
        @assert(size(results,3)==numFolds);
        @assert(size(results,4)==5);
        return (!any(isnan.(results[:,:,:,[1,3,4,5]])) && !any(results[:,:,:,3:end].==0))
    end;
end

function experimentsDoMEFinished(datasetName, strategy)
    fileName = resultsFileName(datasetName, strategy);
    !isfile(fileName) && return (false, []);
    results = load(fileName, "results");
    results = results[1:min(size(results,1),length(MinimumReductionsMSE)), 1:min(size(results,2),length(MaxNumNodes)), :, :];
    return (experimentsDoMEFinished(results), results);
end

function experimentsDoMEFinished(datasetName)
    results = Array{Float64,5}(undef, 4, length(MinimumReductionsMSE), length(MaxNumNodes), numFolds, 5);
    for strategy in 1:4
        (finished, resultsThisStrategy) = experimentsDoMEFinished(datasetName, strategy);
        !finished && return (false, []);
        results[strategy, :, :, :, :] .= resultsThisStrategy;
    end;
    return (true, results);
end

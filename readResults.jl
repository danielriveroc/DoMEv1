
using FileIO
using JLD2
using Statistics
using Printf
using StatsPlots

include("usefulFunctions.jl")

WriteLatexString = false;


# Copy-pasted from "Analytic Continued Fractions for Regression: A Memetic Algorithm Approach"
# The only modification done was to fix the dataset names that were abreviated in the paper
Moscato2021Table3 = "ESL 0.274 0.319 0.272 0.268 fri_c3_1000_10 0.077 0.066 0.062 0.068 SWD 0.390 0.405 0.408 0.432 fri_c0_1000_5 0.042 0.074 0.071 0.106 LEV 0.425 0.424 0.422 0.353 fri_c3_100_5 0.250 0.182 0.093 0.125 ERA 2.514 2.585 2.567 2.450 fri_c1_1000_5 0.049 0.047 0.047 0.050 USCrime 3.93e2 2.57e2 3.78e2 2.20e2 fri_c3_250_5 0.126 0.110 0.108 0.086 FacultySalaries 4.035 8.071 4.111 1.277 fri_c4_250_10 0.148 0.173 0.172 0.133 vineyard 6.010 8.223 7.825 4.222 fri_c4_500_50 0.079 0.141 0.114 0.169 auto_price 5.89e6 3.89e6 4.03e6 6.01e6 fri_c3_500_5 0.103 0.085 0.071 0.078 autoPrice 4.17e6 5.29e6 2.87e6 4.83e6 fri_c3_1000_50 0.068 0.082 0.070 0.120 cloud 0.110 0.208 0.144 0.095 fri_c1_1000_25 0.057 0.070 0.067 0.076 elusage 1.35e2 1.99e2 1.37e2 65.797 fri_c0_100_10 0.149 0.319 0.307 0.226 machine_cpu 3.80e3 2.23e3 2.69e3 2.10e3 fri_c2_1000_50 0.063 0.069 0.076 0.118 analcatdata_vehicle 4.14e4 2.41e4 4.20e4 1.55e4 fri_c4_1000_10 0.050 0.074 0.060 0.063 vinnie 2.287 2.860 2.659 1.934 fri_c0_100_5 0.152 0.193 0.164 0.202 pm10 0.640 0.431 0.399 0.621 fri_c2_500_50 0.049 0.101 0.109 0.130 analcatdata_neavote 1.180 0.818 0.917 0.401 fri_c2_500_10 0.064 0.088 0.097 0.081 analcatdata_election2000 4.33e7 3.40e8 7.72e8 5.09e5 fri_c3_1000_5 0.059 0.049 0.048 0.046 pollution 1.87e3 2.19e3 1.67e3 1.42e3 fri_c1_500_5 0.068 0.077 0.074 0.075 no2 0.272 0.227 0.210 0.295 fri_c0_500_25 0.047 0.132 0.127 0.144 analcatdata_apnea2 1.12e6 9.42e5 7.86e5 6.09e5 fri_c2_100_10 0.750 0.321 0.233 0.320 analcatdata_apnea1 8.16e5 9.98e5 5.28e5 6.96e5 fri_c0_250_10 0.047 0.183 0.170 0.259 cpu 1.75e2 2.36e3 8.83e2 1.64e2 fri_c1_500_50 0.091 0.111 0.133 0.116 fri_c0_250_5 0.065 0.162 0.165 0.176 fri_c1_500_10 0.062 0.086 0.069 0.075 fri_c3_500_25 0.064 0.090 0.097 0.081 fri_c2_500_25 0.077 0.097 0.110 0.121 fri_c1_500_25 0.117 0.108 0.107 0.111 fri_c4_250_25 0.157 0.206 0.178 0.110 fri_c1_1000_50 0.066 0.073 0.077 0.104 fri_c3_500_50 0.101 0.132 0.131 0.156 fri_c4_500_25 0.127 0.111 0.095 0.111 fri_c3_500_10 0.062 0.079 0.066 0.075 fri_c3_1000_25 0.056 0.069 0.059 0.064 fri_c1_250_10 0.114 0.123 0.105 0.130 fri_c4_1000_100 0.047 0.087 0.072 0.227 fri_c1_250_50 0.098 0.171 0.154 0.140 fri_c2_1000_25 0.051 0.064 0.070 0.077 fri_c0_500_5 0.048 0.099 0.106 0.124 fri_c0_1000_50 0.044 0.111 0.113 0.141 fri_c0_500_50 0.053 0.131 0.170 0.211 fri_c1_100_10 0.485 0.245 0.154 0.233 fri_c0_100_25 0.259 0.419 0.383 0.365 fri_c4_1000_25 0.058 0.076 0.066 0.081 fri_c0_250_25 0.064 0.209 0.209 0.224 fri_c1_1000_10 0.047 0.051 0.049 0.059 fri_c0_500_10 0.048 0.153 0.138 0.202 fri_c2_100_5 0.509 0.277 0.255 0.259 fri_c1_100_5 0.490 0.208 0.176 0.241 fri_c0_1000_10 0.047 0.096 0.093 0.155 fri_c2_250_10 0.128 0.123 0.116 0.140 fri_c2_250_5 0.104 0.113 0.094 0.086 fri_c3_250_25 0.116 0.194 0.184 0.111 fri_c2_500_5 0.079 0.068 0.073 0.054 sleuth_ex1714 1.42e6 1.57e6 2.31e6 6.83e5 fri_c0_1000_25 0.039 0.092 0.093 0.074 rabe_266 7.107 7.316 3.043 2.643 fri_c2_1000_5 0.046 0.048 0.045 0.044 sleuth_case2002 75.772 56.240 72.363 41.742 fri_c1_250_5 0.084 0.103 0.093 0.093 rmftsa_ladata 3.012 3.511 3.212 2.764 fri_c3_250_10 0.086 0.150 0.120 0.084 visualizing_environmental 9.624 9.798 9.540 4.700 fri_c0_250_50 0.071 0.238 0.252 0.254 sleuth_ex1605 101.824 98.440 91.986 83.998 fri_c4_500_10 0.065 0.084 0.072 0.072 visualizing_galaxy 3.13e2 2.68e2 2.24e2 4.34e2 fri_c2_250_25 0.179 0.173 0.156 0.149 chatfield_4 2.82e2 3.84e2 2.88e2 1.89e2 fri_c2_1000_10 0.069 0.062 0.058 0.073 sleuth_case1202 3.29e3 3.42e3 3.20e3 1.39e3 fri_c4_1000_50 0.051 0.088 0.068 0.092 chscase_geyser1 38.888 39.902 42.887 31.053";
Moscato2021Table3 = split(Moscato2021Table3, ' ', keepempty=false);
@assert(rem(length(Moscato2021Table3),5)==0);
datasetNamesMoscato2021 = Moscato2021Table3[1:5:end];
datasetNamesMoscato2021 = replace.(datasetNamesMoscato2021, "_"=>" ");
resultsMoscato2021 = parse.(Float64, [Moscato2021Table3[2:5:end] Moscato2021Table3[3:5:end] Moscato2021Table3[4:5:end] Moscato2021Table3[5:5:end]])
function readResultsMoscato2021(datasetName)
    index = findfirst(datasetNamesMoscato2021.==datasetName);
    (index==nothing) && return nothing;
    return resultsMoscato2021[index, :];
end;


resultsDoME = Array{Float64,5}[];
datasets = datasetNames();
finishedExperiments = falses(length(datasets));
for numDatasetName in 1:length(datasets)
    datasetName = datasets[numDatasetName];
    (finished, resultsThisDataset) = experimentsDoMEFinished(datasetName);
    resultsMoscato2021 = readResultsMoscato2021(datasetNameFromFile(datasetName));
    if (finished) && (resultsMoscato2021!=nothing)
        finishedExperiments[numDatasetName] = true;
        push!(resultsDoME, resultsThisDataset);
    else
        println(datasetNameFromFile(datasetName), " ", finished ? "" : "DoME experiments not finished", " ", resultsMoscato2021!=nothing ? "" : "no results found in [Moscato2021]");
    end;
end;
datasets = datasets[finishedExperiments];








function findBestResultsParameters(results, training)
    @assert(!any(isnan.(view(results,:,:,:,:,[1,3,4,5]))));
    meanTrainingResults = mean(results[:,:,:,:,1], dims=4); meanTrainingResults = meanTrainingResults[:,:,:,1];
    stdTrainingResults = std(results[:,:,:,:,1], dims=4); stdTrainingResults = stdTrainingResults[:,:,:,1];
    medianTrainingResults = median(results[:,:,:,:,1], dims=4); medianTrainingResults = medianTrainingResults[:,:,:,1];

    meanTestResults = mean(results[:,:,:,:,2], dims=4); meanTestResults = meanTestResults[:,:,:,1];
    meanTestResults[isnan.(meanTestResults)] .= Inf;
    stdTestResults = std(results[:,:,:,:,2], dims=4); stdTestResults = stdTestResults[:,:,:,1];
    medianTestResults = median(results[:,:,:,:,2], dims=4); medianTestResults = medianTestResults[:,:,:,1];
    medianTestResults[isnan.(medianTestResults)] .= Inf;

    meanHeight = mean(results[:,:,:,:,4], dims=4); meanHeight = meanHeight[:,:,:,1];
    meanNumNodes = mean(results[:,:,:,:,5], dims=4); meanNumNodes = meanNumNodes[:,:,:,1];
    meanExecutionTime = mean(results[:,:,:,:,3], dims=4); meanExecutionTime = meanExecutionTime[:,:,:,1];
    medianHeight = median(results[:,:,:,:,4], dims=4); meanHeight = meanHeight[:,:,:,1];
    medianNumNodes = median(results[:,:,:,:,5], dims=4); meanNumNodes = meanNumNodes[:,:,:,1];
    medianExecutionTime = median(results[:,:,:,:,3], dims=4); meanExecutionTime = meanExecutionTime[:,:,:,1];
    if training
        (bestResult, index) = findmin(medianTrainingResults);
    else
        (bestResult, index) = findmin(medianTestResults);
    end;
    return (medianTrainingResults[index], stdTrainingResults[index], medianTestResults[index], stdTestResults[index], index[1], MinimumReductionsMSE[index[2]], MaxNumNodes[index[3]], medianHeight[index], medianNumNodes[index], medianExecutionTime[index])
end;

findBestTrainingResultsParameters(results) = findBestResultsParameters(results, true);
findBestTestResultsParameters(results) = findBestResultsParameters(results, false);


function fixText(text)
    if !WriteLatexString
        text = replace(text, "\\multirow{4}{*}{" => "");
        text = replace(text, "} & " => "\t");
        text = replace(text, "\\textbf{" => "");
        text = replace(text, " & " => "\t");
        text = replace(text, "& " => "\t\t");
        text = replace(text, "\$" => "");
        text = replace(text, "\\tiny" => "");
        text = replace(text, "\\normalsize" => "");
        text = replace(text, "{" => "(");
        text = replace(text, "}" => ")");
        text = replace(text, ") \\\\" => "");
        text = replace(text, " \\\\" => "");
        text = replace(text, "\\%" => "%");
        text = replace(text, " & " => " ");
    else
        text = replace(text, "e+0" => "e");
    end;
    return text;
end;


println("");
println("Table 5");

for numDataset in 1:length(datasets)
    datasetName = datasetNameFromFile(datasets[numDataset]);
    (medianTrainingResult, stdTrainingResults, medianTestResult, stdTestResults, bestStrategy, bestMinimumReductionsMSE, bestMaxNumNodes, bestHeight, bestNumNodes, bestTime) = findBestTestResultsParameters(resultsDoME[numDataset]);
    testResults = resultsDoME[numDataset][bestStrategy, findfirst(MinimumReductionsMSE.==bestMinimumReductionsMSE), findfirst(MaxNumNodes.==bestMaxNumNodes),:,2];
    stability = 100*median(abs.(testResults .- median(testResults))) / median(testResults);
    trainingResults = resultsDoME[numDataset][bestStrategy, findfirst(MinimumReductionsMSE.==bestMinimumReductionsMSE), findfirst(MaxNumNodes.==bestMaxNumNodes),:,1];
    stability = 100*median(abs.(trainingResults .- median(trainingResults))) / median(trainingResults);
    text = @sprintf("%s & %.3g & %.3g & \$10^{-%d}\$ & %.3g & %.3g & %.2f \\%% \\\\",
        datasetName, medianTestResult, bestStrategy, -log10(bestMinimumReductionsMSE), bestMaxNumNodes, bestTime, stability);
    text = fixText(text);
    println(text);
end;



println("");
println("Figures 5 and 6");

strategyBestMSEResults = zeros(4,length(datasets));
strategyBestTimes = zeros(4,length(datasets));
for numDataset in 1:length(datasets)
    datasetName = datasetNameFromFile(datasets[numDataset]);
    for strategy in 1:4
        newResults = copy(resultsDoME[numDataset]);
        newResults[setdiff(1:4,[strategy]), :, :, :, 2] .= Inf;
        (medianTrainingResult, stdTrainingResults, medianTestResult, stdTestResults, bestStrategy, bestMinimumReductionsMSE, bestMaxNumNodes, bestHeight, bestNumNodes, bestTime) = findBestTestResultsParameters(newResults);
        strategyBestMSEResults[strategy,numDataset] = medianTestResult;
        strategyBestTimes[strategy,numDataset] = bestTime;
    end;
end;

println("Read results from ", size(strategyBestMSEResults,2), " datasets");

strategyBestMSEResults ./= minimum(strategyBestMSEResults, dims=1);
boxplotStrategyMSE = boxplot(strategyBestMSEResults[1,:], label = "Exhaustive");
boxplot!(boxplotStrategyMSE, strategyBestMSEResults[2,:], label = "Exhaustive with constant optimization");
boxplot!(boxplotStrategyMSE, strategyBestMSEResults[3,:], label = "Selective with constant optimization");
boxplot!(boxplotStrategyMSE, strategyBestMSEResults[4,:], label = "Selective");
xaxis!(boxplotStrategyMSE, "Strategy", (0.5, 4.5))
yaxis!(boxplotStrategyMSE, "", (1, 1.75))
display(boxplotStrategyMSE);
savefig(boxplotStrategyMSE,string("results/boxplotStrategyMSE.pdf"));

strategyBestTimes ./= minimum(strategyBestTimes, dims=1);
boxplotStrategyTime = boxplot(strategyBestTimes[1,:], label = "Exhaustive");
boxplot!(boxplotStrategyTime, strategyBestTimes[2,:], label = "Exhaustive with constant optimization");
boxplot!(boxplotStrategyTime, strategyBestTimes[3,:], label = "Selective with constant optimization");
boxplot!(boxplotStrategyTime, strategyBestTimes[4,:], label = "Selective");
xaxis!(boxplotStrategyTime, "Strategy", (0.5, 4.5))
yaxis!(boxplotStrategyTime, "", (1, 100))
display(boxplotStrategyTime);
savefig(boxplotStrategyTime,string("results/boxplotStrategyTime.pdf"));



println("");
println("Table 6");

resultsComparison = zeros(length(datasets),5);
for numDataset in 1:length(datasets)
    datasetName = datasetNameFromFile(datasets[numDataset]);
    (medianTrainingResult, stdTrainingResults, medianTestResult, stdTestResults, bestStrategy, bestMinimumReductionsMSE, bestMaxNumNodes, bestHeight, bestNumNodes, bestTime) = findBestTestResultsParameters(resultsDoME[numDataset]);
    resultsComparisonThisDataset = readResultsMoscato2021(datasetName);
    push!(resultsComparisonThisDataset, medianTestResult);
    resultsComparison[numDataset,:] .= resultsComparisonThisDataset;
    function writeResult(numResult)
        result = resultsComparison[numDataset,numResult];
        if result==minimum(resultsComparisonThisDataset)
            return @sprintf(" & \\textbf{%.3g}", result);
        else
            return @sprintf(" & %.3g", result);
        end;
    end;
    text = @sprintf("%s%s%s%s%s%s \\\\",
        datasetName, writeResult(1), writeResult(2), writeResult(3), writeResult(4), writeResult(5));
    text = fixText(text);
    println(text);
end;



println("");
println("Figure 7");

# https://mirkobunse.github.io/CriticalDifferenceDiagrams.jl/dev/
using CriticalDifferenceDiagrams, DataFrames, PGFPlots

df = DataFrame(dataset_name = repeat(datasets, inner=5), SR_method = repeat(["eplex-1m",  "grad-b",  "xg-b",  "CFR",  "DoME"], size(resultsComparison,1)), MSE = resultsComparison'[:]);

criticalDifferencePlot = CriticalDifferenceDiagrams.plot(
    df,
    :SR_method,
    :dataset_name,
    :MSE;
    maximize_outcome=false,
    reverse_x=true,
    title="Critical Difference Diagram"
)
PGFPlots.save("results/criticalDifference.pdf", criticalDifferencePlot)




println("");
println("Figure 8");

resultsComparison ./= minimum(resultsComparison, dims=2);
boxplotMSEMethods = boxplot(resultsComparison[:,1], label = "eplex-1m");
boxplot!(boxplotMSEMethods, resultsComparison[:,2], label = "grad-b");
boxplot!(boxplotMSEMethods, resultsComparison[:,3], label = "xg-b");
boxplot!(boxplotMSEMethods, resultsComparison[:,4], label = "CFR");
boxplot!(boxplotMSEMethods, resultsComparison[:,5], label = "DoME");
yaxis!(boxplotMSEMethods, "", (1, 3))
display(boxplotMSEMethods);
savefig(boxplotMSEMethods,string("results/boxplotMSEMethods.pdf"));




println("");
println("Table 7");

for numDataset in Int64.(1:round(length(datasets)/2))
    datasetName = datasetNameFromFile(datasets[numDataset]);
    (medianTrainingResult, stdTrainingResults, medianTestResult, stdTestResults, bestStrategy, bestMinimumReductionsMSE, bestMaxNumNodes, bestHeight, bestNumNodes, bestTime) = findBestTestResultsParameters(resultsDoME[numDataset]);
    (inputs,) = loadDataset(datasets[numDataset]); numVariables = size(inputs,2);

    text = @sprintf("%s & %g & %g & %g",
        datasetName, bestNumNodes, bestHeight, numVariables);
    otherDatasetIndex = Int64(numDataset + round(length(datasets)/2));
    if otherDatasetIndex<=length(datasets)
        datasetName = datasetNameFromFile(datasets[otherDatasetIndex]);
        (medianTrainingResult, stdTrainingResults, medianTestResult, stdTestResults, bestStrategy, bestMinimumReductionsMSE, bestMaxNumNodes, bestHeight, bestNumNodes, bestTime) = findBestTestResultsParameters(resultsDoME[otherDatasetIndex]);
        (inputs,) = loadDataset(datasets[otherDatasetIndex]); numVariables = size(inputs,2);
        text = @sprintf("%s & %s & %g & %g & %g",
            text, datasetName, bestNumNodes, bestHeight, numVariables);
    else
        text = @sprintf("%s & & & & ", text);
    end;
    text = @sprintf("%s \\\\", text);
    text = fixText(text);
    println(text);
end;





# Finally, let's take some of the best expressions for some datasets
println("");
println("Table 8");

include("DoME.jl");

datasetNamesFinalExpressions = ["1029_LEV", "1089_USCrime", "192_vineyard", "210_cloud", "228_elusage", "519_vinnie", "522_pm10", "523_analcatdata_neavote", "678_visualizing_environmental", "687_sleuth_ex1605", "712_chscase_geyser1"];

# Redefine this function to write constants with up to three decimals
string(node::Constant) = node.semantic>=0 ? @sprintf("%.3f", node.semantic) : @sprintf("(%.3f)", node.semantic);

for datasetNameFinalExpression in datasetNamesFinalExpressions
    numDataset = findfirst(datasets .== datasetNameFinalExpression);
    (_, _, _, _, bestStrategy, bestMinimumReductionsMSE, bestMaxNumNodes, _, _, _) = findBestTestResultsParameters(resultsDoME[numDataset]);
    if bestMaxNumNodes<=80
        (inputs, targets) = loadDataset(datasets[numDataset]);

        (trainingMSE, validationMSE, testMSE, bestTree) = dome(inputs, targets;
            minimumReductionMSE = bestMinimumReductionsMSE ,
            maximumNodes = bestMaxNumNodes ,
            strategy = [Strategy1, Strategy2, Strategy4, Strategy4][bestStrategy]
        );

        if WriteLatexString
            text = latexString(bestTree);
            # Take the first and last parenthesis
            if (text[1:6]=="\\left(") && (text[end-6:end]=="\\right)")
                text = text[7:end-7];
            end;
            println(datasetNameFromFile(datasets[numDataset]), " & \${", text, "}\$ \\\\");
            println("\\midrule");
        else
            text = string(bestTree);
            println(datasetNameFromFile(datasets[numDataset]), ": ", string(bestTree));
        end;

    end;
end;

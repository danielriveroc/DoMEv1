using FileIO
using DelimitedFiles

# Load the dataset and create a matrix with the inputs and a vector for the targets
dataset = DelimitedFiles.readdlm("datasets/561_cpu.tsv");
inputs  = Float64.(dataset[2:end, 1:end-1]);
targets = Float64.(dataset[2:end, end]);

# Load the DoME system
include("DoME.jl");
# Run DoME with these parameters
(trainingMSE, validationMSE, testMSE, bestTree) = dome(inputs, targets;
   minimumReductionMSE = 1e-6,
   maximumNodes = 30,
   strategy = StrategyExhaustive,
   showText = true
);
# Write the expression on screen
println("Best expression found: ", string(bestTree));
println("Best expression found (written in Latex): ", latexString(bestTree));
# If you want to rename the variable names, one of the easiest way is to do something like:
expression = string(bestTree);
expression = replace(expression, "X1" => "vendor");
expression = replace(expression, "X2" => "MYCT");
expression = replace(expression, "X3" => "MMIN");
expression = replace(expression, "X4" => "MMAX");
expression = replace(expression, "X5" => "CACH");
expression = replace(expression, "X6" => "CHMIN");
expression = replace(expression, "X7" => "CHMAX");
println("Best expression found (with the real names of the variables): ", expression);

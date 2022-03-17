
using Random;
using DelimitedFiles;
using Random: seed!

include("DoME.jl");


seed!(1); # Make the experiment repeatable

numData = 1000;

inputs = Array{Float64,2}(undef, numData, 3);
targets = Array{Float64,1}(undef, numData);
G = 6.67392e-11; # N * (m/kg)^2

minimumMass = 1e23; maximumMass = 1e25;
minimumDist = 1e8;  maximumDist = 1e12;

planet1Mass = 10. .^(log10(minimumMass) .+ rand(numData).*(log10(maximumMass)-log10(minimumMass)));
planet2Mass = 10. .^(log10(minimumMass) .+ rand(numData).*(log10(maximumMass)-log10(minimumMass)));
distance    = 10. .^(log10(maximumDist) .+ rand(numData).*(log10(minimumDist)-log10(maximumDist)));

inputs = [planet1Mass planet2Mass distance];
targets     = G.*planet1Mass.*planet2Mass./(distance.^2);


# As explained in the paper, with so high differences in the variables ranges,
#  in irder to avoid rounding errors, the goal is the mean of the targets
(trainingMSE, _, _, bestTree) = dome(inputs, targets;
    goalMSE = mean(targets) ,
    strategy = Strategy1 ,
    showText = true
    );

expression = string(bestTree);
expression = replace(expression, "X1" => "M1");
expression = replace(expression, "X2" => "M2");
expression = replace(expression, "X3" => "d");
println("Best expression found: ", expression);


########################################################################################################
# A different execution: without so big differences in inputs and targets (not documented in the paper)

println("---------------------------------------------------------------------------------------------")
println("Additional exucution with lower values in masses and distances (not documented in the paper):")

maximumMass/minimumMass
planet1Mass ./= minimumMass;
planet2Mass ./= minimumMass;
distance ./= minimumDist;
G = G*minimumMass/minimumDist;

inputs = [planet1Mass planet2Mass distance];
targets     = G.*planet1Mass.*planet2Mass./(distance.^2);

# In this case, since no so big rounding errors are expected, the MSE goal can be 0
# However, some small rounding erros might happen. For this reason, it is better to have a small MSE goal (1e-10)

(trainingMSE, _, _, bestTree) = dome(inputs, targets;
    goalMSE = 1e-10 ,
    strategy = Strategy1 ,
    showText = true
    );

expression = string(bestTree);
expression = replace(expression, "X1" => "M1");
expression = replace(expression, "X2" => "M2");
expression = replace(expression, "X3" => "d");
println("Best expression found: ", expression);

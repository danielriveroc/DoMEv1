
include("DoME.jl");
include("usefulFunctions.jl");

# Uncomment the desired line to perform the experiments and obtain the results shown in Table 5
datasetName = "1027_ESL";                      Strategy = Strategy2; MinimumReductionMSE = 1e-5; MaximumNodes = 100;
# datasetName = "1028_SWD";                      Strategy = Strategy1; MinimumReductionMSE = 1e-5; MaximumNodes = 60;
# datasetName = "1029_LEV";                      Strategy = Strategy3; MinimumReductionMSE = 1e-4; MaximumNodes = 45;
# datasetName = "1030_ERA";                      Strategy = Strategy2; MinimumReductionMSE = 1e-7; MaximumNodes = 100;
# datasetName = "1089_USCrime";                  Strategy = Strategy3; MinimumReductionMSE = 1e-3; MaximumNodes = 20;
# datasetName = "1096_FacultySalaries";          Strategy = Strategy1; MinimumReductionMSE = 1e-7; MaximumNodes = 155;
# datasetName = "192_vineyard";                  Strategy = Strategy4; MinimumReductionMSE = 1e-3; MaximumNodes = 15;
# datasetName = "195_auto_price";                Strategy = Strategy1; MinimumReductionMSE = 1e-7; MaximumNodes = 65;
# datasetName = "207_autoPrice";                 Strategy = Strategy2; MinimumReductionMSE = 1e-4; MaximumNodes = 160;
# datasetName = "210_cloud";                     Strategy = Strategy4; MinimumReductionMSE = 1e-7; MaximumNodes = 55;
# datasetName = "228_elusage";                   Strategy = Strategy1; MinimumReductionMSE = 1e-3; MaximumNodes = 15;
# datasetName = "230_machine_cpu";               Strategy = Strategy4; MinimumReductionMSE = 1e-7; MaximumNodes = 70;
# datasetName = "485_analcatdata_vehicle";       Strategy = Strategy4; MinimumReductionMSE = 1e-5; MaximumNodes = 145;
# datasetName = "519_vinnie";                    Strategy = Strategy3; MinimumReductionMSE = 1e-4; MaximumNodes = 5;
# datasetName = "522_pm10";                      Strategy = Strategy2; MinimumReductionMSE = 1e-4; MaximumNodes = 30;
# datasetName = "523_analcatdata_neavote";       Strategy = Strategy3; MinimumReductionMSE = 1e-2; MaximumNodes = 10;
# datasetName = "527_analcatdata_election2000";  Strategy = Strategy3; MinimumReductionMSE = 1e-5; MaximumNodes = 100;
# datasetName = "542_pollution";                 Strategy = Strategy4; MinimumReductionMSE = 1e-4; MaximumNodes = 30;
# datasetName = "547_no2";                       Strategy = Strategy1; MinimumReductionMSE = 1e-4; MaximumNodes = 50;
# datasetName = "556_analcatdata_apnea2";        Strategy = Strategy4; MinimumReductionMSE = 1e-4; MaximumNodes = 115;
# datasetName = "557_analcatdata_apnea1";        Strategy = Strategy1; MinimumReductionMSE = 1e-6; MaximumNodes = 125;
# datasetName = "561_cpu";                       Strategy = Strategy3; MinimumReductionMSE = 1e-7; MaximumNodes = 140;
# datasetName = "659_sleuth_ex1714";             Strategy = Strategy3; MinimumReductionMSE = 1e-2; MaximumNodes = 40;
# datasetName = "663_rabe_266";                  Strategy = Strategy1; MinimumReductionMSE = 1e-7; MaximumNodes = 185;
# datasetName = "665_sleuth_case2002";           Strategy = Strategy3; MinimumReductionMSE = 1e-6; MaximumNodes = 30;
# datasetName = "666_rmftsa_ladata";             Strategy = Strategy1; MinimumReductionMSE = 1e-3; MaximumNodes = 70;
# datasetName = "678_visualizing_environmental"; Strategy = Strategy3; MinimumReductionMSE = 1e-5; MaximumNodes = 30;
# datasetName = "687_sleuth_ex1605";             Strategy = Strategy2; MinimumReductionMSE = 1e-7; MaximumNodes = 30;
# datasetName = "690_visualizing_galaxy";        Strategy = Strategy1; MinimumReductionMSE = 1e-7; MaximumNodes = 75;
# datasetName = "695_chatfield_4";               Strategy = Strategy3; MinimumReductionMSE = 1e-6; MaximumNodes = 30;
# datasetName = "706_sleuth_case1202";           Strategy = Strategy3; MinimumReductionMSE = 1e-3; MaximumNodes = 70;
# datasetName = "712_chscase_geyser1";           Strategy = Strategy1; MinimumReductionMSE = 1e-3; MaximumNodes = 45;

# Load the dataset
(inputs, targets) = loadDataset(datasetName);

numFolds = 10;
# Create the same cross-validation indices used in the rest of experiments
indicesKFold = crossvalidationIndices(length(targets), numFolds);

testValues = Array{Float64,1}(undef,numFolds);

println("Dataset \"", datasetNameFromFile(datasetName), "\" with strategy ", Strategy, ", Min. MSE reduction ", MinimumReductionMSE, " and Maximum num. nodes ", MaximumNodes);

for numFold = 1:numFolds

    trainingInputs = inputs[indicesKFold.!=numFold,:];
    trainingTargets = targets[indicesKFold.!=numFold];

    (_, _, _, bestTree) = dome(trainingInputs, trainingTargets;
        minimumReductionMSE = MinimumReductionMSE ,
        maximumNodes = MaximumNodes ,
        strategy = Strategy ,
        showText = false
    );

    # Get this expression as a string with vector operations
    expr = vectorString(bestTree);
    # Convert the variable names "X" to "inputs"
    expr = replace(expr, "X" => "inputs");
    # Evaluate this expression
    outputs = eval(Meta.parse(expr))
    # Calculate MSE in the test samples
    testMSE = mean(((targets .- outputs).^2)[indicesKFold.==numFold]);

    println("   Finished fold $numFold/$numFolds, MSE in test: $testMSE");

    testValues[numFold] = testMSE;

end;

println("Median test MSE: ", median(testValues));

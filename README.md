# DoMEv1

This project contains the source code of the first version of the DoME algorithm for Symbolic Regression, written in Julia. The aim of this code is to be able to repeat the experiments described in the paper. For this reason, this code will not be updated. Instead, a new project will be created with the DoME source code that will incorporate the updates. However, this code is fully functional. Feel free to use this source code to perform your experiments. However, if any publication is generated through this system, please add a citation to the following paper:

# How to use DoME

The easiest way to wun DoME is by calling the function dome. Here is an example of use, in which only the main hyperparameters are set:

	using FileIO
	using DelimitedFiles
	
	# Load the dataset and create a matrix with the inputs and a vector for the targets
	dataset = DelimitedFiles.readdlm("datasets/561_cpu.tsv");
	inputs  = Float64.(dataset[2:end, 1:end-1]);
	targets = Float64.(dataset[2:end, end]);

	# Load the DoME system
	include("DoME.jl");
	# Run DoME with this parameters
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

When calling the function dome, inputs is a NxP matrix of real numbers, and targets is a N-length vector or real numbers (N: number of instances, P: number of variables). Inputs and targets can have Float32 or Float64 values; however, since many constants are generated during the run of the algorithm, it is recommended to use Float64 to have the highest precision. Also, the elements of both inputs and targets must have the same type (Float32 or Float64). The parameters minimumReductionMSE, maximumNodes and strategy are the 3 hyperparameters described in the paper.

The declaration of this function is the following, with the whole set of parameters and their default values:

	function dome(inputs::Array{<:Real,2}, targets::Array{<:Real,1};
	    dataInRows = true,
	    validationIndices::Array{Int64,1} = Array{Int64,1}([]),
	    testIndices::Array{Int64,1} = Array{Int64,1}([]),
	    validationRatio = 0. ,
	    testRatio = 0. ,
	    minimumReductionMSE = 1e-6,
	    maximumHeight = Inf ,
	    maximumNodes = 50 ,
	    strategy = StrategyExhaustive ,
	    goalMSE = 0 ,
	    maxIterations = Inf ,
	    showText = false ,
	    checkForErrors = false
	    )

The description of these parameters is the following, grouped in

	dataInRows -> allows the input matrix to have dimensions NxP when it is set to true (by default) or PxN when it is false (N: number of instances).
	minimumReductionMSE -> A search is found to be successful if the reduction in MSE is positive and higher than the previous MSE value multiplied by this parameter.
	maximumHeight -> maximum height of the tree. As explained in the paper, this parameter is not recommended to be used in the experiments.
	maximumNodes -> maximum number of nodes in the tree.
	strategy -> the strategy used to select which searches are going to be performed on which nodes. The 4 strategies described in the paper are available, with names StrategyExhaustive (by default), StrategyExhaustiveWithConstantOptimization, StrategySelectiveWithConstantOptimization and StrategySelective. They are also called Strategy1, Strategy2, Strategy3, Strategy4 respectively as used in the paper.
	goalMSE -> if the algorithm reaches this MSE value in training, the iterative process is stopped.
	maxIterations -> maximum number of iterations to be performed.
	testIndices -> allows to split the dataset by separating some instances to perform the test, specifying which ones will be used for test.
	testRatio ->   allows to split the dataset by separating some random instances to perform the test, specifying the ratio used for test.
	validationIndices -> allows to split the dataset by separating some instances to perform the validation, specifying which ones will be used for validation.
	validationRatio ->   allows to split the dataset by separating some random instances to perform the validation, specifying the ratio used for validation.
	showText -> if it is set to true, on each iteration some text (iteration number, best tree, MSE in training and test) is shown.
	checkForErrors -> this parameter was used only for debugging, to easily find bugs in the source code. Therefore, it is best to leave it as false.

As it can be sen, this function allows the definition of a validation set. However, it was not used in the experiments of the paper and thus this part of the code may have errors.

To run this, you need to have the following files in the same folder:

	- DoME.jl -> Main algorithm
	- Equation.jl
	- NodePool.jl
	- Tree.jl -> Tree structure creation and some useful functions to operate

An alternative way to run DoME is by creating a DoME struct and calling the function Step! for each iteration. This is automatically done by the previous way to run DoME.

# How to repeat the experiments described in the paper

To obtain the values of tables 1 and 3 and 4 in the paper, run the file examplePaper.jl

To obtain the values of tables 4 in the paper (Newton's law of universal gravitation), run the file experimentNewton.jl

To repeat the experiments describen in section 4, create a folder named "datasets" and store the corresponding datasets from PMLB. When running the experiments, a folder named "results" will be created. The following files were used to perform the experiments:

	- usefulFunctions.jl
	- createScripts.jl
	- experimentsDoME.jl
	- readResults.jl

# How to define your own strategy

Strategies are based on calling the functions PerformSearches! and OptimizeConstants!

The function PerformSearches! allows to specify in the call which nodes are going to be used on each search. To do this, this function has keyword parameters for each search, and in each one the user can specify type of the nodes in which this search will take place. These types are Terminal, Variable, Constant, and NonTerminal. Also, the types Any and Nothing can be used to specify that a search will be performed on all of the nodes, or in none of them respectively. The declaration of this function is the following:

	function PerformSearches!(obj::DoME;
	   whichNodesPerformConstantSearch        ::Union{DataType,Union} = Nothing ,
	   whichNodesPerformVariableSearch        ::Union{DataType,Union} = Nothing ,
	   whichNodesPerformConstantVariableSearch::Union{DataType,Union} = Nothing ,
	   performConstantExpressionSearch        ::Bool = false)

Note that constant-expression search receives a boolean value, because this search is only performed on non-terminal nodes.

This function returns a Boolean value: if it is true, a search has been succesful, otherwise no search was succesful. The strategy function to be defined should also return a Boolean value, with the same interpretation.

An example is the Exhaustive strategy, in which the searches are performed on all of the nodes of the tree:

	function StrategyExhaustive(obj::DoME)
	   changeDone = PerformSearches!(obj;
	      whichNodesPerformConstantSearch=Any ,
	      whichNodesPerformVariableSearch=Any ,
	      whichNodesPerformConstantVariableSearch=Any ,
	      performConstantExpressionSearch=true);
	   return changeDone;
	end;

Another example is the Selective strategy, that performs the searches performs searches sequentially, moving on to the next one only if the previous one has been unsuccessful:

	function StrategySelective(obj::DoME)
	   changeDone =               PerformSearches!(obj; whichNodesPerformConstantSearch=Constant);
	   # Variable search only on constants
	   changeDone = changeDone || PerformSearches!(obj; whichNodesPerformVariableSearch=Constant);
	   changeDone = changeDone || PerformSearches!(obj; performConstantExpressionSearch=true);
	   # Constant-variable search only on terminals
	   changeDone = changeDone || PerformSearches!(obj; whichNodesPerformConstantVariableSearch=Union{Constant,Variable});
	   if (!changeDone)
	      # Constant search on variables and non-terminals, variable seach on variables and non-terminals, and constant-variable search on non-terminals
	      changeDone = PerformSearches!(obj;
	         whichNodesPerformConstantSearch=Union{Variable,NonTerminal} ,
	         whichNodesPerformVariableSearch=Union{Variable,NonTerminal} ,
	         whichNodesPerformConstantVariableSearch=NonTerminal );
	   end;
	   return changeDone;
	end;

	function StrategyExhaustiveWithConstantOptimization(obj::DoME)
	   changeDone = PerformSearches!(obj;
	      whichNodesPerformConstantSearch=Union{Variable,NonTerminal} ,
	      whichNodesPerformVariableSearch=Any ,
	      whichNodesPerformConstantVariableSearch=Any ,
	      performConstantExpressionSearch=true);
	      changeDone && OptimizeConstants(obj);
	   return changeDone;
	end;

	function StrategySelectiveWithConstantOptimization(obj::DoME)
	   # Variable search only on constants
	   changeDone =               PerformSearches!(obj; whichNodesPerformVariableSearch=Constant);
	   changeDone = changeDone || PerformSearches!(obj; performConstantExpressionSearch=true);
	   # Constant-variable search only on terminals
	   changeDone = changeDone || PerformSearches!(obj; whichNodesPerformConstantVariableSearch=Union{Constant,Variable});
	   if (!changeDone)
	      # Constant search on variables and non-terminals, variable seach on variables and non-terminals, and constant-variable search on non-terminals
	      changeDone = PerformSearches!(obj;
	         whichNodesPerformConstantSearch=Union{Variable,NonTerminal} ,
		 whichNodesPerformVariableSearch=Union{Variable,NonTerminal} ,
		 whichNodesPerformConstantVariableSearch=NonTerminal );
	   end;
	   changeDone && OptimizeConstants(obj);
	   return changeDone;
	end;

These four strategies are provided in the file DoME.jl and are available for use.

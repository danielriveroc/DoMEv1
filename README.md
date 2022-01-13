# DoMEv1

This project contains the source code of the first version of the DoME algorithm for Symbolic Regression, written in Julia. The aim of this code is to be able to repeat the experiments described in the paper. For this reason, this code will not be updated. Instead, a new project will be created with the DoME source code that will incorporate the updates. However, this code is fully functional. Feel free to use this source code to perform your experiments. However, if any publication is generated through this system, please add a citation to the following paper:

# How to use DoME

The easiest way to wun DoME is by calling the function dome. Here is an example of use, in which only the main hyperparameters are set:

	(trainingMSE, validationMSE, testMSE, bestTree) = dome(inputs, targets;
	   minimumReductionMSE = 1e-6,
	   maximumNodes = 50 ,
	   strategy = StrategyExhaustive
	);
	println(string(bestTree));

where inputs is a NxP matrix of real numbers, and targets is a N-length vector or real numbers. Inputs and targets can have Float32 or Float64 values; however, since many constants are generated during the run of the algorithm, it is recommended to use Float64 to have the highest precision. Also, the elements of both inputs and targets must have the same type (Float32 or Float64). The parameters minimumReductionMSE, maximumNodes and strategy are the 3 hyperparameters described in the paper.

The declaration of this function is the following, with the whole set of parameters:

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

The parameter dataInRows allows the input matrix to have dimensions NxP when it is set to true (by default) or PxN when it is false (N: number of instances).

Regarding strategies, the 4 strategies described in the paper are available, with names StrategyExhaustive (by default), StrategyExhaustiveWithConstantOptimization, StrategySelectiveWithConstantOptimization and StrategySelective. They are also called Strategy1, Strategy2, Strategy3, Strategy4 respectively as used in the paper.

You may see that the source code allows the definition of a validation set. However, it was not used in the experiments of the paper and thus this part of the code may have errors.

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

# How to repeat the experiments described in the paper

	function StrategyExhaustive(obj::DoME)
	   changeDone = PerformSearches!(obj;
	      whichNodesPerformConstantSearch=Any ,
	      whichNodesPerformVariableSearch=Any ,
	      whichNodesPerformConstantVariableSearch=Any ,
	      performConstantExpressionSearch=true);
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



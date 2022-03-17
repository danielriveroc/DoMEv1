
using Statistics

include("NodePool.jl")

MSE(v1, v2) = mean((v1.-v2).^2);

function createVariableValues(inputs::Array{<:Real,2}; dataInRows=true)
    if (dataInRows)
        inputValues = [inputs[:,numInput] for numInput in 1:size(inputs,2)];
    else
        inputValues = [inputs[numInput,:] for numInput in 1:size(inputs,1)];
    end;
    return inputValues;
end;


mutable struct DoME
    tree::Tree
    mse::Real
    goalMSE::Real
    nodePoolVariables::NodePool
    inputs::Array{<:Real,2}
    targets::Array{<:Real,1}
    initialEquation::NodeEquation
    minimumReductionMSE::Real
    maximumHeight::Number
    maximumNodes::Number
    strategy::Function
    checkForErrors::Bool
    function DoME(inputs::Array{<:Real,2}, targets::Array{<:Real,1};
        dataInRows=true ,
        maximumNodes = Inf ,
        # initialTreeOrder = 1 ,
        checkForErrors = false
    )

        @assert(!any(isnan.(inputs)));
        @assert(!any(isinf.(inputs)));
        @assert(!anyNaN(targets));
        @assert(!anyInf(targets));

        @assert(eltype(inputs)==eltype(targets));
        @assert(eltype(inputs)==Float32 || eltype(inputs)==Float64);

        nodePoolVariables = NodePool();

        if (dataInRows)
            @assert (length(targets)==size(inputs,1))
        else
            @assert (length(targets)==size(inputs,2))
        end;

        # Create variables, add them to the node pool and to the prototypes
        variableValues = createVariableValues(inputs; dataInRows=dataInRows);
        for numVariable in 1:length(variableValues)
            # if (length(unique(variableValues[numVariable]))>1)
            if (!allequal(variableValues[numVariable]))
                nodeVariable = Variable(UInt(numVariable); semantic=variableValues[numVariable]);
                addNodeToPool!(nodePoolVariables, nodeVariable, UInt(1), UInt(1));
            end;
        end;
        if (length(nodePoolVariables.nodes)==0)
            error("No variables have been added");
        end;

        # Build initial tree
        tree = Constant(mean(targets));

        mse = MSE(targets, evaluateTree(tree));

        initialEquation = NodeEquation(eltype(inputs)(1.), targets, eltype(inputs)(0.), eltype(inputs)(-1.));

        maximumHeight = Inf;
        strategy(obj::DoME) = PerformSearches!(obj; whichNodesPerformConstantSearch=Any, whichNodesPerformVariableSearch=Any, whichNodesPerformConstantVariableSearch=Any, performConstantExpressionSearch=true);

        return new(tree, mse, eltype(inputs)(0), nodePoolVariables, inputs, targets, initialEquation, eltype(inputs)(1e-6), maximumHeight, maximumNodes, strategy, checkForErrors);
    end;
end;










#################################################################################################################################################################
#
# The four searches and the function to optimize constants
#

function constantSearch(nodes::Array{Tree,1}, mse::Real, numNodes::Array{UInt,1}, whichNodesPerformConstantSearch::Union{DataType,Union}; checkForErrors=false)
    bestReductionMSE = -Inf;
    bestNumNodes = -Inf;
    constantBest = -Inf;
    numNodeBest = -1;
    for numNode in 2:length(nodes)
        node = nodes[numNode];
        if isa(node, whichNodesPerformConstantSearch) && !isNaN(node.equation.S)
            (constant, reduction) = calculateConstantMinimizeEquation(node.equation, mse; checkForErrors=checkForErrors);
            if (reduction>bestReductionMSE)
                constantBest = constant;
                bestReductionMSE = reduction;
                numNodeBest = numNode;
            end;
        end;
    end;
    return (bestReductionMSE, numNodeBest, constantBest);
end;


function variableSearch(nodes::Array{Tree,1}, nodePoolVariables::NodePool, mse::Real, whichNodesPerformVariableSearch::Union{DataType,Union}; checkForErrors=false)

    maxReductionMSE = -Inf;
    indexNodePoolMaxReductionMSE = nothing;
    numNodeMaxReductionMSE = -1;

    if (checkForErrors)
        @assert(unique(nodePoolVariables.heights)==[1]);
        @assert(unique(nodePoolVariables.numNodes)==[1]);
        @assert(unique(isa.(nodePoolVariables.nodes,Variable))==[true]);
    end;

    for numNode in 1:length(nodes)

        node = nodes[numNode];

        if isa(node, whichNodesPerformVariableSearch) && !isNaN(node.equation.S)

            for numVar in 1:length(nodePoolVariables.nodes)

                # If looking into the variables of the tree, do not look into the same variables
                isa(node,Variable) && (node.variableNumber==nodePoolVariables.nodes[numVar].variableNumber) && continue;

                semantic = nodePoolVariables.semantics[numVar];

                # If the semantic of the variable is not in the domain, do not check it
                semanticNotInDomain(semantic, node.equation) && continue;

                reduction = calculateMSEReduction(semantic, node.equation, mse; checkForErrors=checkForErrors);
                if (reduction>maxReductionMSE)
                    maxReductionMSE = reduction;
                    indexNodePoolMaxReductionMSE = numVar;
                    numNodeMaxReductionMSE = numNode;
                end;
            end;
        end;
    end;

    return (maxReductionMSE, numNodeMaxReductionMSE, indexNodePoolMaxReductionMSE);

end;




function constantVariableSearch(nodes::Array{Tree,1}, nodePool::NodePool, mse::Real, maximumHeights, maximumNodes, whichNodesPerformConstantVariableSearch::Union{DataType,Union}; checkForErrors=false)

    maxReductionMSE = -Inf;
    operationMaxReductionMSE = nothing;
    constantMaxReductionMSE = -Inf;
    variableMaxReductionMSE = nothing;
    numNodeMaxReductionMSE = -1;

    for numNode in 1:length(nodes)

        node = nodes[numNode];

        if isa(node, whichNodesPerformConstantVariableSearch) && !isNaN(node.equation.S)

            # (a,b,c,d) = node.equation;
            eq=node.equation; a=eq.a; b=eq.b; c=eq.c; d=eq.d; S=eq.S;

            for numVar = findall(((nodePool.numNodes.+2).<=maximumNodes[numNode]) .& (nodePool.heights.<=(maximumHeights[numNode]-1)));

                semantic = nodePool.nodes[numVar].semantic;

                if ( (!isa(node,BinaryNode)) || (node.name!="+") || (!isa(node.child1,Constant)) || (!isa(node.child2,Variable)) || (node.child2.variableNumber!=numVar) )
                    equationConstant = equationAddChild1(semantic, a, b, c, d, S);
                    if (checkForErrors) checkEquation(equationConstant); end;
                    # If there are possible values in the domain
                    if (!isNaN(equationConstant.S))
                        (constant, reductionMSE) = calculateConstantMinimizeEquation(equationConstant, mse; checkForErrors=checkForErrors);
                        if (reductionMSE>maxReductionMSE)
                            maxReductionMSE = reductionMSE;
                            constantMaxReductionMSE = constant;
                            operationMaxReductionMSE = AddNode;
                            numNodeMaxReductionMSE = numNode;
                            variableMaxReductionMSE = nodePool.nodes[numVar];
                        end;
                    end;
                end;

                if ( (!isa(node,BinaryNode)) || (node.name!="-") || (!isa(node.child1,Constant)) || (!isa(node.child2,Variable)) || (node.child2.variableNumber!=numVar) )
                    equationConstant = equationSubChild1(semantic, a, b, c, d, S);
                    # b2 = cleanEquation(b2);
                    # d2 = cleanEquation(d2);
                    if (checkForErrors) checkEquation(equationConstant); end;
                    # If there are possible values in the domain
                    if (!isNaN(equationConstant.S))
                        (constant, reductionMSE) = calculateConstantMinimizeEquation(equationConstant, mse; checkForErrors=checkForErrors);
                        if (reductionMSE>maxReductionMSE)
                            maxReductionMSE = reductionMSE;
                            constantMaxReductionMSE = constant;
                            operationMaxReductionMSE = SubNode;
                            numNodeMaxReductionMSE = numNode;
                            variableMaxReductionMSE = nodePool.nodes[numVar];
                        end;
                    end;
                end;

                if ( (!isa(node,BinaryNode)) || (node.name!="*") || (!isa(node.child1,Constant)) || (!isa(node.child2,Variable)) || (node.child2.variableNumber!=numVar) )
                    equationConstant = equationMulChild1(semantic, a, b, c, d, S);
                    if (checkForErrors) checkEquation(equationConstant) end;
                    # If there are possible values in the domain
                    if (!isNaN(equationConstant.S))
                        (constant, reductionMSE) = calculateConstantMinimizeEquation(equationConstant, mse; checkForErrors=checkForErrors);
                        if (reductionMSE>maxReductionMSE)
                            maxReductionMSE = reductionMSE;
                            constantMaxReductionMSE = constant;
                            operationMaxReductionMSE = MulNode;
                            numNodeMaxReductionMSE = numNode;
                            variableMaxReductionMSE = nodePool.nodes[numVar];
                        end;
                    end;
                end;

                if ( (!isa(node,BinaryNode)) || (node.name!="/") || (!isa(node.child1,Constant)) || (!isa(node.child2,Variable)) || (node.child2.variableNumber!=numVar) )
                    if !(nodePool.hasZerosInSemantics[numVar])
                        equationConstant = equationDivChild1(semantic, a, b, c, d, S);
                        # b2 = cleanEquation(b2);
                        # d2 = cleanEquation(d2);
                        if (checkForErrors) checkEquation(equationConstant); end;
                        # If there are possible values in the domain
                        if (!isNaN(equationConstant.S))
                            (constant, reductionMSE) = calculateConstantMinimizeEquation(equationConstant, mse; checkForErrors=checkForErrors);
                            if (reductionMSE>maxReductionMSE)
                                maxReductionMSE = reductionMSE;
                                constantMaxReductionMSE = constant;
                                operationMaxReductionMSE = DivNode;
                                numNodeMaxReductionMSE = numNode;
                                variableMaxReductionMSE = nodePool.nodes[numVar];
                            end;
                        end;
                    end;
                end;

            end;
        end;
    end;

    return (maxReductionMSE, numNodeMaxReductionMSE, operationMaxReductionMSE, constantMaxReductionMSE, variableMaxReductionMSE);
end;



function constantExpressionSearch(nodes::Array{Tree,1}, mse::Real, heights, maximumHeights, numNodes, maximumNodes; checkForErrors=false)

    maxReductionMSE = -Inf;
    subTree = nothing;
    operationMaxReductionMSE = nothing;
    constantMaxReductionMSE = -Inf;
    numNodeMaxReductionMSE = -1;

    # Search only on non terminal nodes
    for numNode in findall(isa.(nodes,NonTerminal) .& (heights.<maximumHeights) .& ((numNodes.+2).<=maximumNodes))
        node = nodes[numNode];
        semantic = node.semantic;
        eq=node.equation; a=eq.a; b=eq.b; c=eq.c; d=eq.d; S=eq.S;

        if !isNaN(S)

            if ( ((node.name!="+") && (node.name!="-")) || ((!isa(node.child1,Constant)) || (!isa(node.child2,Constant))) )
                equationConstant = equationAddChild1(node.semantic, a, b, c, d, S);
                if (checkForErrors) checkEquation(equationConstant); end;
                # If there are possible values in the domain
                if (!isNaN(equationConstant.S))
                    (constant, reductionMSE) = calculateConstantMinimizeEquation(equationConstant, mse; checkForErrors=checkForErrors);
                    if (reductionMSE>maxReductionMSE)
                        maxReductionMSE = reductionMSE;
                        constantMaxReductionMSE = constant;
                        operationMaxReductionMSE = AddNode;
                        numNodeMaxReductionMSE = numNode;
                    end;
                end;
            end

            if ( ((node.name!="*") && (node.name!="/")) || ((!isa(node.child1,Constant)) || (!isa(node.child2,Constant))) )
                equationConstant = equationMulChild1(node.semantic, a, b, c, d, S);
                if (checkForErrors) checkEquation(equationConstant); end;
                # If there are possible values in the domain
                if (!isNaN(equationConstant.S))
                    (constant, reductionMSE) = calculateConstantMinimizeEquation(equationConstant, mse; checkForErrors=checkForErrors);
                    if (reductionMSE>maxReductionMSE)
                        maxReductionMSE = reductionMSE;
                        constantMaxReductionMSE = constant;
                        operationMaxReductionMSE = MulNode;
                        numNodeMaxReductionMSE = numNode;
                    end;
                end;
            end;

        end;
    end;

    return (maxReductionMSE, numNodeMaxReductionMSE, operationMaxReductionMSE, constantMaxReductionMSE);
end;




function OptimizeConstants!(obj::DoME)

    nodes, = iterateTree(obj.tree);

    constantNodes = [isa(node,Constant) for node in nodes];
    nodes = nodes[constantNodes];
    isempty(nodes) && return nothing;
    pathsToConstantNodes = findNode.([obj.tree], nodes);

    if (length(nodes)==1)
        calculateEquations!(obj.tree, obj.initialEquation, path=pathsToConstantNodes[1]; checkForErrors=obj.checkForErrors);
        node = nodes[1];
        minimumReduction = obj.minimumReductionMSE * obj.mse;
        (constant, reduction) = calculateConstantMinimizeEquation(node.equation, obj.mse; checkForErrors=obj.checkForErrors);
        clearEquations!(obj.tree);
        if (reduction>minimumReduction) && (constant!=nodes[1].semantic)
            node.semantic = constant;
            mse = MSE(obj.targets, reevaluatePath(obj.tree, pathsToConstantNodes[1]; checkForErrors=obj.checkForErrors));
            if (obj.checkForErrors)
                @assert(isa(node, Constant))
                @assert(!isnan(mse));
                @assert(equal(mse+reduction, obj.mse));
            end;
            obj.mse = mse;
        end;
        return nothing;
    end;

    numNode = 0;
    numIterationsWithNoImprovement = 0;
    while (true)

        numNode = (numNode>=length(nodes)) ? 1 : numNode+1;

        calculateEquations!(obj.tree, obj.initialEquation; path=pathsToConstantNodes[numNode], checkForErrors=obj.checkForErrors)

        if (obj.checkForErrors)
            checkEquations(obj.tree, obj.mse);
        end;

        minimumReduction = obj.minimumReductionMSE * obj.mse;
        (constant, reduction) = calculateConstantMinimizeEquation(nodes[numNode].equation, obj.mse; checkForErrors=obj.checkForErrors);

        clearEquations!(obj.tree);
        if (reduction>minimumReduction) && (constant!=nodes[numNode].semantic)

            oldConstant = nodes[numNode].semantic;
            nodes[numNode].semantic = constant;

            mse = MSE(obj.targets, reevaluatePath(obj.tree, pathsToConstantNodes[numNode]; checkForErrors=obj.checkForErrors));

            # Check that the resulting MSE is lower than the previous MSE in obj.minimumReductionMSE
            # This is expected to be. However, precision errors may happen
            if (mse >= (obj.mse - (obj.minimumReductionMSE * obj.mse)))

                # Undo the change in the tree
                nodes[numNode].semantic = oldConstant;
                mse = MSE(obj.targets, reevaluatePath(obj.tree, pathsToConstantNodes[numNode]; checkForErrors=obj.checkForErrors));
                obj.checkForErrors && @assert(mse==obj.mse);
                numIterationsWithNoImprovement += 1;
                if (numIterationsWithNoImprovement>=length(nodes))
                    return;
                end;

            else

                if (obj.checkForErrors)
                    @assert(isa(nodes[numNode], Constant))
                    @assert(!isnan(mse));
                    @assert(equal(mse+reduction, obj.mse));
                end;
                obj.mse = mse;
                numIterationsWithNoImprovement = 0;

            end;
        else
            numIterationsWithNoImprovement += 1;
            if (numIterationsWithNoImprovement>=length(nodes))
                return;
            end;
        end;
    end;

end;



#################################################################################################################################################################
#
# The function that performs the searches, called by a strategy
#

function PerformSearches!(obj::DoME;
    whichNodesPerformConstantSearch        ::Union{DataType,Union} = Nothing ,
    whichNodesPerformVariableSearch        ::Union{DataType,Union} = Nothing ,
    whichNodesPerformConstantVariableSearch::Union{DataType,Union} = Nothing ,
    performConstantExpressionSearch::Bool = false)

    if (obj.mse<=obj.goalMSE)
        return false;
    end;

    (nodes, heights, depths, _, numNodesEach) = iterateTree(obj.tree);

    # If the equations of the nodes were already calculated (the last search was unsuccessful), do not calculate them again
    if (obj.tree.equation==nothing)
        calculateEquations!(obj.tree, obj.initialEquation)
    end;

    if (obj.checkForErrors)
        checkEquations(obj.tree, obj.mse);
        # The equation must be correct: none of the semantic sets must be NaN (meaning that there are possible values in the domain)
        for node in nodes
            @assert(!isNaN(node.equation.S))
        end;
    end;

    subTreeInsert = nothing;
    indexNode = -1;
    maximumReduction = -Inf;
    maximumReduction = obj.minimumReductionMSE * obj.mse;

    #########################################################################################
    # Constant search
    if (whichNodesPerformConstantSearch!=Nothing)
        # Check if a constant can be inserted anywhere in the tree
        (maxReductionMSEConstantSearch, numNodeMaxReductionMSEConstantSearch, constantMaxReductionMSEConstantSearch) = constantSearch(nodes, obj.mse, numNodesEach, whichNodesPerformConstantSearch; checkForErrors=obj.checkForErrors);
        # if (maxReductionMSEConstantSearch>maximumReduction) && (maxReductionMSEConstantSearch>(obj.mse*obj.minimumDivReductionMSEConstantSearch))
        if (maxReductionMSEConstantSearch>maximumReduction)
            maximumReduction = maxReductionMSEConstantSearch;
            indexNode = numNodeMaxReductionMSEConstantSearch;
            # A subtree has been found
            subTreeInsert = Constant(constantMaxReductionMSEConstantSearch);
        end;

    end;

    #########################################################################################
    # Variable search
    if (whichNodesPerformVariableSearch!=Nothing)
        (maxReductionVariableSearch, indexNodeVariableSearch, indexNodePoolVariableSearch) = variableSearch(nodes, obj.nodePoolVariables, obj.mse, whichNodesPerformVariableSearch; checkForErrors=obj.checkForErrors);
        if (maxReductionVariableSearch>maximumReduction)
            maximumReduction = maxReductionVariableSearch;
            indexNode = indexNodeVariableSearch;
            # It is necessary to clone the variables; otherwise, they will be the same node in different parts of the tree, and will share the equation
            subTreeInsert = clone(obj.nodePoolVariables.nodes[indexNodePoolVariableSearch]);
        end;
    end;

    if (whichNodesPerformConstantVariableSearch!=Nothing) || (performConstantExpressionSearch)
        maximumHeightSubTree = (obj.maximumHeight+1).-depths;
        maximumNodesSubTree = (obj.maximumNodes - numNodesEach[1]) .+ numNodesEach;
    end;

    #########################################################################################
    # Constant-Variable search
    if (whichNodesPerformConstantVariableSearch!=Nothing)
        (maxReductionMSEConstantVariableSearch, numNodeMaxReductionMSEConstantVariableSearch, operationConstantVariableSearch, constantMaxReductionMSEConstantVariableSearch, variableMaxReductionMSEConstantVariableSearch) = constantVariableSearch(nodes, obj.nodePoolVariables, obj.mse, maximumHeightSubTree, maximumNodesSubTree, whichNodesPerformConstantVariableSearch; checkForErrors=obj.checkForErrors);
        if (maxReductionMSEConstantVariableSearch>maximumReduction)
            maximumReduction = maxReductionMSEConstantVariableSearch;
            indexNode = numNodeMaxReductionMSEConstantVariableSearch;
            # A subtree has been found
            subTreeInsert = operationConstantVariableSearch(Constant(constantMaxReductionMSEConstantVariableSearch), clone(variableMaxReductionMSEConstantVariableSearch));
        end;
    end;


    #########################################################################################
    # Constant-Expression search
    if (performConstantExpressionSearch)
        (maxReductionMSEConstantExpressionSearch, numNodeMaxReductionMSEConstantExpressionSearch, operationConstantExpressionSearch, constantMaxReductionMSEConstantExpressionSearch) = constantExpressionSearch(nodes, obj.mse, heights, maximumHeightSubTree, numNodesEach, maximumNodesSubTree; checkForErrors=obj.checkForErrors);
        if (maxReductionMSEConstantExpressionSearch>maximumReduction)
            maximumReduction = maxReductionMSEConstantExpressionSearch;
            indexNode = numNodeMaxReductionMSEConstantExpressionSearch;
            # A subtree has been found
            subTreeInsert = operationConstantExpressionSearch(Constant(constantMaxReductionMSEConstantExpressionSearch), nodes[indexNode]);
        end;
    end;

    # Do not clear equations if the search was unsuccessful: they can be used for the next search

    if (subTreeInsert==nothing)
        return false;
    end;

    # A subtree that can be inserted (whether it is constant, variable or a subtree) has been found

    clearEquations!(obj.tree);

    if (indexNode==1)
        obj.tree = subTreeInsert;
    else
        evaluateTree(subTreeInsert; checkForErrors=obj.checkForErrors);
        replaceSubtree!(obj.tree, nodes[indexNode], subTreeInsert);
    end;

    mse = MSE(obj.targets, evaluateTree(obj.tree; checkForErrors=obj.checkForErrors));

    # Check that the resulting MSE is lower than the previous MSE in obj.minimumReductionMSE
    # This is expected to be. However, precision errors may happen
    if (mse >= (obj.mse - (obj.minimumReductionMSE * obj.mse)))
        # Undo the change in the tree
        if (indexNode==1)
            obj.tree = nodes[1];
        else
            replaceSubtree!(obj.tree, subTreeInsert, nodes[indexNode]);
        end;
        mse = MSE(obj.targets, evaluateTree(obj.tree; checkForErrors=obj.checkForErrors));
        obj.checkForErrors && @assert(mse==obj.mse);
        return false;
    end;


    if (obj.checkForErrors)
        @assert(!isnan(mse));
        @assert(!isinf(mse));
        @assert(equal(mse+maximumReduction, obj.mse));
        if (!isinf(obj.maximumHeight))
            @assert (height(obj.tree)<=obj.maximumHeight);
        end;
        if (!isinf(obj.maximumNodes))
            @assert (numNodes(obj.tree)<=obj.maximumNodes);
        end;
        checkIntegrity(obj.nodePoolVariables);
    end;
    obj.mse = mse;

    return true;
end;








##################################################################################################
#
# Strategies
#    Return true if change has been done or false otherwise
#

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
    changeDone && OptimizeConstants!(obj);
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
    changeDone && OptimizeConstants!(obj);
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

Strategy1 = StrategyExhaustive;
Strategy2 = StrategyExhaustiveWithConstantOptimization;
Strategy3 = StrategySelectiveWithConstantOptimization;
Strategy4 = StrategySelective;




#################################################################################################################################################################
#
# The function Step!, that performs one iteration of the algorithm
#

function Step!(obj::DoME)
    obj.strategy(obj) && return true;
    # No change was done in the first tree: build different initial trees with a higher order
    if isa(obj.tree, Constant)
        constantValue = obj.tree.semantic;
        numNodesPreviousTree = 0;
        numNodesThisTree = numNodes(obj.tree);
        order = 0;
        while (numNodesPreviousTree!=numNodesThisTree)
            order = order + 1;
            obj.tree = buildInitialTree(order, convert(Array{Variable,1},obj.nodePoolVariables.nodes), constantValue, obj.maximumNodes);
            obj.mse = MSE(obj.targets, evaluateTree(obj.tree));
            obj.strategy(obj) && return true;
            numNodesPreviousTree = numNodesThisTree;
            numNodesThisTree = numNodes(obj.tree);
        end;
    end;
    return false;
end;





#################################################################################################################################################################
#
# The function is an interface that creates the DoME object and calls the function Step! on each iteration
#

function dome(inputs::Array{<:Real,2}, targets::Array{<:Real,1};
    # Each instance in inputs is in a row or in a column
    dataInRows          ::Bool           = true,
    # Hyperparameters of the algorithm
    minimumReductionMSE ::Real           = 1e-6,
    maximumNodes        ::Int64          = 50 ,
    strategy            ::Function       = StrategyExhaustive ,
    # Other hyperparameter that the user might find useful
    maximumHeight       ::Real           = Inf ,
    # Stopping criteria
    goalMSE             ::Real           = 0 ,
    maxIterations       ::Real           = Inf ,
    # Indices of the instances used for validation and test (validation has not been tested and may fail)
    validationIndices   ::Array{Int64,1} = Int64[],
    testIndices         ::Array{Int64,1} = Int64[],
    # Validation and test ratios (these parameters have not been tested and may fail)
    validationRatio     ::Real           = 0. ,
    testRatio           ::Real           = 0. ,
    # If you want to see the iterations on screen. This makes the execution slower
    showText            ::Bool           = false ,
    # This parameter was used only for development. If it is set to true, the execution becomes much slower
    checkForErrors      ::Bool           = false
    )


    if (dataInRows)
        @assert (length(targets)==size(inputs,1))
    else
        @assert (length(targets)==size(inputs,2))
    end;

    @assert(isempty(validationIndices) | (validationRatio==0.))
    @assert(isempty(testIndices) | (testRatio==0.))
    @assert(validationRatio+testRatio<1.);
    if ((!isempty(validationIndices)) & (!isempty(testIndices)))
        @assert(isempty(intersect(validationIndices,testIndices)));
    end;

    # Hold-out
    randomIndices = randperm(length(targets));
    if (isempty(validationIndices))
        numValidation = Int(round(length(targets)*validationRatio));
        validationIndices = randomIndices[1:numValidation];
        indices = sort(randomIndices[(numValidation+1):end]);
    else
        numValidation = length(validationIndices);
        randomIndices = setdiff(randomIndices, validationIndices);
    end;
    if (isempty(testIndices))
        numTest = Int(round(length(targets)*testRatio));
        testIndices = sort(randomIndices[1:numTest]);
        randomIndices = randomIndices[(numTest+1):end];
    else
        numTest = length(testIndices);
        randomIndices = setdiff(randomIndices, testIndices);
    end;
    numTraining = length(targets) - numValidation - numTest;
    trainingIndices = sort(randomIndices);
    @assert(numTraining == length(trainingIndices));
    @assert(isempty(intersect(trainingIndices,validationIndices)) & isempty(intersect(trainingIndices,testIndices)) & isempty(intersect(validationIndices,testIndices)));
    @assert(!isempty(trainingIndices));

    # Create variable values for validation and test sets
    TrainingVariableValues   = createVariableValues(dataInRows ? inputs[trainingIndices,:]   : inputs[:,trainingIndices];   dataInRows=dataInRows);
    ValidationVariableValues = createVariableValues(dataInRows ? inputs[validationIndices,:] : inputs[:,validationIndices]; dataInRows=dataInRows);
    TestVariableValues       = createVariableValues(dataInRows ? inputs[testIndices,:]       : inputs[:,testIndices];       dataInRows=dataInRows);

    sr = DoME(dataInRows ? inputs[trainingIndices,:] : inputs[:,trainingIndices], targets[trainingIndices]; dataInRows=dataInRows, maximumNodes=maximumNodes)

    sr.minimumReductionMSE = eltype(inputs)(minimumReductionMSE);
    sr.maximumHeight = maximumHeight;
    sr.checkForErrors = checkForErrors;
    sr.goalMSE = eltype(inputs)(goalMSE);
    sr.strategy = strategy;

    iteration = 0;

    BestTrainingMSE = Inf; BestValidationMSE = Inf; BestTestMSE = Inf; BestExpression = ""; BestTree =  nothing;

    function evaluateTreeIndices(tree, indices, variableValues)
        setVariableValues!(tree, variableValues);
        outputs = evaluateTree(tree);
        return MSE(targets[indices], outputs);
    end;

    function evaluateIteration()
        if (showText)
            println(" Iteration ", iteration);
            println(" Expression: ", string(sr.tree));
            println("   Nuber of nodes: ", numNodes(sr.tree), " - Height: ", height(sr.tree));
            println("   Current training:   MSE: ", sr.mse);
        end;

        # Perform validation and test of the expression, only if text is shown on screen (to make the algorithm faster)
        if (!isempty(testIndices) || !isempty(validationIndices)) && showText
            treeCopy = clone(sr.tree);

            # If there is validation set, the expression improves if it improves in the validation set
            if !isempty(validationIndices)
                mseValidation = evaluateTreeIndices(treeCopy, validationIndices, ValidationVariableValues);
                if (mseValidation<BestValidationMSE)
                    BestValidationMSE = mseValidation;
                    improves = true;
                else
                    improves = false;
                end;
            else
                # If there is no validation set, the expression is always improved
                improves = true;
            end;

            # If this expression is improved, calculate the result on the test set
            if improves
                BestTree = treeCopy;
                BestTrainingMSE = sr.mse;
                BestTestMSE = evaluateTreeIndices(treeCopy, testIndices, TestVariableValues);
            end;

            if showText
                # Print the overall results
                if !isempty(validationIndices)
                    println("   Overall training:   MSE: ", BestTrainingMSE);
                    println("           validation: MSE: ", BestValidationMSE);
                end;
                println("           test:       MSE: ", BestTestMSE);
                println("           Nuber of nodes: ", numNodes(BestTree), " - Height: ", height(BestTree));
                # clearEquations!(sr.tree); calculateEquations!(sr.tree, sr.initialEquation); println("   Max memory usage: ", Base.summarysize(sr)/(2^20), " MiB"); clearEquations!(sr.tree); # <- missing memory of the equations in each node
            end;
        end;
        if showText
            println("-------------------------------------------------------------------------------------")
        end;
    end;


    # Evaluate iteration 0
    evaluateIteration();

    stoppingCriteria = false;
    while !stoppingCriteria

        stoppingCriteria = !Step!(sr);

        if !stoppingCriteria

            iteration += 1;

            # Evaluate each iteration
            evaluateIteration();

            # Stopping criteria
            if (sr.mse<=goalMSE) || (iteration>=maxIterations)
                stoppingCriteria = true;
            end;

        end;

    end;

    if BestTree==nothing
        BestTree = sr.tree;
        BestTrainingMSE = sr.mse;
        BestValidationMSE = evaluateTreeIndices(sr.tree, validationIndices, ValidationVariableValues);
        BestTestMSE = evaluateTreeIndices(sr.tree, testIndices, TestVariableValues);
    end;

    if showText
        println("Best expression found: ", string(BestTree));
    end;

    return (BestTrainingMSE, BestValidationMSE, BestTestMSE, BestTree);

end;

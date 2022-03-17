
using Random

include("Equation.jl")

macro def(name, definition)
  return quote
      macro $(esc(name))()
          esc($(Expr(:quote, definition)))
      end
  end
end

@def addNodeFields begin
    equation::Union{Nothing,NodeEquation}
end



abstract type Tree end

abstract type Terminal <: Tree end

mutable struct Variable <: Terminal
    @addNodeFields
    variableNumber::UInt
    semantic::Array{<:Real,1}
    Variable(num::UInt; semantic=[]) = new(nothing, num, semantic)
end

mutable struct Constant <: Terminal
    @addNodeFields
    semantic::Real
    Constant(value::Real) = new(nothing, value)
end

mutable struct RandomConstant <: Terminal
    @addNodeFields
    lowerLimit::Real
    upperLimit::Real
    RandomConstant(lowerLimit::Real, upperLimit::Real) = new(nothing, lowerLimit, upperLimit)
end


@def addNonTerminalFields begin
    name::String
    semantic::Union{Nothing,Semantic}
    evalFunction::Function
end

abstract type NonTerminal <: Tree end

mutable struct BinaryNode <: NonTerminal
    @addNodeFields
    @addNonTerminalFields
    child1::Tree
    child2::Tree
    functionEquationChild1::Function
    functionEquationChild2::Function
    function BinaryNode(name::String, evalFunction::Function, child1::Tree, child2::Tree, functionEquationChild1::Function, functionEquationChild2::Function; semantic=nothing)
        return new(nothing, name, semantic, evalFunction, child1, child2, functionEquationChild1, functionEquationChild2);
    end;
end

mutable struct NonBinaryNode <: NonTerminal
    @addNodeFields
    @addNonTerminalFields
    children::Array{Tree,1}
    functionEquationChildren::Array{Function,1}
    function NonBinaryNode(name::String, evalFunction::Function, functionEquationChildren::Array{Function,1}; semantic=nothing, children::Array{Tree,1}=[])
        return new(nothing, name, semantic, evalFunction, children, functionEquationChildren);
    end;
    function NonBinaryNode(type::Int, name::String, evalFunction::Function, functionEquationChildren::Array{Function,1}, children::Array{Tree,1}; semantic=nothing)
        return new(nothing, name, semantic, evalFunction, children, functionEquationChildren);
    end;
end


AddNode(child1::Tree, child2::Tree) = BinaryNode("+", addSemantics, child1, child2, equationAddChild1, equationAddChild2);
SubNode(child1::Tree, child2::Tree) = BinaryNode("-", subSemantics, child1, child2, equationSubChild1, equationSubChild2);
MulNode(child1::Tree, child2::Tree) = BinaryNode("*", mulSemantics, child1, child2, equationMulChild1, equationMulChild2);
DivNode(child1::Tree, child2::Tree) = BinaryNode("/", divSemantics, child1, child2, equationDivChild1, equationDivChild2);



clearEvaluationValues!(node::Constant) = nothing;
clearEvaluationValues!(node::RandomConstant) = error("");
clearEvaluationValues!(node::Variable) = nothing;
function clearEvaluationValues!(node::BinaryNode)
    node.semantic = nothing;
    clearEvaluationValues!(node.child1);
    clearEvaluationValues!(node.child2);
    return nothing;
end
function clearEvaluationValues!(node::NonBinaryNode)
    node.semantic = nothing;
    for childNode in node.children
        clearEvaluationValues!(childNode);
    end;
    return nothing;
end


# These 3 functions clear the equations except in the given path
clearEquations!(node::Terminal) = ( node.equation = nothing; )
function clearEquations!(node::BinaryNode; exceptInPath=nothing)
    if (exceptInPath==nothing)
        node.equation = nothing;
        clearEquations!(node.child1);
        clearEquations!(node.child2);
    else
        if (isempty(exceptInPath))
            clearEquations!(node.child1);
            clearEquations!(node.child2);
        elseif (exceptInPath[1]==1)
            clearEquations!(node.child1, exceptInPath[2:end]);
            clearEquations!(node.child2);
        else
            clearEquations!(node.child1);
            clearEquations!(node.child2, exceptInPath[2:end]);
        end;
    end;
end
function clearEquations!(node::NonBinaryNode; exceptInPath=nothing)
    if (exceptInPath==nothing)
        node.equation = nothing;
        for child in node.children
            clearEquations!(child);
        end;
    else
        if (isempty(exceptInPath))
            for child in node.children
                clearEquations!(child);
            end;
        else
            nextStep = exceptInPath[1];
            for numChild in 1:length(node.children)
                if (numChild==nextStep)
                    clearEvaluationValues!(node.children[numChild], exceptInPath[2:end]);
                else
                    clearEvaluationValues!(node.children[numChild]);
                end;
            end;
        end;
    end;
end

function clearEquations!(node::BinaryNode; exceptInPaths=nothing)
    if (exceptInPaths==nothing)
        node.equation = nothing;
        clearEquations!(node.child1);
        clearEquations!(node.child2);
    else
        if (isempty(exceptInPaths))
            clearEquations!(node.child1);
            clearEquations!(node.child2);
        else
            firstSteps = [path[1] for path in exceptInPaths];
            has1 = anyequal(1,firstSteps);
            has2 = anyequal(2,firstSteps);
            if (has1) && (has2)
                followingSteps = exceptInPaths[firstSteps.==1];
                followingSteps = [path[2:end] for path in followingSteps];
                followingSteps = followingSteps[!isempty.(followingSteps)];
                clearEquations!(node.child1, exceptInPaths=followingSteps);
                followingSteps = exceptInPaths[firstSteps.==2];
                followingSteps = [path[2:end] for path in followingSteps];
                followingSteps = followingSteps[!isempty.(followingSteps)];
                clearEquations!(node.child1, exceptInPaths=followingSteps);
            elseif (has1)
                followingSteps = [path[2:end] for path in exceptInPaths];
                followingSteps = followingSteps[!isempty.(followingSteps)];
                clearEquations!(node.child1, exceptInPaths=followingSteps);
                clearEquations!(node.child2);
            elseif (has2)
                followingSteps = [path[2:end] for path in exceptInPaths];
                followingSteps = followingSteps[!isempty.(followingSteps)];
                clearEquations!(node.child1);
                clearEquations!(node.child2, exceptInPaths=followingSteps);
            else
                error("Don't know which path to go")
            end;
        end;
    end;
end




import Base.string
string(node::Constant) = node.semantic>=0 ? string(node.semantic) : string("(",node.semantic,")");
string(node::RandomConstant) = error("");
string(node::Variable) = string("X", node.variableNumber);
string(node::BinaryNode) = string("(",string(node.child1),node.name,string(node.child2),")")
function string(node::NonBinaryNode)
    if (length(node.children)==2)
        text = string("(",string(node.children[1]),node.name,string(node.children[2]),")")
    else
        text = string(node.name, "(");
        for numChildren = 1:length(node.children)
            text = string(text, string(node.children[numChildren]));
            if (numChildren!=length(node.children))
                text = string(text, ",");
            end;
        end;
        text = string(text, ")");
    end;
    return text;
end;


vectorString(node::Constant; dataInRows=true) = node.semantic>=0 ? string(node.semantic) : string("(",node.semantic,")");
vectorString(node::RandomConstant; dataInRows=true) = error("");
vectorString(node::Variable; dataInRows=true) = dataInRows ? string("X[:,", node.variableNumber,"]") : string("X[", node.variableNumber,",:]");
vectorString(node::BinaryNode; dataInRows=true) = string("(", vectorString(node.child1; dataInRows=dataInRows), " .", node.name, " ", vectorString(node.child2; dataInRows=dataInRows), ")")
function vectorString(node::NonBinaryNode; dataInRows=true)
    if (length(node.children)==2)
        text = string("(", vectorString(node.children[1]; dataInRows=dataInRows), " .", node.name, " ", vectorString(node.children[2]; dataInRows=dataInRows), ")")
    else
        text = string(node.name, ".(");
        for numChildren = 1:length(node.children)
            text = string(text, vectorString(node.children[numChildren]; dataInRows=dataInRows));
            if (numChildren!=length(node.children))
                text = string(text, ",");
            end;
        end;
        text = string(text, ")");
    end;
    return text;
end;


latexString(node::Constant) = string(node);
latexString(node::RandomConstant) = string(node);
latexString(node::Variable) = string("X_{", node.variableNumber,"}");
latexString(node::BinaryNode) = (node.name=="/") ? string("\\frac{",latexString(node.child1),"}{",latexString(node.child2),"}") : string("\\left(",latexString(node.child1),(node.name=="*") ? "\\cdot " : node.name,latexString(node.child2),"\\right)")
latexString(node::NonBinaryNode) = string(node);
function latexString(node::NonBinaryNode)
    if (length(node.children)==2)
        text = string("\\left(",string(node.children[1]),node.name,string(node.children[2]),"\\right)")
    else
        text = string(node.name, "\\left(");
        for numChildren = 1:length(node.children)
            text = string(text, string(node.children[numChildren]));
            if (numChildren!=length(node.children))
                text = string(text, ',');
            end;
        end;
        text = string(text, "\\right)");
    end;
    return text;
end;



evaluateTree(tree::Terminal; checkForErrors=false) = tree.semantic;
evaluateTree(tree::RandomConstant; checkForErrors=false) = error("");
function evaluateTree(tree::BinaryNode; checkForErrors=false)
    if (tree.semantic==nothing)
        tree.semantic = tree.evalFunction( evaluateTree(tree.child1;checkForErrors=checkForErrors), evaluateTree(tree.child2;checkForErrors=checkForErrors), nothing );
    end;
    return tree.semantic;
end;
function evaluateTree(tree::NonBinaryNode; checkForErrors=false)
    if (tree.semantic==nothing)
        evaluationChildren = [evaluateTree(child; checkForErrors=checkForErrors) for child in tree.children];
        tree.semantic = tree.evalFunction( evaluationChildren... );
    end;
    if (checkForErrors)
        @assert(!any(isinf.(tree.semantic)));
        @assert(!any(isnan.(tree.semantic)));
    end;
    return tree.semantic;
end;

reevaluatePath(tree::Terminal, path::Array{Int64,1}, indexPath::Int64; checkForErrors=false) = tree.semantic;
reevaluatePath(tree::RandomConstant, path::Array{Int64,1}, indexPath::Int64; checkForErrors=false) = error("");
function reevaluatePath(tree::BinaryNode, path::Array{Int64,1}, indexPath::Int64=1; checkForErrors=false)
    if (indexPath>length(path))
        semanticChild1 =   evaluateTree(tree.child1;            checkForErrors=checkForErrors);
        semanticChild2 =   evaluateTree(tree.child2;            checkForErrors=checkForErrors);
    elseif (path[indexPath]==1)
        semanticChild1 = reevaluatePath(tree.child1, path, indexPath+1; checkForErrors=checkForErrors);
        semanticChild2 =   evaluateTree(tree.child2;            checkForErrors=checkForErrors);
    elseif (path[indexPath]==2)
        semanticChild1 =   evaluateTree(tree.child1;            checkForErrors=checkForErrors);
        semanticChild2 = reevaluatePath(tree.child2, path, indexPath+1; checkForErrors=checkForErrors);
    else
        error("Don't know which path to go")
    end;
    tree.semantic = tree.evalFunction( semanticChild1, semanticChild2, tree.semantic );
    if (checkForErrors)
        @assert(!any(isinf.(tree.semantic)));
        @assert(!any(isnan.(tree.semantic)));
    end;
    return tree.semantic;
end;




numNodes(node::Terminal) = UInt(1);
numNodes(node::BinaryNode) = UInt(1 + numNodes(node.child1) + numNodes(node.child2));
numNodes(node::NonBinaryNode) = UInt(1 + sum([numNodes(child) for child in node.children]));

height(node::Terminal) = UInt(1);
height(node::BinaryNode) = UInt(1) + max(height(node.child1), height(node.child2));
height(node::NonBinaryNode) = UInt(1) + maximum([height(child) for child in node.children]);

clone(node::Constant) = Constant(node.semantic);
clone(node::RandomConstant) = Constant(rand()*(node.upperLimit-node.lowerLimit) + node.lowerLimit, node.type);
clone(node::Variable) = Variable(node.variableNumber; semantic=((node.semantic==nothing) ? node.semantic : copy(node.semantic)));
clone(node::BinaryNode) = BinaryNode(node.name, node.evalFunction, clone(node.child1), clone(node.child2), node.functionEquationChild1, node.functionEquationChild2; semantic=((node.semantic==nothing) ? node.semantic : copy(node.semantic)));
clone(node::NonBinaryNode) = NonBinaryNode(node.name, node.evalFunction, node.equationsFunction; semantic=((node.semantic==nothing) ? node.semantic : copy(node.semantic)), children=convert(Array{Tree,1},[clone(child) for child in node.children]));

iterateTree(tree::Terminal) = (convert(Array{Tree,1},[tree]), [0x0000000000000001], [0x0000000000000001], [0x0000000000000000], [0x0000000000000001]);
function iterateTree(tree::BinaryNode)
    (nodesChild1, heightsChild1, depthsChild1, numChildrenChild1, numNodesChild1) = iterateTree(tree.child1);
    (nodesChild2, heightsChild2, depthsChild2, numChildrenChild2, numNodesChild2) = iterateTree(tree.child2);
    nodes = [nodesChild1; nodesChild2]; pushfirst!(nodes, tree);
    heights = [heightsChild1; heightsChild2]; pushfirst!(heights, 0x0000000000000000);
    depths = [depthsChild1; depthsChild2]; pushfirst!(depths, 0x0000000000000001);
    numChildren = [numChildrenChild1; numChildrenChild2]; pushfirst!(numChildren, 0x0000000000000001);
    numNodes = [numNodesChild1; numNodesChild2]; pushfirst!(numNodes, 0x0000000000000000);
    heights[1] = maximum(heights[2:end])+1;
    numNodes[1] = length(numNodes);
    depths[2:end] .+= 1;
    return (nodes, heights, depths, numChildren, numNodes);
end
function iterateTree(tree::NonBinaryNode)
    nodes = [tree];
    heights = [0x0000000000000000];
    depths = [0x0000000000000001];
    numChildren = [length(tree.children)];
    numNodes = [0x0000000000000000];
    for child in tree.children
        (nodesThisChild, tiposThisChild, heightsThisChild, depthsThisChild, numChildrenThisChild, numNodesThisChild) = iterateTree(child);
        nodes = [nodes; nodesThisChild];
        # hcat(tipos, tiposEsteHijo);
        heights = [heights; heightsThisChild];
        depths = [depths; depthsThisChild.+1];
        numChildren = [numChildren; numChildrenThisChild];
        numNodes = [numNodes; numNodesThisChild];
    end;
    heights[1] = maximum(heights[2:end])+1;
    numNodes[1] = length(numNodes);
    depths[2:end] .+= 1;
    return (nodes, heights, depths, numChildren, numNodes);
end



calculateEquations!(tree::Terminal, equation::NodeEquation; path=nothing, indexPath=1, checkForErrors=false) = ( tree.equation = equation; );
function calculateEquations!(tree::BinaryNode, equation::NodeEquation; path=nothing, indexPath=1, checkForErrors=false)
    if (checkForErrors)
        @assert (tree.semantic != nothing);
        checkEquation(equation);
    end;
    tree.equation = equation;

    if (path==nothing)
        calculateEquations!(tree.child1, tree.functionEquationChild1( evaluateTree(tree.child2; checkForErrors=checkForErrors), equation); checkForErrors=checkForErrors );
        calculateEquations!(tree.child2, tree.functionEquationChild2( evaluateTree(tree.child1; checkForErrors=checkForErrors), equation); checkForErrors=checkForErrors );
    else
        (indexPath>length(path)) && return
        if (path[indexPath]==1)
            calculateEquations!(tree.child1, tree.functionEquationChild1( evaluateTree(tree.child2; checkForErrors=checkForErrors), equation); path=path, indexPath=indexPath+1, checkForErrors=checkForErrors );
        elseif (path[indexPath]==2)
            calculateEquations!(tree.child2, tree.functionEquationChild2( evaluateTree(tree.child1; checkForErrors=checkForErrors), equation); path=path, indexPath=indexPath+1, checkForErrors=checkForErrors );
        else
            error("Don't know which path to go");
        end;
    end;
end;






setVariableValues!(tree::Terminal, variableValues::Array{<:Array{<:Real,1}}) = ()
setVariableValues!(tree::Variable, variableValues::Array{<:Array{<:Real,1}}) = (tree.semantic = variableValues[tree.variableNumber]; )
setVariableValues!(tree::BinaryNode, variableValues::Array{<:Array{<:Real,1}}) = (setVariableValues!(tree.child1, variableValues); setVariableValues!(tree.child2, variableValues); tree.semantic = nothing; )
setVariableValues!(tree::NonBinaryNode, variableValues::Array{<:Array{<:Real,1}}) = (for child in tree.children setVariableValues!(child, variableValues); end; tree.semantic = nothing; )




findNode(tree::Terminal, node::Tree) = (tree==node) ? Array{Int,1}([]) : nothing
function findNode(tree::BinaryNode, node::Tree)
    (tree==node) && return Array{Int,1}([]);
    path = findNode(tree.child1, node);
    if (path!=nothing)
        pushfirst!(path, 1);
        return path;
    end;
    path = findNode(tree.child2, node);
    if (path!=nothing)
        pushfirst!(path, 2);
        return path;
    end;
    return nothing;
end;
function findNode(tree::NonBinaryNode, node::Tree)
    (tree==node) && return Array{Int,1}([]);
	for numChild in 1:length(tree.children)
        path = findNode(tree.children[numChild], node);
        if (path!=nothing)
            pushfirst!(path, numChild);
            return path;
        end;
    end;
    return nothing;
end;



replaceSubtree!(tree::Terminal, oldSubTree::Tree, newSubTree::Tree) = false;
function replaceSubtree!(tree::BinaryNode, oldSubTree::Tree, newSubTree::Tree)
    replaced = false;
    if (tree.child1==oldSubTree)
        tree.child1 = newSubTree;
        replaced = true;
    elseif (tree.child2==oldSubTree)
        tree.child2 = newSubTree;
        replaced = true;
    else
        replaced = replaceSubtree!(tree.child1, oldSubTree, newSubTree);
        if !replaced
            replaced = replaceSubtree!(tree.child2, oldSubTree, newSubTree);
        end;
    end;
    if (replaced)
        tree.semantic = reevaluatePath(tree, Int64[]);
        return true;
    end;
    return false;
end
function replaceSubtree!(tree::NonBinaryNode, oldSubTree::Tree, newSubTree::Tree)
	for numChild in 1:length(tree.children)
        child = tree.children[numChild]
        if (child==oldSubTree)
            tree.children[numChild] = newSubTree;
            tree.semantic = reevaluatePath(tree, Int64[]);
            return true;
        else
            if (replaceSubtree!(child, oldSubTree, newSubTree))
                tree.semantic = reevaluatePath(tree, Int64[]);
                return true;
            end;
        end;
    end;
    return false;
end



function checkEquations(tree::Tree, mse::Real)
    @assert(!any(isnan.(tree.semantic)));
    @assert(!any(isinf.(tree.semantic)));
    if (tree.equation!=nothing)
        checkEquation(tree.equation);
        @assert(equal(mse, calculateMSEFromEquation(tree.semantic, tree.equation; checkForErrors=true)));
        @assert(equal(mse, calculateMSEFromEquationEfficient(tree.semantic, tree.equation; checkForErrors=true)));
        @assert(equal(mse, calculateMSEFromEquationMoreEfficient(tree.semantic, tree.equation; checkForErrors=true)));
    end;
    if (isa(tree,BinaryNode))
        checkEquations(tree.child1, mse)
        checkEquations(tree.child2, mse)
    elseif (isa(tree,NonBinaryNode))
        for numChild in 1:length(tree.children)
            checkEquations(tree.children[numChild], mse)
        end;
    end;
end;





function buildInitialTree(order::Int64, variables::Array{Variable,1}, constant::Real, maximumNodes::Number)

    function buildTreeOrder(order::Int64, variables::AbstractArray{Variable,1})

        function buildVariableCombinationTreeList(thisOrder::Int64, variablesThisOrder::AbstractArray{Variable,1})
            (thisOrder==0) && return [Constant(eltype(constant)(0.))];
            treeListThisOrder = [];
            for numVariable in 1:length(variablesThisOrder)
                subTrees = buildVariableCombinationTreeList(thisOrder-1, view(variablesThisOrder, numVariable:length(variablesThisOrder)));
                for subTree in subTrees
                    push!(treeListThisOrder, MulNode(subTree, clone(variablesThisOrder[numVariable])) );
                end;
            end;
            return treeListThisOrder;
        end;

        treeList = buildVariableCombinationTreeList(order, variables);
        tree = treeList[1];
        for subTree in view(treeList, 2:length(treeList))
            tree = AddNode( tree, subTree );
        end;
        return tree;
    end;

    outputTree = Constant(constant);


    numNodesCurrent = 1;

    for i in 1:order
        treeOrder = buildTreeOrder(i, variables);
        numNodesTreeOrder = numNodes(treeOrder);
        if (numNodesCurrent+numNodesTreeOrder+1<=maximumNodes)
            outputTree = AddNode(outputTree, treeOrder);
            numNodesCurrent += numNodesTreeOrder + 1;
        else
            # The rest of the trees of higher order are going to have even more nodes,
            #  so exit this function
            return outputTree;
        end;
    end;

    return outputTree;
end;

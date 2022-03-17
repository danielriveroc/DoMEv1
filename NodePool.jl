
include("Tree.jl")

mutable struct NodePool
    population::Array{Tree,1}
    nodes::Array{Tree,1}
    numNodes::Array{UInt,1}
    heights::Array{UInt,1}
    semantics::Array{Array{<:Real,1}}
    hasZerosInSemantics::Array{Bool,1}
    toleranceValue::Real
    NodePool() = new([],[],[],[],[],[], 0.);
end


function checkIntegrity(sol::NodePool)
    @assert (length(sol.nodes)==length(sol.numNodes));
    @assert (length(sol.nodes)==length(sol.heights));
    @assert (length(sol.nodes)==length(sol.semantics));
    for numNode = 1:length(sol.nodes)
        @assert (numNodes(sol.nodes[numNode])==sol.numNodes[numNode]);
        @assert (height(sol.nodes[numNode])==sol.heights[numNode]);
    end;
end;

function removeNodesFromPool!(sol::NodePool, indexNodes::Array{Int,1})
    numberOfNodes = length(sol.nodes);
    indexRemaining = setdiff(1:numberOfNodes, indexNodes);

    sol.nodes = sol.nodes[indexRemaining];
    sol.numNodes = sol.numNodes[indexRemaining];
    sol.heights = sol.heights[indexRemaining];
    sol.semantics = sol.semantics[indexRemaining];
    return indexRemaining;
end;

function addNodeToPool!(sol::NodePool, node::Tree, numNodes::UInt, height::UInt)
    semantic = node.semantic;
    @assert (semantic!=nothing);
    if (sol.toleranceValue!=0.)
        add = all( sqrt.( (x->sum(x.^2)).([semantic] .- sol.semantics)) .> sol.toleranceValue );
    else
        add = all( [semantic] .!= sol.semantics );
    end;

    if (add)
        push!(sol.nodes, node);
        push!(sol.numNodes, numNodes);
        push!(sol.heights, height);
        push!(sol.semantics, semantic);
        # push!(sol.hasZerosInSemantics, any(semantic.==0));
        push!(sol.hasZerosInSemantics, any0(semantic));
        return true;
    else
        return false;
    end;
end;


include("Tree.jl");

x3 = Variable(UInt(3), semantic=[1.,4.,-2.]);

tree = DivNode(
    Constant(2.) ,
    AddNode(
        DivNode(
            Constant(3.) ,
            SubNode(
                Constant(1.) ,
                DivNode(
                    Constant(2.),
                    x3
                )
            )
        ) ,
        Constant(1.)
    )
);

println("Figure 1: ", string(tree));

targets = [5.,4.,1.];
evaluateTree(tree);

initialEquation = NodeEquation(1., targets, 0., -1.);

calculateEquations!(tree,initialEquation);

println("Table 1:");

nodes, = iterateTree(tree);
for numNode in 1:length(nodes)
    println("-----------------------------------------------------------------------------------------")
    if (isa(nodes[numNode], NonTerminal))
        println("Node ", numNode, ": ", nodes[numNode].name);
    else
        println("Node ", numNode, ": ", string(nodes[numNode]));
    end;
    showEquation(nodes[numNode].equation, nodes[numNode].semantic);
end;




println("-----------------------------------------------------------------------------------------")
println("-----------------------------------------------------------------------------------------")
println("-----------------------------------------------------------------------------------------")

x1 = Variable(UInt(1), semantic=[3.,2.,0.]);
x2 = Variable(UInt(2), semantic=[0.,1.,0.]);
x3 = Variable(UInt(3), semantic=[1.,0.,3.]);
# x4 = Variable(UInt(4), semantic=[2.,0.,4.]);

tree = DivNode(
    Constant(2.) ,
    AddNode(
        DivNode(
            x1,
            Constant(2.)
        ) ,
        MulNode(
            x2,
            SubNode(
                MulNode(
                    x3,
                    # DivNode(
                        Constant(2.) ,
                        # x4
                    # ) ,
                ) ,
                Constant(1.)
            )
        )
    )
);
println("Figure 3: ", string(tree));

evaluateTree(tree);

targets = [5.,4.,1.];
initialEquation = NodeEquation(1., targets, 0., -1.);

calculateEquations!(tree,initialEquation);

println("Table 3: ");

nodes, = iterateTree(tree);
for numNode in 1:length(nodes)
    println("-----------------------------------------------------------------------------------------")
    if (isa(nodes[numNode], NonTerminal))
        println("Node ", numNode, ": ", nodes[numNode].name);
    else
        println("Node ", numNode, ": ", string(nodes[numNode]));
    end;
    showEquation(nodes[numNode].equation, nodes[numNode].semantic);
end;


using Statistics

Semantic = Union{<:Real,Array{<:Real,1}}
SemanticSet = Union{<:Real,Array{Semantic,1}}

@inline valueSemantic(semantic::Real, index::Int64) = semantic;
@inline valueSemantic(semantic::Array{<:Real,1}, index::Int64) = @inbounds semantic[index];

# Use "mutable" here so this struct is allocated in the heap, and the garbage collector can free the memory stored on a,b,c,d,S
mutable struct NodeEquation
    a::Semantic
    b::Semantic
    c::Semantic
    d::Semantic
    S::SemanticSet
    NodeEquation(a,b,c,d) = new(a,b,c,d,Semantic[])
    NodeEquation(a,b,c,d,S) = new(a,b,c,d,S)
end

#############################################################################################################
#
# Useful functions for evaluating the semantics
#

@inline allequal(x::Real) = true
@inline function allequal(x::Array{<:Real,1})
    length(x) < 2 && return true
    v = x[1]
    @inbounds for i=2:length(x)
        x[i] == v || return false
    end
    return true
end

@inline any0(x::Real) = (x==0)
@inline function any0(x::Array{<:Real,1})
    @inbounds for i=1:length(x)
        x[i] == 0 && return true
    end
    return false
end

# If vectors x and y have any zero in the same position
@inline any0(x::Real,y::Real) = (x==0) && (y==0)
@inline any0(x::Real,y::Array{<:Real,1}) = x==0 ? any0(y) : false
@inline any0(x::Array{<:Real,1},y::Real) = any0(y,x)
@inline function any0(x::Array{<:Real,1},y::Array{<:Real,1})
    @inbounds for i=1:length(x)
        (x[i]==0) && (y[i]==0) && return true
    end
    return false
end

@inline anyInf(x::Real) = isinf(x)
@inline function anyInf(x::Array{<:Real,1})
    @inbounds for i=1:length(x)
        isinf(x[i]) && return true
    end
    return false
end

@inline anyNaN(x::Real) = isnan(x)
@inline function anyNaN(x::Array{<:Real,1})
    @inbounds for i=1:length(x)
        isnan(x[i]) && return true
    end
    return false
end

@inline all0(x::Real) = (x==0)
@inline function all0(x::Array{<:Real,1})
    @inbounds for i=1:length(x)
        x[i]!=0 && return false
    end
    return true
end

@inline all1(x::Real) = (x==1)
@inline function all1(x::Array{<:Real,1})
    @inbounds for i=1:length(x)
        x[i]!=1 && return false
    end
    return true
end

@inline anyequal(x::Real, y::Real) = (x==y)
@inline function anyequal(x::Real, y::Array{<:Real,1})
    @inbounds for i=1:length(y)
        y[i]==x && return true
    end
    return false
end
@inline anyequal(x::Int64, y::Int64) = (x==y)
@inline function anyequal(x::Int64, y::Array{Int64,1})
    @inbounds for i=1:length(y)
        y[i]==x && return true
    end
    return false
end
@inline anyequal(x::Array{<:Real,1}, y::Real) = anyequal(y,x)
@inline function anyequal(x::Array{<:Real,1}, y::Array{<:Real,1})
    @inbounds for i=1:length(x)
        x[i]==y[i] && return true
    end
    return false
end



# Only for checking errors
equal(x,y; tolerance=NaN) = (x_y = abs.(x.-y); indNot0 = isa(x_y,Number) ? ((x_y!=0) ? 1 : []) : findall(x_y.!=0); if isnan(tolerance) tolerance = eltype(x)==Float32 ? 1e-3 : 1e-5; end; isempty(indNot0) ? true : all( ((x_y.<tolerance) .| (x_y.<abs.(x.*tolerance)) .| (x_y.<abs.(y.*tolerance)))[indNot0]) );

@inline cleanVector(v::Real) = v;
@inline cleanVector(v::Array{<:Real,1}) = (allequal(v)) ? v[1] : v


#############################################################################################################
#
# Functions for the operators
#

# These equations make the same computations, but allow to specify a buffer to write the results
# This avoids memory allocations and thus are faster
@inline addSemantics(x, y, buffer::Array{<:Real,1}) = (isa(x,Number) && isa(y,Number)) ? x+y : ( buffer .= x; buffer .+= y; cleanVector(buffer); )
@inline subSemantics(x, y, buffer::Array{<:Real,1}) = (isa(x,Number) && isa(y,Number)) ? x-y : ( buffer .= x; buffer .-= y; cleanVector(buffer); )
@inline mulSemantics(x, y, buffer::Array{<:Real,1}) = (isa(x,Number) && isa(y,Number)) ? x*y : ( buffer .= x; buffer .*= y; cleanVector(buffer); )
@inline divSemantics(x, y, buffer::Array{<:Real,1}) = (isa(x,Number) && isa(y,Number)) ? x/y : ( buffer .= x; buffer ./= y; cleanVector(buffer); )
@inline addSemantics(x, y, buffer::Nothing) = x.+y;
@inline subSemantics(x, y, buffer::Nothing) = x.-y;
@inline mulSemantics(x, y, buffer::Nothing) = x.*y;
@inline divSemantics(x, y, buffer::Nothing) = x./y;
@inline addSemantics(x, y, buffer::Real) = x.+y;
@inline subSemantics(x, y, buffer::Real) = x.-y;
@inline mulSemantics(x, y, buffer::Real) = x.*y;
@inline divSemantics(x, y, buffer::Real) = x./y;
@inline addSemantics(x, y) = x.+y;
@inline subSemantics(x, y) = x.-y;
@inline mulSemantics(x, y) = x.*y;
@inline divSemantics(x, y) = x./y;

#############################################################################################################
#
# Functions for calculating the equation coefficients a,b,c,d, and S
#

isNaN(S::SemanticSet) = isa(S,Number) && isnan(S)
semanticSetOpSemantic(S::SemanticSet, f::Function, semantic::Semantic) = (isempty(S) || isNaN(S)) ? S : Semantic[cleanVector(f.(s,semantic)) for s in S]
semanticOpSemanticSet(semantic::Semantic, f::Function, S::SemanticSet) = (isempty(S) || isNaN(S)) ? S : Semantic[cleanVector(f.(semantic,s)) for s in S]
function semanticSetDivSemantic(S::SemanticSet, semantic::Semantic)
    (isempty(S) || isNaN(S)) && return S;
    for s in S
        any0(s,semantic) && return NaN;
    end;
    return Semantic[cleanVector(s./semantic) for s in S];
end;
function semanticDivSemanticSet(semantic::Semantic, S::SemanticSet)
    isempty(S) && return Semantic[0.];
    isNaN(S) && return S;
    for s in S
        any0(s,semantic) && return NaN;
    end;
    newSemanticSet = Semantic[cleanVector(semantic./s) for s in S];
    for s in newSemanticSet
        s==0 && return newSemanticSet
    end;
    pushfirst!(newSemanticSet,0.);
    return newSemanticSet
end;


@inline equationAddChild1(y, a, b, c, d, S) = NodeEquation(                               a , all1(a) ? cleanVector(b.-y) : cleanVector(b.-(a.*y)) ,                               c , all0(c) ?  d : cleanVector(d.-(c.*y)), semanticSetOpSemantic(S,-,y) );
@inline equationAddChild2(x, a, b, c, d, S) = NodeEquation(                               a , all1(a) ? cleanVector(b.-x) : cleanVector(b.-(a.*x)) ,                               c , all0(c) ?  d : cleanVector(d.-(c.*x)), semanticSetOpSemantic(S,-,x) );
@inline equationSubChild1(y, a, b, c, d, S) = NodeEquation(                               a , all1(a) ? cleanVector(b.+y) : cleanVector(b.+(a.*y)) ,                               c , all0(c) ?  d : cleanVector(d.+(c.*y)), semanticSetOpSemantic(S,+,y) );
@inline equationSubChild2(x, a, b, c, d, S) = NodeEquation(                               a , all1(a) ? cleanVector(x.-b) : cleanVector((a.*x).-b) ,                               c , all0(c) ? -d : cleanVector((c.*x).-d), semanticOpSemanticSet(x,-,S) );
@inline equationMulChild1(y, a, b, c, d, S) = NodeEquation( all1(a) ? y : cleanVector(a.*y) ,                                                    b , all0(c) ? c : cleanVector(c.*y) ,                                     d, semanticSetDivSemantic(S,y) );
@inline equationMulChild2(x, a, b, c, d, S) = NodeEquation( all1(a) ? x : cleanVector(a.*x) ,                                                    b , all0(c) ? c : cleanVector(c.*x) ,                                     d, semanticSetDivSemantic(S,x) );
@inline equationDivChild1(y, a, b, c, d, S) = NodeEquation(                               a ,                                    cleanVector(b.*y) ,                               c , all0(d) ?  d :      cleanVector(d.*y), semanticSetOpSemantic(S,*,y) );
@inline equationDivChild2(x, a, b, c, d, S) = NodeEquation(                               b , all1(a) ?                 x :      cleanVector(a.*x) ,                               d , all0(c) ?  c :      cleanVector(c.*x), semanticDivSemanticSet(x,S) );

@inline equationAddChild1(y, eq::NodeEquation) = equationAddChild1(y, eq.a, eq.b, eq.c, eq.d, eq.S)
@inline equationAddChild2(x, eq::NodeEquation) = equationAddChild2(x, eq.a, eq.b, eq.c, eq.d, eq.S)
@inline equationSubChild1(y, eq::NodeEquation) = equationSubChild1(y, eq.a, eq.b, eq.c, eq.d, eq.S)
@inline equationSubChild2(x, eq::NodeEquation) = equationSubChild2(x, eq.a, eq.b, eq.c, eq.d, eq.S)
@inline equationMulChild1(y, eq::NodeEquation) = equationMulChild1(y, eq.a, eq.b, eq.c, eq.d, eq.S)
@inline equationMulChild2(x, eq::NodeEquation) = equationMulChild2(x, eq.a, eq.b, eq.c, eq.d, eq.S)
@inline equationDivChild1(y, eq::NodeEquation) = equationDivChild1(y, eq.a, eq.b, eq.c, eq.d, eq.S)
@inline equationDivChild2(x, eq::NodeEquation) = equationDivChild2(x, eq.a, eq.b, eq.c, eq.d, eq.S)


function semanticNotInDomain(semantic::Semantic, S::SemanticSet)
    isempty(S) && return false;
    isNaN(S) && return true;
    for s in S
        # anynan(s) && return false;
        anyequal(semantic,s) && return true;
    end;
    return false;
end;
semanticNotInDomain(semantic::Semantic, eq::NodeEquation) = semanticNotInDomain(semantic, eq.S)


function showEquation(eq::NodeEquation, semantic=nothing)
    println("a= ", eltype(eq.a), ".(", eq.a, ");");
    println("b= ", eltype(eq.b), ".(", eq.b, ");");
    println("c= ", eltype(eq.c), ".(", eq.c, ");");
    println("d= ", eltype(eq.d), ".(", eq.d, ");");
    # println("S= ", isempty(eq.S) ? "∅" : eq.S , ";");
    if isempty(eq.S)
        println("S= ", isempty(eq.S) ? "∅" : eq.S , ";");
    elseif isa(eq.S, Real) && isnan(eq.S)
        println("S= NaN;");
    else
        text = string(eq.S[1]);
        for i = 2:length(eq.S)
            text = string(text, ", ", string(eq.S[i]));
        end;

        println("S= {", text , "}");
    end;
    (semantic!=nothing) && println("semantic= ", semantic, ";");
end;



function checkEquation(eq::NodeEquation)
    @assert(isa(eq.a,Vector) || isa(eq.b,Vector));
    @assert(eltype(eq.a)==eltype(eq.b)==eltype(eq.c)==eltype(eq.d));
    # checkVector(v) = (isa(v,Number) || !allequal(v)) && !anyNaN(v) && !anyInf(v);
    checkVector(v) = !anyNaN(v) && !anyInf(v);
    @assert(checkVector(eq.a));
    @assert(checkVector(eq.b));
    @assert(checkVector(eq.c));
    @assert(checkVector(eq.d));
    if (isempty(eq.S))
        @assert(all0(eq.c))
    else
        if (isa(eq.S,Number))
            @assert(isnan(eq.S));
        else
            @assert(isa(eq.S,Array));
            @assert(sum(eq.S .== 0)<=1);
            for s in eq.S
                @assert(!anyNaN(s));
            end;
            poles = eq.S[end];
            if (isa(poles,Number) && isinf(poles))
                @assert(all0(eq.c));
            elseif (isa(poles,Number) && (poles==0))
                @assert(all0(eq.d));
            else
                if (anyInf(poles))
                    poles = poles[.!isinf.(poles)]
                end;
                if (!isempty(poles))
                    poles = cleanVector(poles);
                end;
                @assert(!anyNaN(poles));

                @assert(!any0(eq.c,eq.d));
                poles2 = eq.d./eq.c;

                if (isa(poles2,Vector))
                    if (anyInf(poles2))
                        poles2 = poles2[.!isinf.(poles2)]
                    end;
                    if (!isempty(poles2))
                        poles2 = cleanVector(poles2);
                    end;
                end;
                @assert(!anyNaN(poles2));
                @assert(equal(poles,poles2));
            end;
        end;
    end;
end;



function horizontalAsymptote(a,b,c,d,S; checkForErrors=false)
    if (isa(c,Number))
        checkForErrors && @assert(c!=0)
        return mean(a.^2)/(c^2);
    elseif (isa(a,Number))
        if (a==0)
            asymptote = b./d;
            if (isa(asymptote,Number))
                return (asymptote^2)*mean(c.==0)
            else
                asymptote[c.!=0] .= 0;
                return mean(asymptote.^2);
            end;
        else
            any0(c) && return NaN;
            return mean((a./c).^2);
        end;
    else
        # a and c are vectors
        asymptote = a./c;
        for i in 1:length(c)
            if (c[i]==0)
                (a[i]!=0) && return NaN
                asymptote[i] = (isa(b,Vector) ? b[i] : b) / (isa(d,Vector) ? d[i] : d);
            else
                if (a[i]==0)
                    asymptote[i] = 0;
                end;
            end;
        end;
    end;
    if (checkForErrors)
        @assert(!anyNaN(asymptote))
        @assert(!anyInf(asymptote))
    end;
    return mean(asymptote.^2);
end;
horizontalAsymptote(eq::NodeEquation; checkForErrors=false) = horizontalAsymptote(eq.a,eq.b,eq.c,eq.d,eq.S; checkForErrors=checkForErrors);


#############################################################################################################
#
# Functions for calculating the MSE from the equation values and the semantics of this node
#

calculateMSEFromEquation(semantic::Semantic, eq::NodeEquation; checkForErrors=false) = mean(( ((semantic.*eq.a).-eq.b)./((semantic.*eq.c).-eq.d) ).^2);

function calculateMSEFromEquationEfficient(semantic::Semantic, eq::NodeEquation; checkForErrors=false)
    a=eq.a; b=eq.b; c=eq.c; d=eq.d; S=eq.S;
    # Simplified common form of this equation:
    if all0(c)
        numerator = (semantic.*a).-b;
        if any0(d)
            any0(numerator,d) && return NaN;
            return Inf;
        end;
        if (isa(d,Number))
             # d is a constant
            if (checkForErrors) @assert(d!=0); end;
            result = mean(numerator.^2)/(d^2);
            # result = mean((mulFunction(semantic,a) .- b).^2) ./ (d^2)
        else
            # d is a vector
            result = mean((numerator./d).^2);
            # result = mean(((mulFunction(semantic,a).-b)./d).^2);
        end;
    elseif all1(c) && all0(d)
        if any0(semantic)
            any0(semantic,b) && return NaN;
            return Inf;
        end;
        # result = mean(((a.*semantic .- b)./semantic).^2);
        result = mean((a .- (b./semantic)).^2);
    elseif all0(d)
        if any0(semantic)
            any0(semantic,b) && return NaN;
            return Inf;
        end;
        any0(c) && return Inf;
        result = mean(((a.*semantic .- b)./(c.*semantic)).^2);
    else
        numerator = (semantic.*a).-b;
        denominator = (semantic.*c).-d;
        any0(numerator,denominator) && return NaN;
        any0(denominator) && return Inf;
        result = mean( (numerator./denominator).^2 );
    end;
    return result;
end;


function calculateMSEFromEquationMoreEfficient(semantic::Semantic, a,b,c,d,S; checkForErrors=false)
    # This function is similar to the one above, except for the fact that NaNs and Infs are returned as NaNs
    # Also, vectorized operations have been replaced by loops, so no memory allocations are needed
    #  Therefore, since less checks and memory allocations are done, it is slightly faster
    if all0(c)
        (any0(d)) && return NaN;
        if allequal(d)
            # d is a constant
            d = d[1];
            # return mean(((semantic.*a).-b).^2)/(d^2);

            n = max(length(a),length(b)); result = eltype(a)(0.);
            @inbounds for i in 1:n
                result += (valueSemantic(semantic,i)*valueSemantic(a,i) - valueSemantic(b,i))^2;
            end;
            result /= n*(d^2);

            if (checkForErrors) @assert(equal(result, mean(((semantic.*a).-b).^2)/(d^2))); end;
            return result;

            # result = mean((mulFunction(semantic,a) .- b).^2) ./ (d^2)
        else
            # d is a vector
            # return mean((((semantic.*a).-b)./d).^2);

            n = length(b); result = eltype(a)(0.);
            @inbounds for i in 1:n
                result += ((valueSemantic(semantic,i)*valueSemantic(a,i) - valueSemantic(b,i)) / valueSemantic(d,i)) ^2;
            end;
            result /= n;

            if (checkForErrors) @assert(equal(result, mean((((semantic.*a).-b)./d).^2))); end;
            return result;

            # result = mean(((mulFunction(semantic,a).-b)./d).^2);
        end;
    elseif all1(c) && all0(d)
        # any0(semantic) && return NaN;
        # result = mean(((a.*semantic .- b)./semantic).^2);
        # return mean((a .- (b./semantic)).^2);

        n = max(length(a),length(b)); result = eltype(a)(0.);
        @inbounds for i in 1:n
            semantic_i = valueSemantic(semantic,i);
            (semantic_i==0) && return NaN;
            result += (valueSemantic(a,i) - valueSemantic(b,i)/semantic_i) ^2;
        end;
        result /= n;

        if (checkForErrors) @assert(equal(result, mean((a .- (b./semantic)).^2))); end;
        return result;

    elseif all0(d)
        # any0(semantic) && return NaN;
        # any0(c) && return NaN;
        # return mean((((a.*semantic) .- b)./(c.*semantic)).^2);

        n = max(length(a),length(b)); result = eltype(a)(0.);
        @inbounds for i in 1:n
            semantic_i = valueSemantic(semantic,i);
            (semantic_i==0) && return NaN;
            c_i = valueSemantic(c,i);
            (c_i==0) && return NaN;
            result += ((valueSemantic(a,i)*semantic_i - valueSemantic(b,i))/(c_i*semantic_i) ) ^2;
        end;
        result /= n;

        if (checkForErrors) @assert(equal(result, mean((((a.*semantic) .- b)./(c.*semantic)).^2))); end;
        return result;
    else
        # denominator = (semantic.*c).-d;
        # any0(denominator) && return NaN;
        # return mean( (((semantic.*a).-b)./denominator).^2 );

        n = max(length(a),length(b),length(c),length(d)); result = eltype(a)(0.);
        @inbounds for i in 1:n
            denominator = (valueSemantic(c,i)*valueSemantic(semantic,i) - valueSemantic(d,i));
            (denominator==0) && return NaN;
            result += ((valueSemantic(a,i)*valueSemantic(semantic,i) - valueSemantic(b,i))/denominator)^2;
        end;
        result /= n;

        denominator = (semantic.*c).-d;
        any0(denominator) && return NaN;
        if (checkForErrors)
            denominator = (semantic.*c).-d;
            @assert(!any0(denominator))
            @assert(equal(result, mean( (((semantic.*a).-b)./denominator).^2 )));
        end;
        return result;
    end;
    return nothing;
end;
calculateMSEFromEquationMoreEfficient(semantic::Semantic, eq::NodeEquation; checkForErrors=false) = calculateMSEFromEquationMoreEfficient(semantic, eq.a,eq.b,eq.c,eq.d,eq.S; checkForErrors=checkForErrors);



#############################################################################################################
#
# Function for calculating the derivative of the equation on a value x
#

function derivativeEquation(eq::NodeEquation, x::Real)
    # return 2 .* mean( (b.*c .- a.*d) .* (a.*x .- b ) ./ ((c.*x .- d).^3) );
    a=eq.a; b=eq.b; c=eq.c; d=eq.d; S=eq.S;

    if anyDenominatorHas0Coefficients(c,d)
        return Inf;
    end;
    return 2 .* mean( (b.*c .- a.*d) .* (a.*x .- b ) ./ ((c.*x .- d).^3) );

    (a,b,c,d) = equation;
    bc_ad = (b.*c .- a.*d);
    result = bc_ad .* (a.*x .- b ) ./ ((c.*x .- d).^3);
    result[bc_ad.==0] .= 0.;
    return 2*mean(result);
end;

#############################################################################################################
#
# Function for calculating the reduction in MSE
#

function calculateMSEReduction(semantic::Semantic, equation::NodeEquation, mse::Real; checkForErrors=false)
    if (checkForErrors)
        @assert(!all0(equation.c) || !all0(equation.d));
        @assert(!anyNaN(semantic));
        @assert(!anyInf(semantic));
        @assert(!isnan(mse) && !isinf(mse));
    end;
    newMSE = calculateMSEFromEquationMoreEfficient(semantic, equation; checkForErrors=checkForErrors);
    if (checkForErrors)

        @assert(!isinf(newMSE));
        newMSE2 = calculateMSEFromEquationEfficient(semantic, equation; checkForErrors=checkForErrors);
        @assert(isnan(newMSE) == (isnan(newMSE2) || isinf(newMSE2)))
        if (!isnan(newMSE))
            @assert(equal(newMSE,newMSE2));
        end;
        newMSE3 = calculateMSEFromEquation(semantic, equation; checkForErrors=checkForErrors);
        @assert(isnan(newMSE2) == isnan(newMSE3));
        @assert(isinf(newMSE2) == isinf(newMSE3));
        if isinf(newMSE2)
            @assert(sign(newMSE2)==sign(newMSE3))
        elseif (!isnan(newMSE))
            @assert(equal(newMSE,newMSE2));
        end;
    end;
    return mse - newMSE;
end;


#############################################################################################################
#
# Functions for calculating the best constant for a node
#  Return the constant and the reduction in MSE
#

function calculateRoots(c,d; checkForErrors=false)
    (c==0) && return Array{eltype(c),1}();
    roots = d./c;
    if any0(c)
        # c must be a vector
        (checkForErrors) && @assert(isa(c,Vector));
        roots=roots[c.!=0];
    end;
    roots = unique(roots);
    if (checkForErrors)
        @assert(!any(isnan.(roots)));
        @assert(!any(isinf.(roots)));
        @assert(isa(roots,Vector));
    end;
    return roots;
end;



function calculateConstantGeneralCase(a,b,c,d,S; checkForErrors=false)

    if (checkForErrors)
        @assert(!isempty(S))
        @assert(!isNaN(S))
        @assert(!anyNaN(S[end]))
    end;

    zeros = calculateRoots(a,b; checkForErrors=checkForErrors);
    # These are not zeros of the whole equation, but zeros in each term

    isempty(zeros) && return NaN;
    zeros = zeros[.!semanticNotInDomain.(zeros, [S])];

    if (checkForErrors)
        @assert(!any(isnan.(zeros)));
    end;

    isempty(zeros) && return NaN;


    # Avoid those zeros which are too close to a root of the denominator
    # The roots of the denominator should already be calculated in the last term of S
    poles = S[end];
    if (isa(poles,Number) && isinf(poles))
        # In this case, any value is valid (there are no roots in the denominator)
        if (checkForErrors)
            @assert(all0(c));
            # @assert(isempty(cleanEquation(calculateRoots(c,d; checkForErrors=checkForErrors))));
            @assert(isempty(calculateRoots(c,d; checkForErrors=checkForErrors)));
        end;
    else

        # However, in order to avoid possible rounding errors, calculate the roots of the denominator instead of S[end]
        poles = calculateRoots(c,d; checkForErrors=checkForErrors);

        if (checkForErrors);
            @assert(!isempty(poles));
            @assert(!anyInf(poles));
            @assert(!anyNaN(poles));
        end;

        # poles = cleanVector(poles);

        function isNotClose(zero, poles, tolerance)
            for pole in poles
                limit = abs(pole*tolerance);
                (zero>=pole-limit) && (zero<=pole+limit) && return false;
            end;
            return true;
        end;
        zeros = zeros[isNotClose.(zeros,[poles],[eltype(a)(1e-2)])];

        isempty(zeros) && return NaN;

    end;

    yZeros = calculateMSEFromEquationMoreEfficient.(zeros, [a],[b],[c],[d],[S]; checkForErrors=checkForErrors);

    if (checkForErrors)
        @assert(!anyNaN(yZeros));
        @assert(!anyInf(yZeros));
    end;
    (_,index) = findmin(yZeros);
    zero = zeros[index];

    return zero;

end;
calculateConstantGeneralCase(eq::NodeEquation; checkForErrors=false) = ( a=eq.a; b=eq.b; c=eq.c; d=eq.d; S=eq.S; calculateConstantGeneralCase(a,b,c,d,S; checkForErrors=checkForErrors); )


function calculateConstant(a,b,c,d,S; checkForErrors=false)
    any0(c,d) && return NaN;
    lengthArray = maximum([length(a), length(b), length(c), length(d)]);
    # anyDenominatorHas0Coefficients(c,d) && return NaN;
    createVector(v::Vector) = v;
    createVector(x::Number) = ( v = Array{eltype(a),1}(undef,lengthArray); v.=x; return v; )
    a = createVector(a);
    b = createVector(b);
    c = createVector(c);
    d = createVector(d);

    if all0(c)
        any0(d) && return NaN;
        if (isa(a,Vector))
            index0 = (a.==0);
            num = (b.*a)./(d.^2);
            num[index0] .= eltype(a)(0.);
            den = ((a./d).^2);
            den[index0] .= eltype(a)(0.);
            constant = sum(num)./sum(den);
        else
            constant = sum((b.*a)./(d.^2))/sum((a./d).^2);
        end;
        semanticNotInDomain(constant, S) && return NaN;
    elseif all0(d)
        any0(c) && return NaN;
        if (isa(b,Vector))
            index0 = (b.==0);
            num = ((b./c).^2);
            num[index0] .= eltype(a)(0.);
            den = ((a.*b)./(c.^2));
            den[index0] .= eltype(a)(0.);
            constant = sum(num)/sum(den);
        else
            constant = sum((b./c).^2) ./ sum(a.*b./(c.^2));
        end;
        semanticNotInDomain(constant, S) && return NaN;
    elseif allequal(c) && allequal(d)
        mult = (b.*c .- a.*d);
        constant = sum(mult.*b)/sum(mult.*a);
        semanticNotInDomain(constant, S) && return NaN;
    else
        constant = calculateConstantGeneralCase(a,b,c,d,S; checkForErrors=checkForErrors);
    end;
    return constant;
end;
calculateConstant(eq::NodeEquation; checkForErrors=false) = ( a=eq.a; b=eq.b; c=eq.c; d=eq.d; S=eq.S; calculateConstant(a,b,c,d,S; checkForErrors=checkForErrors); )


function calculateConstantEfficient(a,b,c,d,S; checkForErrors=false)
    # Vectorial operations have been replaced by loops, so no memory allocation is needed and all the operations are much faster
    if all0(c)
        all0(a) && return NaN;
        any0(d) && return NaN;
        if allequal(a)
            # a is constant
            a = a[1];
            if (checkForErrors) @assert(isa(b,Vector)); end;
            # if isa(d,Number)
            if allequal(d)
                # d is constant
                d = d[1];
                # constant = mean(b)/a;
                constant = eltype(a)(0.);
                @inbounds for i in 1:length(b)
                    constant += b[i];
                end;
                constant /= (length(b)*a);

                if (checkForErrors) @assert(equal(constant,sum((b.*a)./(d.^2))/(length(b)*((a./d).^2)))); end;
            else
                # d is vector
                # inv_d = (1. ./ d).^2; constant = sum(b.*inv_d)/ (a*sum(inv_d));

                numerator = eltype(a)(0.); denominator = eltype(a)(0.);
                @inbounds for i in 1:length(d)
                    v = (1/d[i])^2;
                    numerator += valueSemantic(b,i)*v;
                    denominator += v;
                end;
                constant = numerator / (a*denominator);

                if (checkForErrors) @assert(equal(constant,sum((b.*a)./(d.^2))/sum((a./d).^2))); end;
            end;
        else
            # a is vector
            if allequal(d)
                # d is constant
                d = d[1];
                # constant = sum(b.*a)/sum(a.^2);

                numerator = eltype(a)(0.); denominator = eltype(a)(0.);
                @inbounds for i in 1:length(a)
                    a_i = valueSemantic(a,i);
                    numerator += a_i*valueSemantic(b,i);
                    denominator += (a_i^2);
                end;
                constant = numerator / denominator;

                if (checkForErrors) @assert(equal(constant,sum((b.*a)./(d.^2))/sum((a./d).^2))); end;
            else
                # d is vector
                # constant = sum((b.*a)./(d.^2))/sum((a./d).^2);

                numerator = eltype(a)(0.); denominator = eltype(a)(0.);
                @inbounds for i in 1:length(a)
                    a_i = valueSemantic(a,i);
                    d_i = valueSemantic(d,i);
                    numerator += a_i*valueSemantic(b,i)/(d_i^2);
                    denominator += (a_i/d_i)^2;
                end;
                constant = numerator / denominator;

                if (checkForErrors) @assert(equal(constant,sum((b.*a)./(d.^2))/sum((a./d).^2))); end;

            end;
        end;
        semanticNotInDomain(constant, S) && return NaN;

    elseif all0(d)
        any0(c) && return NaN;
        if allequal(b)
            # b is constant
            b = b[1];
            if (checkForErrors) @assert(isa(a,Vector)); end;
            # if isa(c,Number)
            if allequal(c)
                # c is constant
                c = c[1];
                # constant = b./mean(a);

                constant = eltype(a)(0.);
                @inbounds for i in 1:length(a)
                    constant += a[i];
                end;
                constant = length(a)*b/constant;

                if (checkForErrors) @assert(equal(constant,length(a)*((b./c).^2) ./ sum(a.*b./(c.^2)))); end;
            else
                # c is vector
                # inv_c = (1. ./ c).^2;
                # constant = b.*sum(inv_c)/sum(a.*inv_c);

                numerator = eltype(a)(0.); denominator = eltype(a)(0.);
                @inbounds for i in 1:length(c)
                    inv_c_i = (1/c[i])^2;
                    numerator += inv_c_i;
                    denominator += (valueSemantic(a,i)*inv_c_i);
                end;
                constant = b * numerator / denominator;

                if (checkForErrors) @assert(equal(constant,sum((b./c).^2) ./ sum(a.*b./(c.^2)))); end;
            end;
            semanticNotInDomain(constant, S) && return NaN;
        else
            # b is vector
            if allequal(c)
                # c is constant
                c = c[1];
                # constant = sum(b.^2)/sum(a.*b);

                numerator = eltype(a)(0.); denominator = eltype(a)(0.);
                @inbounds for i in 1:length(b)
                    b_i = b[i];
                    numerator += b_i^2;
                    denominator += (valueSemantic(a,i)*b_i);
                end;
                constant = numerator / denominator;

                if (checkForErrors) @assert(equal(constant,sum((b./c).^2) ./ sum(a.*b./(c.^2)))); end;
            else
                # c is vector
                # constant = sum((b./c).^2) ./ sum(a.*b./(c.^2));

                numerator = eltype(a)(0.); denominator = eltype(a)(0.);
                @inbounds for i in 1:length(c)
                    b_i = valueSemantic(b,i);
                    numerator += (b_i/c[i])^2;
                    denominator += (valueSemantic(a,i)*b_i/(c[i]^2));
                end;
                constant = numerator / denominator;

                if (checkForErrors) @assert(equal(constant,sum((b./c).^2) ./ sum(a.*b./(c.^2)))); end;
            end;
            semanticNotInDomain(constant, S) && return NaN;
        end;

    elseif allequal(c) && allequal(d)
        # constant = sum(b.*b.*c .- a.*b.*d)/sum(a.*b.*c .- a.*a.*d);

        # mult = (b.*c .- a.*d);
        # constant = sum(mult.*b)/sum(mult.*a);

        c = c[1]; d = d[1];
        numerator = eltype(a)(0.); denominator = eltype(a)(0.);
        @inbounds for i in 1:max(length(a),length(b))
            a_i = valueSemantic(a,i);
            b_i = valueSemantic(b,i);
            mult = b_i*c - a_i*d;
            numerator += b_i*mult;
            denominator += a_i*mult;
        end;
        constant = numerator / denominator;

        if (checkForErrors)
            mult = (b.*c .- a.*d);
            @assert(!all0(a));
            @assert(!all0(mult));
            @assert(equal(constant,sum(mult.*b)/sum(mult.*a)));
        end;

        semanticNotInDomain(constant, S) && return NaN;

    else
        any0(c,d) && return NaN;
        return calculateConstantGeneralCase(a,b,c,d,S; checkForErrors=checkForErrors);
    end;

    if (checkForErrors)
        @assert(!isnan(constant));
        @assert(!isinf(constant));
    end;
    return constant;
end;
calculateConstantEfficient(eq::NodeEquation; checkForErrors=false) = calculateConstantEfficient(eq.a,eq.b,eq.c,eq.d,eq.S; checkForErrors=checkForErrors);



function calculateConstantMinimizeEquation(equation::NodeEquation, mse::Real; checkForErrors=false)
    constant = calculateConstantEfficient(equation; checkForErrors=checkForErrors);
    if (checkForErrors)
        @assert(!isnan(mse) && !isinf(mse));
        @assert(!isinf(constant));
        newConstant = calculateConstant(equation; checkForErrors=checkForErrors);
        @assert(isnan(constant)==isnan(newConstant));
        if (!isnan(constant))
            @assert(equal(constant,newConstant));
        end;
        @assert(!isinf(constant))
    end;
    if (isnan(constant))
        reduction = -Inf;
    else
        reduction = calculateMSEReduction(constant, equation, mse; checkForErrors=checkForErrors);
    end;
    reduction = isnan(reduction) ? -Inf : reduction;
    return (constant, reduction);
end;

mutable struct foo
    a::Float64
    b::Float64
    c::Float64 #a+b 

    foo(a,b) = new(a,b,a+b^2)

    function foo(a::Float64, b::Float64)
        if a> b 
            k = a^10

        else
            k = b^10

        end
        return new(a,b,k)
    end
end

 a = foo(1,2)
 a = foo(1.0,2.0)
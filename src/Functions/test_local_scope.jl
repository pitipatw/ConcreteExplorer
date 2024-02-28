a = 1 
for i in 1:4
    a = 2
end

println(a)


function foo()
    a = 1 
    for i in 1:4
        local a = 2
    end
    println(a)
return nothing
end

foo()
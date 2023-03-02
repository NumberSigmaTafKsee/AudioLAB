function objective_function(v)
    local sum = 0
    for i=1,#v do
        sum = sum + v[i]^2
    end
    return sum
end

function rand_in_bounds(min,max)
    return min + (max-min)*math.random()
end

function random_vector(size,min,max)
    local v = {}
    for i=1,size do
        v[i] = min + (max-min)*math.random()
    end
    return v
end

function Max(a,b)
    if( a > b ) then return a end
    return b
end

function Min(a,b)
    if( a < b ) then return a end
    return b
end

function take_step(min,max,current,step_size)
    local position = {}
    for i=1,#current do
        local _min = Max(min,current[i] - step_size)
        local _max = Min(max,current[i] + step_size)
        position[i] = rand_in_bounds(_min,_max)
    end
    return position
end

function large_step_size(iter,step_size,s_factor,l_factor,iter_mul)
    if(iter > 0 and iter % iter_mul == 0) then
        return step_size * l_factor
    else
        return step_size * s_factor
    end
end

function take_steps(min,max,current,step_size,big_stepsize)
    local step = {}
    local big_step = {}
    step.vector = take_step(min,max,current.vector,step_size)
    step.cost   = objective_function(step.vector)
    big_step.vector = take_step(min,max,current.vector,big_stepsize)
    big_step.cost   = objective_function(big_step.vector)
    return step,big_step
end

function search(max_iter,bounds,init_factor,s_factor,l_factor,iter_mult,max_no_impr)
    local step_size = (bounds[2] - bounds[1]) * init_factor
    local current   = {}
    local count     = 0
    current.vector = random_vector(2,bounds[1],bounds[2])
    current.cost   = objective_function(current.vector)
    for iter=1,max_iter do
        local big_stepsize = large_step_size(iter,step_size,s_factor,l_factor,iter_mult)    
        local step,big_step = take_steps(bounds[1],bounds[2],current,step_size,big_stepsize)
        if(step.cost <= current.cost or big_step.cost <= current.cost) then
            if(big_step.cost <= step.cost) then
                step_size,current = big_stepsize, big_step
            else
                current = step
            end
            count = 0
        else
            count = count + 1
            if(count >= max_no_impr) then                
                count,step_size = 0,step_size/s_factor
            end            
        end
        print("iteration #",iter-1," best=", current.cost)
    end
    return current
end

function print_vector(v)
    for i=1,#v do
        io.write(v[i],",")    
    end
    print()
end

local problem_size=2
local bounds = {-5,5}
local max_iter=1000
local init_factor = 0.05
local s_factor = 1.3
local l_factor = 3.0
local iter_mult=10
local max_no_impr=30
best = search(max_iter,bounds,init_factor,s_factor,l_factor,iter_mult,max_no_impr)
print("Best solution c=",best.cost)
print_vector(best.vector)
    



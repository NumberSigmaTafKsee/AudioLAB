
function objective_function(v)
    local sum = 0
    for i=1,#v do
        sum = sum + v[i]^2
    end
    return sum
end

function random_vector(size,min,max)
    local v = {}
    for i=1,size do
        v[i] = min + (max-min)*math.random()
    end
    return v
end

function search(search_space,max_iter,show_log)
    local best = nil
    local old  = nil
    for iter=1,max_iter do
        local candidate = {}        
        candidate.vector = random_vector(search_space.size,search_space.min,search_space.max)
        candidate.cost  = objective_function(candidate.vector)
        if(best == nil or candidate.cost < best.cost) then 
            best = candidate
        end
        if(show_log == true) then 
            print("iteration:",iter,"best=",best.cost)
        end
    end
    return best
end

function print_vector(v)
    for i=1,#v do
        io.write(v[i],",")    
    end
    print()
end

problem_size = 2
search_space = {}
search_space.size = 2
search_space.min  = -5
search_space.max  = 5
max_iter = 100
show_log = true
best = search(search_space,max_iter,show_log)
if(show_log) then
    print("Best solution cost=",best.cost," value=")
    print_vector(best.vector)
end
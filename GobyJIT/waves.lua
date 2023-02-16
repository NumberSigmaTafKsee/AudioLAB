function sine_wave(v)
    -- sine wave
    for i=1,4096 do 
        v[i-1] = math.sin(2*math.pi*phase/4096.0)
        phase = phase + 1
    end 
end 

function triangle_wave(f,v)
    -- triangle wave
    pi8 = 8.0 / (math.pi*math.pi)
    phase = 0.0    
    for i = 1,4096 do 
        x = 0.0
        -- 0 to max harmonic must be calculated from frequency
        for n = 0,(44100/f) do             
            q = math.pow(-1.0,n)            
            h = 2*n+1
            r = math.pow(h,-2.0)
            x = x + q*r*math.sin(2*math.pi*h*phase/4096.0)
        end 
        v[i-1] = pi8*x 
        phase = phase + 1 
    end
end

function saw_wave(f,v)    
    phase = 0.0    
    for i = 1,4096 do 
        x = 0.0
        -- 0 to max harmonic must be calculated from frequency
        for n = 1,(44100/f) do             
            q = math.pow(-1.0,n)                        
            x = x + q*(math.sin(2*math.pi*n*phase/4096.0)/n)
        end 
        v[i-1] = 0.5-(1/math.pi)*x
        phase = phase + 1 
    end
end

function reverse_saw_wave(f,v)    
    phase = 0.0    
    for i = 1,4096 do 
        x = 0.0
        -- 0 to max harmonic must be calculated from frequency
        for n = 1,(44100/f) do             
            q = math.pow(-1.0,n)                        
            x = x + q*(math.sin(2*math.pi*n*phase/4096.0)/n)
        end 
        v[i-1] = (1/math.pi)*x 
        phase = phase + 1 
    end
end

function square_wave(f,v)
    phase = 0.0    
    for i = 1,4096 do 
        x = 0.0
        -- 0 to max harmonic must be calculated from frequency
        for n = 1,(44100/f) do          
            h = 2*n-1               
            x = x + math.sin(2*math.pi*h*phase/4096.0)/h
        end 
        v[i-1] = (4/math.pi)*x 
        phase = phase + 1 
    end
end
    
function pulse_wave(f,v,t)
    phase = 0.0    
    for i = 1,4096 do 
        x = 0.0
        -- 0 to max harmonic must be calculated from frequency
        for n = 1,(44100/f) do          
            q = n * t
            h = math.sin(math.pi*q)/(math.pi*q)
            x = x + math.sin(h)*math.cos(2*math.pi*n*phase/4096.0)
        end 
        v[i-1] = 1+2*x 
        phase = phase + 1 
    end
end
    
function gauss_white_noise()
    local R1 = math.random()
    local R2 = math.random()
    return math.sqrt( -2.0 * math.log( R1 )) * math.cos( 2.0 * math.pi * R2 );
end 


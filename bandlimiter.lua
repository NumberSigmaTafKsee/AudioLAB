-- neural bandlimiter 

require('luapa')
require('Cranium')
require('LuaPolyBLEP')

NUM=8
sin = LuaPolyBLEP.PolyBLEP(44100,LuaPolyBLEP.PolyBLEP.SAWTOOTH)
network = Network_New(NUM, 1, {{size=2*NUM, activation=cranium.active_tanh }}, 1, cranium.active_tanh)
data = {} 
for i = 1,NUM do data[i] = 0 end 
example = DataSet_New(1, NUM, {data})
classes = DataSet_New(1, 1,{{0}})
 
local params = cranium.ParameterSet()    
params.network = network.network
params.data     = example.ds
params.classes  = classes.ds
params.lossFunction = 1
params.batchSize    = 1
params.learningRate = .1
params.searchTime   = 1000
params.regularizationStrength = 1.0e-9
params.momentumFactor = 0.9
params.maxIters = 10
params.shuffle  = 0
params.verbose  = 0

f = 220
sin:setFrequency(f)

function train()
    local cur1=0,cur2
        
    for j=1,NUM do                             
        cur1 = sin:getAndInc()
        data[j] = cur1            
    end        
        
    for i=-10,10 do           
        for j=1,NUM do 
            example:set(0,j-1,data[j])
        end            
        cur2 = sin:getAndInc()
        classes:set(0,0,cur2)
        cranium.optimize(params)                                                                                              
        table.remove(data)        
        table.insert(data,1,cur2)
        --table.insert(data,1,i/100)        
    end 
end 

train() 
f =110
sin:setFrequency(f)

function noise(input,output,frames)                
    for i = 0,frames-1 do                          
        for j=1,NUM do 
            example:set(0,j-1,data[j])
        end                        
        network:forwardPassDataSet(example)        
        local o = network:getOutput()        
        local x = o:get_index(0)/4        
        table.remove(data)                
        table.insert(data,1,2.0*(math.random()-0.5))        
        luapa.float_set(output,2*i,x)
        luapa.float_set(output,2*i+1,x)        
    end     
end 

function sys()
    print(x)
end 
luapa.set_audio_func(noise)    

luapa.InitAudioDevice(10,-1,2,44100,128)
luapa.RunAudio()
luapa.StopAudio()




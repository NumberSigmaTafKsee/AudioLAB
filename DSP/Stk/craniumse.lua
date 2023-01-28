#!~/bin/luajit
-- requires simpleeigen(se)
rng = require('random')
JSON= require('JSON')
csv = require('csv')
profiler=require('profiler')
require('se')

FLT_MIN=1.17549435e-38

use_profiler = false

INPUT=1
OUTPUT=2
HIDDEN=3

CROSS_ENTROPY_LOSS=1
MEAN_SQUARED_ERROR=2

math.randomseed(os.time(nil))

function matrix_new(rows,cols,data,foo)    
    local m = se.FloatMatrix(rows,cols)
    foo = foo or false
    for i=0,rows-1 do         
        for j=0,cols-1 do 
            local x
            if(foo == true) then x = data[i+1][j+1] 
            else  x = data[i][j] end             
            m:set(i,j,x)
        end 
    end     
    return m 
end 
function matrix_create(rows,cols)    
    local m =  se.FloatMatrix(rows,cols)
    m:zero() 
    return m
end

function createMatrixZeros(rows,cols)
    return matrix_create(rows,cols)    
end    

function matrix_getMatrix(m,row,col)    
    return m:get(row,col)
end 

function matrix_setMatrix(m,row,col,val)    
    m:set(row,col,val)
end 

function matrix_copyValuesInto(from,to)    
    to = from:eval()
    return to 
end 

function matrix_printMatrix(input)
    for i=0,input:rows()-1 do
        print() 
        for j=0,input:cols()-1 do 
            io.write( tostring(input:get(i,j)) )
            io.write(',')
        end 
        print() 
    end 
end 

function matrix_zero(m)
    m:zero()
    return m
end 

function zeros(rows,cols)
    local v = {} 
    for i = 0,rows*cols-1 do 
        v[i] = 0.0 
    end 
    return v 
end 

function matrix_transpose(m)
    return m:transpose()
end 
 
function matrix_transposeInto(m,origT)    
    m:transposeInto(origT)
    return origT
end 


 
function MAX(a,b)
    if( a > b) then return a 
    else return b 
    end
end 

function sigmoidFunc(input)
    return 1.0 / (1.0 + math.exp(-1 * input))
end 

function reluFunc(input)
    return MAX(0, input);
end 

function tanHFunc(input)    
    return math.tanh(input)
end 

function relu(input)    
    se.relu_float(input)
end 

function reluDeriv(reluInput)
    se.relu_deriv_float(reluInput)
end 

function sigmoid(input)
    se.sigmoid_float(input)
end 

function sigmoidDeriv(sigmoidInput)    
    se.sigmoid_deriv_float(sigmoidInput)
end 

function tanH(input)    
    se.tanh_float(input)
end 

function tanHDeriv(input)
    se.tanh_deriv_float(input)    
end 

function softmax(input)    
    se.softmax_float(input)    
end 

function linear(input)
    
end 

function linearDeriv(linearInput)
    
end 

function getFunctionName(f)
    if(f == sigmoid) then return "sigmoid"
    elseif(f == relu) then return "relu"
    elseif(f == tanH) then return "tanH"
    elseif(f == softmax) then return "softmax"
    end 
    return "linear"
end

function getFunctionByName(f)    
    if(f == "sigmoid") then return sigmoid
    elseif(f == "relu") then return relu
    elseif(f == "tanH") then return tanH
    elseif(f == "softmax") then return softmax
    end 
    return linear
end

function activationDerivative(func)
    if(func == sigmoid) then return sigmoidDeriv
    elseif(func == relu) then return reluDeriv 
    elseif(func == tanH) then return tanHDeriv 
    end 
    return linearDeriv
end 

function box_muller_new()
    local bm = {} 
    bm.z0 = 0 
    bm.z1 = 0
    bm.generate = 0 
    bm.box_muller = box_muller
    return bm 
end 

function box_muller(bm)
    local epsilon = FLT_MIN 
    local two_pi  = 2 * math.pi 
    if(bm.generate == 1) then bm.generate = 0 
    else bm.generate = 1 
    end
    if(bm.generate == 0) then return bm.z1 end 
    local u1 = math.random()
    local u2 = math.random()
    while(u1 <= epsilon) do 
        u1 = math.random()        
    end   
    while(u2 <= epsilon) do 
        u2 = math.random()        
    end   
    bm.z0 = math.sqrt(-2.0 * math.log(u1)) * math.cos(two_pi * u2);
    bm.z1 = math.sqrt(-2.0 * math.log(u1)) * math.sin(two_pi * u2);
    return bm.z0;
end

boxmuller = box_muller_new()

function layer_new(type, size, activation)
    local layer = {} 
    layer.type = type 
    layer.size = size 
    layer.activation = activation          
    layer.input = matrix_create(1,size)
    layer.activate = layerivateLayer
    return layer 
end 

function connection_InitializeConnection(connection)
    connection.bias:fill(1)    
    --connection.weights:random()
    for i=0,connection.weights:rows()-1 do 
        for j=0,connection.weights:cols()-1 do 
            local neuronsIn = connection.weights:rows()            
            connection.weights:set(i,j,boxmuller:box_muller() / math.sqrt(neuronsIn))
        end        
    end    
end 

function connection_new(from,to)
    local connection = {} 
    connection.from = from 
    connection.to = to     
    connection.weights = matrix_create(from.size,to.size)
    connection.bias = matrix_create(1, to.size)
    connection.initializeConnection = connection_InitializeConnection    
    return connection 
end 

function layerivateLayer(layer)
    if(layer.activation ~= nil) then 
        layer.activation(layer.input)
    end 
end 

function parameterset_new(network,data,classes, lossFunction, batchSize, learningRate, searchTime, regularizationStrength, momentumFactor,maxIters,shuffle,verbose)
    local ps = {} 
    ps.network = network or nil 
    ps.data = data or nil 
    ps.classes = classes or nil
    ps.lossFunction = lossFunction or MEAN_SQUARED_ERROR
    ps.batchSize = batchSize or 1
    ps.learningRate = learningRate or 0.001
    ps.searchTime = searchTime or 1000
    ps.regularizationStrength = regularizationStrength or 1e-6
    ps.momentumFactor = momentumFactor or 0.9 
    ps.maxIters = maxIters or 1000
    ps.shuffle = shuffle or false 
    ps.verbose = verbose or true 
    return ps 
end 

function network_new(numFeatures, numHiddenLayers, hiddenSizes, hiddenActivations, numOutputs, outputActivation)
    local net = network_createNetwork(numFeatures, numHiddenLayers, hiddenSizes, hiddenActivations, numOutputs, outputActivation)
    net.saveNetwork = network_saveNetwork
    net.readNetwork = network_readNetwork    
    net.accuracy = network_accuracy
    net.predict = network_predict
    net.getOutput = network_getOutput    
    net.meanSquaredError = network_meanSquaredError
    net.crossEntropyLoss = network_crossEntropyLoss
    net.forwardPassDataSet = network_forwardPassDataSet
    net.forwardPass = network_forwardPass
    net.optimize = network_optimize
    return net
end 

function network_createNetwork(numFeatures, numHiddenLayers, hiddenSizes, hiddenActivations, numOutputs, outputActivation)
    assert(numFeatures > 0 and numHiddenLayers >= 0 and numOutputs > 0);
    local network = {} 
    network.numLayers = 2 + numHiddenLayers
    local layers = {}     
    for i=0,network.numLayers-1 do layers[i] = nil end 
    for i=0,network.numLayers-1 do 
        if(i == 0) then layers[i] = layer_new(INPUT,numFeatures,nil) 
        elseif(i == network.numLayers-1) then 
            layers[i] = layer_new(OUTPUT,numOutputs,outputActivation)
        else                         
            layers[i] = layer_new(HIDDEN, hiddenSizes[i], hiddenActivations[i])
        end 
    end 
    network.layers = layers 
    network.numConnections = network.numLayers-1
    local connections = {} 
    for i=0,network.numConnections do connections[i] = nil end 
    for i=0,network.numConnections-1 do 
        connections[i] = connection_new(network.layers[i],network.layers[i+1])        
        connections[i]:initializeConnection()
    end 
    network.connections = connections 
    return network 
end


function network_forwardPass(network,input)               
    assert(input:cols() == network.layers[0].input:cols())        
    network.layers[0].input = input:eval()        
    
    local tmp,tmp2     
    for i=0,network.numConnections-1 do                                          
        tmp  = network.layers[i].input * network.connections[i].weights                        
        tmp2 = tmp:addToEachRow(network.connections[i].bias)                       
        network.connections[i].to.input = tmp2:eval()
        network.connections[i].to:activate()            
    end 
end 

function network_forwardPassDataSet(network,input)
    --ocal dataMatrix = input:toMatrix()
    network:forwardPass(input)
end 

function network_crossEntropyLoss(network, prediction, actual, regularizationStrength)
    --assert(prediction.rows == actual.rows)
    --assert(prediction.cols == actual.cols)
    local total_err = 0;
    
    for i=0,prediction:rows()-1 do
        local cur_err = 0;
        for j=0,prediction:cols()-1 do        
            cur_err =  cur_err + actual:get(i,j) * math.log(MAX(FLT_MIN, prediction:get(i, j)))        
        end
        total_err = total_err + cur_err;
    end
    local reg_err = 0;
    for i=0,network.numConnections-1 do    
        local weights = network.connections[i].weights;
        for j=0,weights:rows()-1 do
            for k=0,weights:cols()-1 do            
                reg_err = reg_err + weights:get(j, k) * weights:get(j, k);
            end
        end
    end

    return ((-1.0 / actual:rows()) * total_err) + (regularizationStrength * .5 * reg_err);
end 

function network_meanSquaredError(network, prediction, actual, regularizationStrength)
    --assert(prediction:rows() == actual:rows())
    --assert(prediction:cols() == actual:cols())
    local total_err = 0
    local reg_err = 0 
    local tmp = actual - prediction    
    total_err = tmp:cwiseProduct(tmp:eval()):sum()
    for i=0,network.numConnections - 1 do 
        local w = network.connections[i].weights    
        reg_err = reg_err + w:cwiseProduct(w):sum()
    end 
    return ((0.5 / actual:rows()) * total_err) + (regularizationStrength*0.5*reg_err)
end 

function network_getOutput(network)
    return network.layers[network.numLayers-1].input:eval()
end 

function network_predict(network)
    local max = 0 
    local outputLayer = network.layers[network.numLayers-1]
    local prediction = {}
    for i=0,outputLayer.input:rows()-1 do 
        max = 0
        for j=0, outputLayer.input:cols()-1 do             
            if(outputLayer.input[i][j] > outputLayer.input[i][max]) then 
                max = j 
            end 
        end 
        prediction[i] = max 
    end 
    return prediction
end 

function network_accuracy(network, data, classes)
    assert(data:rows() == classes:rows())
    assert(classes:cols() == network.layers[network.numLayers-1].size)
    network:forwardPassDataSet(data)
    local p = network:predict() 
    local numCorrect = 0 
    for i=0, data:rows()-1 do 
        if(classes[i][p[i]] == 1) then 
            numCorrect = numCorrect + 1 
        end 
    end 
    return 100*numCorrect/(classes:rows())
end 



function network_save(network, path)
    local fp = io.open(path,"w")
    fp:write(network.numLayers,"\n")
    for i=0,network.numLayers-1 do 
        fp:write(network.layers[i].size,"\n")
    end 
    for i=0,network.numLayers-1 do         
        local x = getFunctionName(network.layers[i].activation)        
        fp:write(x,"\n")
    end 
    for k=0,network.numConnections-1 do 
        local con = network.connections[k]
        for i=0,con.weights:rows()-1 do 
            for j=0,con.weights:cols()-1 do                 
                    fp:write(con.weights:get(i,j),"\n")
            end 
        end 
    end 
    for k=0,network.numConnections-1 do 
        local con = network.connections[k]
        for i=0,con.bias:cols()-1 do             
            fp:write(con.bias:get(0,i),'\n')
        end 
    end 
    fp:close() 
end 

function network_load(path)
    local fp = io.open(path,"r")    
    local numLayers = tonumber(fp:read("*line"))
    
    local layerSizes={}
    for i=1,numLayers do 
        layerSizes[i] = tonumber(fp:read("*line"))
    end 
    local funcs = {} 
    for i=1,numLayers do 
        local x= fp:read('*line')        
        funcs[i] = getFunctionByName(x)          
    end 
    local network 
    local inputSize = layerSizes[1] 
    local outputSize = layerSizes[numLayers]
    local numHiddenLayers = numLayers - 2 
    local outputFunc = funcs[numLayers]
    local hiddenSizes = {}     
    local hiddenFuncs = {} 
    if(numHiddenLayers > 0) then        
        for i=1,numLayers-2 do             
            hiddenSizes[i] = layerSizes[i]                    
            hiddenFuncs[i] = funcs[i]            
        end         
    end 
    print(inputSize,numHiddenLayers,hiddenSizes[1],hiddenFuncs[1],outputSize,outputFunc)
    network = network_new(inputSize, numHiddenLayers, hiddenSizes, hiddenFuncs, outputSize, outputFunc)
    for k=0,network.numConnections-1 do 
        local con = network.connections[k]
        for i=0,con.weights:rows()-1 do 
            for j=0,con.weights:cols()-1 do 
                con.weights:set(i,j,tonumber(fp:read("*line")))
            end 
        end 
    end 

    for k=0,network.numConnections-1 do 
        local con = network.connections[k]        
        for i=0,con.bias:cols()-1 do 
            con.bias:set(0,i,tonumber(fp:read("*line")))
        end 
    end 
    fp:close()
    return network 
end

function network_optimize(net,params)
    batchGradientDescent(params.network, params.data, params.classes, params.lossFunction, params.batchSize, params.learningRate, params.searchTime, params.regularizationStrength, params.momentumFactor, params.maxIters, params.shuffle, params.verbose)
end


function generate_index(batches,batchSize, data,classes,shuffle)
    local r = {}     
    local rows = data:rows() 
    for i=0,batchSize-1 do         
        r[i] = {}
        r[i][1] = se.FloatMatrix(data:row(i))
        r[i][2] = se.FloatMatrix(classes:row(i))                
    end 
    if(shuffle) then 
        r = rng.shuffle(r)
    end 
    return r 
end 

function generate_batch(numBatches,batchSize, data,classes,shuffle)
    local r = {} 
    local ct = 1
    local rows = data:rows()
    local rc = 0
    for i=0,numBatches-1 do 
        local l = {}   
        local curBatchSize = batchSize      
        if(i == numBatches) then            
            if(data:rows() % batchSize ~= 0) then 
                curBatchSize = data:rows() % batchSize 
            end
        end                                
        for j=0,curBatchSize-1 do                         
            l[j] = {}                    
            local example = data:row(rc)
            local class   = classes:row(rc)
            l[j][1] = se.FloatMatrix(example)
            l[j][2] = se.FloatMatrix(class)
            rc = rc + 1 
            rc = rc % data:rows()                
        end         
        r[i] = l
    end 
    if(shuffle) then 
        for i=1,#r do 
            r[i] = rng.shuffle(r[i])
        end
    end 
    return r 
end 


function batchGradientDescent(network,data,classes,lossFunction, batchSize, learningRate, searchTime, regularizationStrength, momentumFactor,maxIters,shuffle,verbose)    
    assert(network.layers[0].size == data:cols());
    assert(data.rows == classes.rows);
    assert(network.layers[network.numLayers-1].size == classes:cols());
    assert(batchSize <= data:rows());
    assert(maxIters >= 1);

    local errori = {}
    local dWi = {} 
    local dbi = {} 
    local regi = {} 
    local beforeOutputT = createMatrixZeros(network.layers[network.numLayers-2].size,1)
    
    for i=0,network.numConnections-1 do         
        errori[i] = createMatrixZeros(1,network.layers[i].size)                
        dWi[i] = createMatrixZeros(network.connections[i].weights:rows(), network.connections[i].weights:cols())
        dbi[i] = createMatrixZeros(1,network.connections[i].bias:cols())
        regi[i] = createMatrixZeros(network.connections[i].weights:rows(), network.connections[i].weights:cols())
    end 
    errori[network.numConnections] = createMatrixZeros(1, network.layers[network.numConnections].size);
    
    local numHidden = network.numLayers - 2 
    local wTi, errorLastTi,fprimei,inputTi
    
    wTi = {} 
    errorLastTi = {} 
    fprimei = {} 
    inputTi = {} 
    for k=0,numHidden-1 do 
        wTi[k] = createMatrixZeros(network.connections[k+1].weights:cols(),network.connections[k+1].weights:rows())
        errorLastTi[k] = createMatrixZeros(1,wTi[k]:cols())
        fprimei[k] = createMatrixZeros(1,network.connections[k].to.size)
        inputTi[k] = createMatrixZeros(network.connections[k].from.size,1)
    end
    
    local dWi_avg = {} 
    local dbi_avg = {} 
    local dWi_last = {} 
    local dbi_last = {}
    for i=0,network.numConnections-1 do 
        dWi_avg[i] = createMatrixZeros(network.connections[i].weights:rows(),network.connections[i].weights:cols())
        dbi_avg[i] = createMatrixZeros(1,network.connections[i].bias:cols())
        dWi_last[i] = createMatrixZeros(network.connections[i].weights:rows(),network.connections[i].weights:cols())
        dbi_last[i] = createMatrixZeros(1,network.connections[i].bias:cols())
    end 
    local numBatches = math.ceil(data:rows() / batchSize)
    if(data:rows() % batchSize ~=0) then numBatches = numBatches + 1  end 
    local training,batch,epoch,layer 
    local dataBatches = {} 
    local classBatches= {}

    epoch = 0
    batch = 0
        
    local index_list = generate_batch(numBatches,batchSize,data,classes,shuffle)
    while(epoch <= maxIters) do                         
        if(shuffle) then 
            for i=1,#index_list do 
                index_list[i] = rng.shuffle(index_list[i])
            end
        end 
        for batch=0,numBatches-1 do 
            local curBatchSize = batchSize            
            epoch = epoch + 1
            if(epoch > maxIters) then goto label_done end
            if(batch == numBatches) then
                if(data:rows() % batchSize ~= 0) then 
                    curBatchSize = data:rows() % batchSize 
                end
            end                                
            for training=0,curBatchSize-1 do                            
                local x = index_list[batch][training]                
                local example = x[1]
                local target  = x[2]         
                network:forwardPass(example)                      
                
                for layer=network.numLayers-1,1,-1 do                                         
                    local to = network.layers[layer]                                        
                    local con = network.connections[layer-1]
                    if(layer == network.numLayers-1) then                                                                                                                          
                        errori[layer] = to.input - target                                                                           
                        con.from.input:transposeInto(beforeOutputT)                                                                                                    
                        dWi[layer-1] = beforeOutputT * errori[layer]       
                        dbi[layer-1] = errori[layer]:eval()                                                                           
                    else 
                        local hiddenLayer = layer-1                                                          
                        network.connections[layer].weights:transposeInto(wTi[hiddenLayer])
                        errorLastTi[hiddenLayer] =  errori[layer+1]*wTi[hiddenLayer]                        
                        fprimei[hiddenLayer] = con.to.input:eval()                                                                        
                        local d = activationDerivative(con.to.activation)                                                    
                        d(fprimei[hiddenLayer])                                                                             
                        errori[layer] = se.hadamard_float(errorLastTi[hiddenLayer],fprimei[hiddenLayer])                                                                    
                        con.from.input:transposeInto(inputTi[hiddenLayer])
                        dWi[hiddenLayer] = inputTi[hiddenLayer] * errori[layer]                        
                        dbi[hiddenLayer] = errori[layer]:eval()                                               
                    end                     
                end                 
                for i=0,network.numConnections-1 do                     
                    dWi_avg[i] = dWi[i] + dWi_avg[i] 
                    dbi_avg[i] = dbi[i] + dbi_avg[i]                                                                                                                         
                end                                                                         
            end                        
            local currentLearningRate = learningRate 
            if(searchTime ~= 0) then 
                currentLearningRate = learningRate / (1 + (epoch / searchTime))
            end             
            local clr = currentLearningRate/data:rows()
            for i=0,network.numConnections-1 do 
                dWi_avg[i] = dWi_avg[i] * clr 
                dbi_avg[i] = dbi_avg[i] * clr             
                regi[i] = network.connections[i].weights*regularizationStrength
                dWi_avg[i] = regi[i] + dWi_avg[i]            
                dWi_last[i] = dWi_last[i] * momentumFactor               
                dbi_last[i] = dbi_last[i] * momentumFactor               
                dWi_avg[i] = (dWi_last[i] + dWi_avg[i])*-1
                dbi_avg[i] = (dbi_last[i] + dbi_avg[i])*-1
                network.connections[i].weights = dWi_avg[i] + network.connections[i].weights
                network.connections[i].bias = dbi_avg[i] + network.connections[i].bias                
                dWi_last[i] = dWi_avg[i] * -1
                dbi_last[i] = dbi_avg[i] * -1
                dWi_avg[i]:zero() 
                dbi_avg[i]:zero() 
                regi[i]:zero()
            end                                     
            
            if(verbose == true) then 
                if(epoch % 250 == 0 or epoch ==1 ) then 
                    network:forwardPassDataSet(data)
                    if(lossFunction == CROSS_ENTROPY_LOSS) then 
                        print("EPOCH: ", epoch, " loss is ",  network:crossEntropyLoss(network:getOutput(),classes,regularizationStrength))
                    else 
                        print("EPOCH: ", epoch, " loss is ",  network:meanSquaredError(network:getOutput(),classes,regularizationStrength))
                    end
                end 
            end
        end
    end 
    ::label_done::    
    network:forwardPassDataSet(data)
    if(lossFunction == CROSS_ENTROPY_LOSS) then 
        print("final loss is ",  network:crossEntropyLoss(network:getOutput(),classes,regularizationStrength))
    else 
        print("final loss is ",  network:meanSquaredError(network:getOutput(),classes,regularizationStrength))
    end
end    


function XOR(f)
    examples = { {0,0},{0,1},{1,0},{1,1}}
    training = { {0},{1},{1},{0}}

    examples_bp = { {-1,-1},{-1,1},{1,-1},{1,1}}
    training_bp = { {-1},{1},{1},{-1}}

    e = matrix_new(4,2,examples,true)
    t = matrix_new(4,1,training,true)
    net = network_new(2,1,{16},{f},1,f)
    p = parameterset_new(net,e,t,MEAN_SQUARED_ERROR,1,0.1,0,1e-6,0.9,1000,false,true)
    p.verbose=true
    p.batchSize=4
    p.searchTime = 0 
    p.regularizationStrength = 0.0001
    p.learningRate = 0.1
    p.momentumFactor = 0.9
    p.shuffle = true
    print("CraniumSE Online")
    if(use_profiler == true) then profiler.start() end 
    net:optimize(p)
    print("READY")
    if(use_profile == true) then 
        profiler.stop() 
        profiler.report("cranium5.log")
    end
    net:forwardPassDataSet(e)
    output = net:getOutput()
    print(output:rows(),output:cols())
    for i=0,3 do
        if( output:get(i,0) > 0.5 ) then
            io.write(' 1.0, ')
        else
            io.write(' 0.0, ')
        end
    end
    print()
    for i=0,3 do    
            io.write( tostring(output:get(i,0)) .. ',')    
    end
    print()
    print("%accuracy=",net:accuracy(e,t))
    print()
end

function test()
    a = se.FloatMatrix(3,3)
    b = se.FloatMatrix(3,3)
    a:fill(-1)
    softmax(a)
    a:print()
end



--test()
XOR(sigmoid)
XOR(tanH)
XOR(relu)
os.exit()

require('cranium')
rand=require('random')
JSON=require('JSON')
csv = require('csv')

function CreateNetwork(num_features, num_hidden, hidden_data, num_outputs, activate)
    local hidden = cranium.hiddensize_new(num_hidden)
    for i=0,num_hidden-1 do
        cranium.hiddensize_set(hidden,i,hidden_data[i+1].size)
    end
    local activations = cranium.activate_new(num_hidden)
    for i=0,num_hidden-1 do
        cranium.activate_set(activations,i,hidden_data[i+1].activation)
    end
    network = cranium.createNetwork(num_features, num_hidden, hidden, activations, num_outputs, activate)
    return network
end

function NewData(rows, cols, data)
    local example = cranium.data_new(rows,cols)       
    for i=1,rows do        
        local row = cranium.data_get_row(example,i-1)
        for j=1,cols do               
            cranium.data_set_col(row,j-1,data[i][j])
        end
    end
    return example
end



function LoadParameters(json)
    local data = JSON:decode(json)
    local params = cranium.ParameterSet()    
    params.lossFunction = data.lossFunction
    params.batchSize    = data.batchSize
    params.learningRate = data.learningRate
    params.searchTime   = data.searchTime
    params.regularizationStrength = data.regularizationStrength
    params.momentumFactor = data.momentumFactor
    params.maxIters = data.maxIters
    params.shuffle  = data.shuffle
    params.verbose  = data.verbose
    return params,data
end

function Matrix_New(mat) --rows,cols)
    local m = {} 
    --local data = cranium.float_buffer(rows*cols)
    --m.data = data
    --m.matrix = cranium.createMatrix(rows,cols,data) 
    m.matrix = mat

    m.print = function(m) cranium.printMatrix(m.matrix) end 
    m.zero  = function(m) cranium.zeroMatrix(m.matrix) end 
    m.transpose = function(m) return cranium.transpose(m.matrix) end 
    m.transposeInto = function(m,ot) cranium.transposeInto(m.matrix,ot.matrix) end 
    m.add = function(m,b) return cranium.add(m.matrix,b.matrix) end 
    m.addTo = function(m,b) cranium.addTo(m.matrix,b.matrix) end 
    m.addToEachRow = function(m,b) return cranium.addToEachRow(m.matrix,b.matrix) end 
    m.scalarMultiply = function(m,c) cranium.scalarMultiply(m.matrix,c) end 
    m.multiply = function(m,b) return cranium.multiply(m.matrix,b.matrix) end 
    m.multiplyInto = function(m,b,into) cranium.multiplyInto(m.matrix,b.matrix,into.matrix) end 
    m.hadamard = function(m,b) return cranium.hadamard(m.matrix,b.matrix) end 
    m.hadamardInto = function(m,b,into) cranium.hadamardInto(m.matrix,b.matrix,into.matrix) end 
    m.copy = function(m) return cranium.copy(m.matrix) end 
    m.equals = function(m,b) return cranium.eqals(m.matrix,b.matrix) end 
    
    m.get_index = function(m,index) return cranium.matrix_get_index(m.matrix,index) end     
    m.set_index = function(m,index, val) cranium.matrix_set_index(m.matrix, index,val) end 
    m.get = function(m,row,col) return cranium.matrix_get(m.matrix,row,col) end 
    m.set = function(m,row,col,val) cranium.matrix_set(m.matrix,row,col,value) end 

    local meta = {} 
    
    meta.__gc = function(m) 
        --cranium.float_free(m.data)
        cranium.destroyMatrix(m.matrix)
    end 
    setmetatable(m,meta)
    return m 
end 

function Network_New(num_inputs, num_hidden, hidden,num_outputs, activate_out)
    local net = {} 
    net.network = CreateNetwork(num_inputs, num_hidden, hidden, num_outputs, activate_out)
    net.forwardPass = function(n,m) cranium.forwardPass(n.network, m.matrix) end 
    net.forwardPassDataSet = function(n,input) cranium.forwardPassDataSet(n.network, input.ds) end 
    net.crossEntropyLoss = function(net,prediction,actual, reg) return cranium.crossEntropyLoss(net.network, prediction.matrix, actual.ds, reg) end 
    net.meanSquareError = function(net,prediction,actual, reg) return cranium.meanSquareError(net.network, prediction.matrix, actual.ds, reg) end 
    net.getOutput = function(net) return Matrix_New(cranium.getOuput(net.network)) end 
    net.predict = function(net) return cranium.predict(net.network) end 
    net.accuracy = function(net,data,classes) return cranium.accuracy(net.network,data.ds, classes.ds) end 
    net.saveNetwork = function(net,path) cranium.saveNetwork(net.network,path) end 
    net.readNetwork = function(path) net.network = cranium.readNetwork(path) end
    net.optimize = function(n,p) cranium.optimize(p) end
    local meta = {}
    meta.__gc = function(net) cranium.destroyNetwork(net.network) end 
    setmetatable(net,meta)
    return net 
end 


   

function DataSet_New(rows,cols,data)
    local ds = {}
    ds.data = NewData(rows,cols,data)
    ds.ds   = cranium.createDataSet(rows,cols,ds.data)
    ds.get = function(ds,row,col) return cranium.dataset_get(ds.ds,row,col) end 
    ds.set = function(ds,row,col,v) cranium.dataset_set(ds.ds,row,col,v) end 
    local meta = {}         
    meta.__gc = function(d) 
        cranium.destroyDataSet(d.ds)
    end 
    setmetatable(ds,meta)
    return ds 
end

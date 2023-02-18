dnn=require('minidnn')
require('eigen')

math.randomseed(123)
x = eigen.MatrixXd(4,2)

data = {{0,0},{0,1},{1,0},{1,1}}
for i=1,4 do 
    for j=1,2 do 
        x:set(i-1,j-1,data[i][j])
    end 
end 
train = {{0},{1},{1},{0}}
y = eigen.MatrixXd(4,1)
for i=1,4 do 
    y:set(i-1,0,train[i][1])
end 


net = dnn.Network() 

layer1 = dnn.FullyConnectedLayer(dnn.IDENTITY,2,16)
layer2 = dnn.FullyConnectedLayer(dnn.SIGMOID,16,16)
layer3 = dnn.FullyConnectedLayer(dnn.IDENTITY,16,1)

net:add_layer(layer1)
net:add_layer(layer2)
net:add_layer(layer3)

net:set_output(dnn.RegressionMSEOutput())
opt = dnn.SGDOptimizer(0.1)

net:init(0,0.9)

x:transposeInPlace()
y:transposeInPlace()

net:fit(opt,x,y,1,100)

r=net:predict(x)

for i=1,4 do 
    print(r:get(0,i-1))
end


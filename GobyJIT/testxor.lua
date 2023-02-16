require('cranium')
require('eigen')
_examples = {0,0,0,1,1,0,1,1};
_training = {0,1,1,0};
_examples_bp = {-1,-1,-1,1,1,-1,1,1};
_training_bp = {-1,1,1,-1};

examples = cranium.double_vector(8)
for i=1,#_examples do examples[i] = _examples[i] end
training = cranium.double_vector(4)
for i=1,#_training do training[i] = _training[i] end

e = cranium.matrix_new(4,2,examples);
t = cranium.matrix_new(4,1,training);

net = cranium.NeuralNetwork()
input = cranium.Layer(cranium.INPUT,2,cranium.LINEAR);
hidden= cranium.Layer(cranium.HIDDEN,16,cranium.TANH);
output= cranium.Layer(cranium.OUTPUT,1,cranium.LINEAR);

net:addLayer(input);
net:addLayer(hidden);
net:addLayer(output);
net:connect();

p   = cranium.ParameterSet(e,t,1000,4);

p.learning_rate = 0.1;
p.momentum_factor = 0.9;
p.regularization_strength = 0 --1e-6;
p.search_time=1000;
p.verbose = true;
p.shuffle = true;
p.optimizer = cranium.GD_OPTIMIZER;
net:train(p,cranium.STOCHASTIC);
    
net:ForwardPass(e);
output = net:GetOutput();
for i=0,output:rows()-1 do
    for j=0,output:cols()-1  do 
        io.write(output:get(i,j)..',') 
    end
end

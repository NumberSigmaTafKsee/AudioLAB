-- run with octopus
require("analog_ladder_filter_2")
require('stdsamples')
require('plot')

samples = 64
f = analog_ladder_filter_2.LadderFilter(2)
f:setSampleRate(44100.0)
f:setCutoffFreq(5000.0)
f:setResonance(0.0)

a = stdsamples.sample_vector(samples)
b = stdsamples.complex_vector(samples)
r = stdsamples.sample_vector(samples)
gains = stdsamples.sample_vector(11025)
phases= stdsamples.sample_vector(11025)

function mean(x)
    local r = 0
    for i=1,x:size() do
        r = r + x[i]
    end
    return r/x:size()
end
gains:fill(0)
r:fill(0)
fft = stdsamples.FFTPlanRealDouble(samples);
-- i dont know if this is right but it's in the ballpark
for i=1,11025 do
    v = stdsamples.generate_sin(i,44100.0,samples)    
    for j=1,samples do
        a[j] = f:process(v[j])
    end                
    stdsamples.fft(fft,a:data(),b:data())
    
    local mean = 0        
    for j=1,samples/2 do    
         mean = mean + stdsamples.cabsf(b[j])
    end    
    mean = mean / (b:size()/2)
    gains[i] = math.pow(10,mean/20)
end
print('go')
p = plot.Plot_Double()
p:setstyle("lines")
p:plot_x(gains:data(),gains:size(),"wtf?")
os.execute("sleep 15;")
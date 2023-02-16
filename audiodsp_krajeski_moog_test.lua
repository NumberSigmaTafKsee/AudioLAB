-- run with octopus
require("analog_krajeski_moog")
require('stdsamples')
require('plot')

samples = 256
f = analog_krajeski_moog.MoogLadder(44100.0)
f:setFrequency(1000.0)
f:setQ(0.5)

a = stdsamples.sample_vector(samples)
b = stdsamples.complex_vector(samples)
r = stdsamples.sample_vector(samples)
bumper = 100
s = 11025/bumper
gains = stdsamples.sample_vector(s)
phases= stdsamples.sample_vector(s)
w = stdsamples.window(stdsamples.window.HAMMING,samples)

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

-- dont know why it ripples
for i=1,11025 do
    v = stdsamples.generate_sin(i,44100.0,samples)
    for j=1,samples do
        a[j] = 2*f:process(v[j]) * w[j]
    end                

    stdsamples.fft(fft,a:data(),b:data())
    
    local mean = 0        
    for j=2,samples/2-2 do    
         mean = mean + stdsamples.cabsf(b[j])
    end    
    mean = mean / samples
    gains[(i / bumper) + 1] =mean/2 + gains[(i / bumper) + 1]
end
x = stdsamples.sample_vector(s)
for i=1,s do
    x[i] = i * (22050/s)
    gains[i] = 20*math.log10(gains[i])
end
p = plot.Plot_Double()
--p:setstyle("lines")
p:plot_xy(x:data(),gains:data(),gains:size(),"wtf?")
os.execute("sleep 15;")
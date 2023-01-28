-- run with octopus
require("analog_ladder_filter")
require('stdsamples')
require('plot')

samples = 256
f = analog_ladder_filter.LadderFilter(2)
f:setSampleRate(44100.0)
f:setCutoffFreq(5000.0)
f:setResonance(0.0)

a = stdsamples.sample_vector(samples)
b = stdsamples.complex_vector(samples)
r = stdsamples.sample_vector(samples/2)
r:fill(0)
fft = stdsamples.FFTPlanRealDouble(samples);

-- i dont know if this is right but it's in the ballpark
for i=10,22050,25 do
    v = stdsamples.generate_cos(i,44100.0,samples)    
    for j=1,samples do
        a[j] = f:process(v[j])
    end                
    stdsamples.fft(fft,a:data(),b:data())
    for j=1,samples do
        b[j]:real(b[j]:real()/samples)
        b[j]:imag(b[j]:imag()/samples)
    end                
    for j=1,samples/2 do    
        r[j] = r[j] + stdsamples.cabsf(b[j])
    end    
end
--[[    
for i=1,samples do
    b[i] = stdsamples.complex(math.cos(2*math.pi*i*samples/44100),-1*math.sin(2*math.pi*i*samples/44100))
end
stdsamples.ifft(fft,b:data(),a:data())
for i=1,samples/2 do
    r[i] = math.sqrt(a[i]*a[i])
end
]]
x = stdsamples.sample_vector(samples/2)
for i=1,samples/2 do
    x[i] = math.log((i-1)*44100.0/256.0)
end
p = plot.Plot_Double()
p:setstyle("lines")
p:plot_xy(x:data(),r:data(),r:size(),"wtf?")
os.execute("sleep 15;")
require('audio_function_generator')
require('stdsamples')
require('plot')

print(1)
samples = 256
p = plot.Plot_Double()
s = audio_function_generator.SineGenerator(220,44100)
v = stdsamples.sample_vector(samples)
--s.polarity=audio_function_generator.Function.NEGATIVE
for i=1,samples do
    v[i] = s:Tick()
end
p:setstyle("lines")
p:plot_x(v:data(),v:size(),"Sin")
s = audio_function_generator.CosGenerator(220,44100)
for i=1,samples do
    v[i] = s:Tick()
end
p:plot_x(v:data(),v:size(),"Cos")
s = audio_function_generator.PhasorGenerator(220,44100)
for i=1,samples do
    v[i] = s:Tick()
end
p:plot_x(v:data(),v:size(),"Phasor")
s = audio_function_generator.SquareGenerator(220,44100)
for i=1,samples do
    v[i] = s:Tick()
end
p:plot_x(v:data(),v:size(),"Saw")
s = audio_function_generator.CosGenerator(220,44100)
for i=1,samples do
    v[i] = s:Tick()
end
p:plot_x(v:data(),v:size(),"Saw")
s = audio_function_generator.TriangleGenerator(220,44100)
for i=1,samples do
    v[i] = s:Tick()
end
p:plot_x(v:data(),v:size(),"Triangle")
os.execute('sleep 15;')

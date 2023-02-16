require('luapa')
require('luastk')

inputs = {}
for i = 1,4 do 
    inputs[i] = luastk.SineWave()
    inputs[i]:setFrequency(220.0 * i)
end

function noise(input,output,frames)
    for i = 0,2*frames-1,2 do
        local x = 0
        for j = 1,4 do 
            x = x + inputs[j]:tick()
        end 
        x = x / 4.0
        luapa.float_set(output,i,x)
        luapa.float_set(output,i+1,x)
    end
    
end 

luapa.set_audio_func(noise)
luapa.InitAudio(44100,64)
luapa.RunAudio()
luapa.StopAudio()
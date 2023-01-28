require('luapa')

function noise(input,output,frames)
    for i = 0,2*frames-1,2 do
        luapa.float_set(output,i,math.random())
        luapa.float_set(output,i+1,math.random())
    end
    
end 

luapa.set_audio_func(noise)
luapa.InitAudio(44100,64)
luapa.RunAudio()
luapa.StopAudio()
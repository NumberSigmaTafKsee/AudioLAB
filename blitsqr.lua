require('luapa')
require('luastk')

osc = luastk.BlitSquare()
adsr = luastk.ADSR()
adsr:setAllTimes(10.0,5.0,0.5,10.0)
adsr:keyOn()
function noise(input,output,frames)
    for i = 0,2*frames-1,2 do
        local x = osc:tick() * adsr:tick()       
        luapa.float_set(output,i,x)
        luapa.float_set(output,i+1,x)
    end
    
end 

luapa.set_audio_func(noise)
luapa.InitAudio(44100,64)
luapa.RunAudio()
luapa.StopAudio()
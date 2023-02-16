require('luapa')
require('LuaPolyBLEP')
require('luamoog')

osc1 = LuaPolyBLEP.PolyBLEP(44100)
osc1:setWaveform(LuaPolyBLEP.PolyBLEP.SAWTOOTH)
osc2 = LuaPolyBLEP.PolyBLEP(44100)
filt = luamoog.StilsonMoog(44100)
cut = 700.0
filt:SetCutoff(cut)
b = luapa.float_new(1)
function noise(input,output,frames)
    for i = 0,2*frames-1,2 do
        local x = osc1:getAndInc()         
        x=filt:tick(x)        
        luapa.float_set(output,i,x)
        luapa.float_set(output,i+1,x)        
    end
    
end 

luapa.set_audio_func(noise)
luapa.InitAudio(44100,64)
luapa.RunAudio()
luapa.StopAudio()
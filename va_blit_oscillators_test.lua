require('audiodsp_hirestimer')
require('va_blit_oscillators')
require('stdsamples')

timer = audiodsp_hirestimer.HiResTimer()
osc   = va_blit_oscillators.blitSaw()
v     = stdsamples.sample_vector(256)
timer:Start()
    for j=1,10 do
    for i=1,256 do
        v[i] = osc:Tick()
    end
end
time1 = timer:Stop()
timer:Start()
for j=1,10 do
    osc:ProcessBlock(256,v:data(),v:data())
end
time2 = timer:Stop()
timer:Start()
for j=1,10 do
    osc:ProcessSIMD(256,v:data())
end
time3 = timer:Stop()
print("Time1=",time1/10)
print("Time2=",time2/10)
print("Time3=",time3/10)

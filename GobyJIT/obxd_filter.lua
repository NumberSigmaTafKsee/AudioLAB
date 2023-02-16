

function ObxFilter_new()
    local Filter = {} 
    Filter.selfOscPush = false;
    Filter.bandPassSw  = false;
    Filter.mm = 0
    Filter.mmch = 0
    Filter.mmt = 0
    Filter.s1 = 0
    Filter.s2 = 0
    Filter.s3 = 0
    Filter.s4 = 0
    Filter.SampleRate = 44000
    Filter.sampleRateInv = 1.0 / 44000.0
    Filter.rcor = 500.0 / 44000.0
    Filter.rcordInv = 1 / Filter.rcor
    Filter.rcor24 = 970 / 44000.0 
    Filter.rcor24Inv = 1.0 / Filter.rcor24 
    Filter.R = 1.0
    Filter.R24 = 0.0 
    return Filter
end 

function ObxFilter_setMultimode(filter,m)
    filter.mm = m
    filter.mmch = math.floor(filter.mm*3)
    filter.mmt  = filter.mm*3-filter.mmch 
end 

function ObxFilter_setResonance(filter, q)
    filter.R = 1.0 - q 
    filter.R24= (3.5 * q)
end 

function ObxFilter_setSampleRate(filter, sr)
    filter.SampleRate = sr 
    filter.sampleRateInv = 1/sr 
    local rcrate = math.sqrt(44000/sr)
    filter.rcor = (500.0 / 44000) * rcrate 
    filter.rcor24 = (970.0 / 44000) * rcrate 
    filter.rcorInv = 1 / filter.rcor 
    filter.rcor24Inv = 1 / filter.rcor24 
end 

function ObxFilter_diodePairResistanceApprox(x)
    return (((((0.0103592)*x + 0.00920833)*x + 0.185)*x + 0.05 )*x + 1.0)
end

function ObxFilter_NR(filter,sample,g)
	local tCfb;
	if(filter.selfOscPush == false) then
		tCfb = ObxFilter_diodePairResistanceApprox(filter.s1*0.0876) - 1.0;
	else
		tCfb = ObxFilter_diodePairResistanceApprox(filter.s1*0.0876) - 1.035;
    end      
    return ((sample - 2*(filter.s1*(filter.R+tCfb)) - g*filter.s1  - filter.s2)/(1+ g*(2*(filter.R+tCfb)+ g)))
end

function ObxFilter_Apply(filter,sample,g)
    local gpw = math.tan(g * filter.sampleRateInv * math.pi)
    g = gpw 
    local v = ObxFilter_NR(filter,sample,g)    
    local y1= v*g + filter.s1
    filter.s1 = v*g + y1 
    local y2 = y1*g + filter.s2
    filter.s2 = y1*g + y2
    local mc 
    if(filter.bandPassSw == false) then 
        mc = (1-filter.mm)*y2 + filter.mm*v 
    else 
        if( filter.mm < 0.5) then mc = 2*((0.5-filter.mm)*y2 + filter.mm*y1)
        else mc = 2*((1-filter.mm)*y1 + (filter.mm-0.5)*v) 
        end
    end 
    return mc
end

function ObxFilter_NR24(filter, sample, g, lpc)
    local m1 = 1 / (1+g)
    local S = (lpc*(lpc*(lpc*filter.s1 + filter.s2) + filter.s3) + filter.s4)*m1
    local G = lpc*lpc*lpc*lpc;
    return (sample - filter.R24 * S) / (1 + filter.R24*G)
end 

function tptpc(state,inp,cutoff)
    local v = (inp - state) * cutoff / (1 + cutoff)
    local res = v + state
    return res,res+v
end 

function ObxFilter_Apply4Pole(filter, sample, g)
    local g1 = math.tan(g*filter.sampleRateInv * math.pi)
    g = g1 
    local lpc = g/(1+g)
    local y0  = ObxFilter_NR24(filter,sample,g,lpc)
    local v = (y0 - filter.s1) * lpc 
    local res = v + filter.s1
    filter.s1 = res + v
    filter.s1 = math.atan(filter.s1*filter.rcor24)*filter.rcor24Inv
    local y1 = res 
    local y2,y3,y4
    y2,filter.s2 = tptpc(filter.s2,y1,g)
    y3,filter.s3 = tptpc(filter.s3,y2,g)
    y4,filter.s4 = tptpc(filter.s4,y3,g)
    local mc 
    if(filter.mmch == 0) then 
        mc = ((1-filter.mmt) * y4 + filter.mmt * y3)
    elseif(filter.mmch == 1) then 
        mc = ((1-filter.mmt) * y3 + filter.mmt * y2)
    elseif(filter.mmch == 2) then 
        mc = ((1-filter.mmt) * y2 + filter.mmt * y1)
    elseif(filter.mmch == 3) then 
        mc = y1 
    else
        mc = 0
    end 
    return mc * (1 + filter.R24 * 0.45)
end


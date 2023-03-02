function x = fir1lp(order,fc)
    x = fir1(order,fc)
end

function x = fir1lp_cheby(order,fc,ripple)
    x = fir1(order,fc,chebwin(order+1,ripple))
end

function x = fir1lp_tukey(order,fc)
    x = fir1(order,fc,tukeywin(order+1))
end

function x = fir1lp_kaiser(order,fc,beta)
    x = fir1(order,fc,kaiser(order+1,beta))
end

function x = fir1lp_hann(order,fc,beta)
    x = fir1(order,fc,hann(order+1,beta))
end

function x = fir1hp(order,fc)
    x = fir1(order,fc,'high')
end

function x = fir1hp_cheby(order,fc,ripple)
    x = fir1(order,fc,'high',chebwin(order+1,ripple))
end

function x = fir1hp_tukey(order,fc)
    x = fir1(order,fc,'high',tukeywin(order+1))
end

function x = fir1hp_kaiser(order,fc,beta)
    x = fir1(order,fc,'high',kaiser(order+1,beta))
end

function x = fir1hp_hann(order,fc,beta)
    x = fir1(order,fc,'high',hann(order+1,beta))
end

function x = fir1bp(order,f1,f2)
    x = fir1(order,[f1 f2],'bandpass')
end

function x = fir1bp_cheby(order,f1,f2,ripple)
    x = fir1(order,[f1 f2],'bandpass',chebwin(order+1,ripple))
end

function x = fir1bp_tukey(order,f1,f2)
    x = fir1(order,[f1 f2],'bandpass',tukeywin(order+1))
end

function x = fir1bp_kaiser(order,f1,f2,beta)
    x = fir1(order,[f1 f2],'bandpass',kaiser(order+1,beta))
end

function x = fir1bp_hann(order,f1,f2,beta)
    x = fir1(order,[f1 f2],'bandpass',hann(order+1,beta))
end

function x = fir1bs(order,f1,f2)
    x = fir1(order,[f1 f2],'stop')
end

function x = fir1bs_cheby(order,f1,f2,ripple)
    x = fir1(order,[f1 f2],'stop',chebwin(order+1,ripple))
end

function x = fir1bs_tukey(order,f1,f2)
    x = fir1(order,[f1 f2],'stop',tukeywin(order+1))
end

function x = fir1bs_kaiser(order,f1,f2,beta)
    x = fir1(order,[f1 f2],'stop',kaiser(order+1,beta))
end

function x = fir1bs_hann(order,f1,f2,beta)
    x = fir1(order,[f1 f2],'stop',hann(order+1,beta))
end
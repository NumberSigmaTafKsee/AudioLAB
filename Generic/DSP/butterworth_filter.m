function [sos,g]=butterworth_filter( order, wc )
    [Z,P,K] = butter(order,wc,'s')
    [sos,g] = zp2sos(Z,P,K)       
end
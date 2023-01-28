function [sos,g]=butterworthlpf( order, wc )
    [Z,P,K] = butter(order,wc,'s')
    [sos,g] = zp2sos(Z,P,K)
end

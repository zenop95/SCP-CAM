function [out] = prcLimit(data,prc_up,prc_lo)
    lim_up = prctile(data,prc_up);
    lim_lo = prctile(data,prc_lo);
    data1  = data(data<=lim_up);
    out    = data1(data1>=lim_lo);
end
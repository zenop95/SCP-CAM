function dcm6 = rot6(dcm,w)
dcm6 = [dcm zeros(3); -skew(w) dcm];
end
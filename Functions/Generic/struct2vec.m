function vec = struct2vec(struc,fieldName)
n   = size(struc,2);
vec = nan(n,1);
for i = 1:n
    vec(i) = struc(i).(fieldName);
end
end
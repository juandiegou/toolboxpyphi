function all_bipartitions = all_possible_bipartitions(n)

%-----------------------------------------------------------------------
% This function gives you all possible bipartitions of a system with n
% nodes. Each row in the output is a community assignment vector for one
% possible bipartition. This function can get really slow in systems with 
% more than 16 nodes (because the number of possible bipartitions explodes
% super-exponentially).
%-----------------------------------------------------------------------

all_bipartitions = [];
for m=1:2^(n-1)-1 
 
    binstr=dec2bin(m,n);
    
    vec=[];

    for l=1:n
        vec=cat(2,vec, str2num(binstr(l)));
    end
    all_bipartitions(m,:)=vec+1;
end
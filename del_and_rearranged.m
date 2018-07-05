function [rearranged_array] = del_and_rearranged(array, set) 
[m, n] = size(array);
del_array = zeros(m, n);
for i=1:length(set)
    del_array = del_array + (array > set(i));
end
rearranged_array = array - del_array;

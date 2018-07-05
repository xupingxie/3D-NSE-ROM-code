%delte centers and substract 1 from those greater then these centers
function [array1] = del_centers(array, set)
[array1] = del_and_rearranged(array, set);
%substract centers
for i = 1:length(array1)-length(set)-1
    if array1(i)== array1(i+1)
        array1 = [array1(1:i), array1(i+2:end)];
    end
end
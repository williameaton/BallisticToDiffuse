
% Function provided by Surya 
%this function computes the moving average of the given data vector
%input: data: data vector
%          s: scale factor
%output: mmean: moving average         %         scale = 2                scale = 3
function mmean = movingmean(data,s)    %  u1  u2   u3   u4   u5    u6     u1  u2  u3 u4 u5 u6  
n = length(data);                      %  \   /\   /\   /\   / \   /        \ | \/| \/|\/| /
for i = 1:n-s+1                        %   \ /  \ /  \ /  \ /   \ /          \| /\| /\|/\|/
 mmean(i) = mean(data(i:i+s-1));       %    d1   d2   d3   d3   d4            d1  d2  d3 d4
end
end

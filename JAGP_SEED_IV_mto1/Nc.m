function [N_c] = Nc(Y)
%NTC 此处显示有关此函数的摘要
%   此处显示详细说明
% Y:label matrix, n*c
n_c=sum(Y);
N_c=diag(1./n_c);
end


function [N_c] = Nc(Y)
%NTC �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
% Y:label matrix, n*c
n_c=sum(Y);
N_c=diag(1./n_c);
end


function [S]= Initialization_S(Xs,Xt)
%��ʼ����ͼ�ǻ���ԭʼԴ���Ŀ�������ݽ��й���
X=[Xs,Xt];
options = [];
options.NeighborMode = 'KNN';
options.k = 5;
options.WeightMode = 'HeatKernel';
options.t =1;
%��ʼ������S(����ԭʼ��Դ���Ŀ��������)
S = constructW(X',options);
S = full(S);
end
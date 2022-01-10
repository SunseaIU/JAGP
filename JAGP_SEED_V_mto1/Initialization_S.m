function [S]= Initialization_S(Xs,Xt)
%初始化的图是基于原始源域和目标域数据进行构造
X=[Xs,Xt];
options = [];
options.NeighborMode = 'KNN';
options.k = 5;
options.WeightMode = 'HeatKernel';
options.t =1;
%初始化构造S(基于原始的源域和目标域数据)
S = constructW(X',options);
S = full(S);
end
%JAGP_SEED4_many_to_one:  session1_subject1(target data)
clear; 
clc; 
close all;
% Input£º
% Xs: the data of source domain, d*ns;
% Xt: the data of traget domain, d*nt;
% Ys: the label of source domain, ns*c;
% Yt: the label of target domain, nt*c;
% options: the parameter of JAGP model;
% Ouput:
% Ft: soft label of target domain, nt*c;
% W: projection matrix ,d*p;
% S: Structured graph of aligned source and target domain data, n*n;
% load the labeled and unlabeled data
load('session1_sub1.mat');
Xs = X_src; 
Ys = Y_src;
Xt = X_tar;
Yt = Y_tar;
[~,ns]=size(Xs);
[~,nt]=size(Xt);
n=ns+nt;
c=max(Ys);
options.p=32;
Xs=Xs-repmat(mean(Xs,2),[1,ns]);
Xt=Xt-repmat(mean(Xt,2),[1,nt]);

lambdalib=0.00119;
gammalib=0.009;
alphalib=0.00107;
Ft_init=ones(nt,c)./c;
S_init=Initialization_S(Xs,Xt);
for i1 = 1:length(lambdalib)
    options.lambda=lambdalib(i1);
    for i2 = 1 : length(gammalib)
        options.gamma=gammalib(i2);
        for i3 = 1: length(alphalib)
            options.alpha=alphalib(i3);
            [W,S,Ft,acc,obj]=JAGP(Xs,Ys,Xt,Yt,Ft_init,S_init,options);
            fprintf('dim=%d,lambda=%.4f,gamma=%0.4f,alpha=%0.4f,acc=%0.4f \n',options.p,options.lambda,options.gamma,options.alpha,acc);
            if i1==1&&i2==1&&i3==1
                acc_best=acc;
                Obj=obj;
                Ft_pre=Ft;
                S_pre=S;
                W_pre=W;
                option_best=options;
            elseif acc_best<acc
                acc_best=acc;
                Obj=obj;
                Ft_pre=Ft;
                S_pre=S;
                W_pre=W;
                option_best=options;
            end
        end
    end
end
fprintf('The best acc=%0.4f\n',acc_best);


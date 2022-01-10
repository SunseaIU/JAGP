function [W_pre,S_pre,Ft_pre,acc_best,obj]=JAGP(Xs,Ys,Xt,Yt,Ft_init,S_init,options)
p = options.p;
alpha = options.alpha;
lambda=options.lambda;
gamma=options.gamma;
[d,ns]=size(Xs);
[~,nt]=size(Xt);
c=max(Ys);
n=ns+nt;
X=[Xs,Xt]; 
Ys=onehot(Ys,c);

Ns=Nc(Ys);
Nss=diag([1/ns,sum(Ns)]);%ns*c
Yss=[ones(ns,1),Ys];%ns*(c+1)
H=eye(n)-1/n * ones(n,n);
Ft=Ft_init;
S= S_init;
D=diag(sum(S,2));
L=D-S;
Nt=Nc(Ft);
Ntt=diag([1/nt,sum(Nt)]);%(c+1)*(c+1)
Ftt=[ones(nt,1),Ft];%nt*(c+1)
A=Xs*Yss*Nss-Xt*Ftt*Ntt;
[W, ~, ~]=eigs(A*A'+2*alpha*X*L*X'+0.1*eye(d),X*H*X', p, 'SM');
F=[Ys;Ft];
NITER=30;
obj=[];
for iter=1:NITER
    dist_wx=L2_distance_1(W'*X,W'*X); 
    dist_f=L2_distance_1(F',F'); 
    dw=alpha*dist_wx+0.5*lambda*dist_f; 
    % update S
    for i=1:1:n
        m=-dw(i,:)./(2*gamma);
        [S(:,i),~]=EProjSimplex_new(m);
    end
    clear i;
    S=(S+S')/2;
    D=diag(sum(S,2));
    %Laplacian matrix
    L=D-S;
    Lst=L(1:ns,ns+1:n);
    Dtt=D(ns+1:n,ns+1:n); 
    Stt=S(ns+1:n,ns+1:n);
    N=Nt*Nt'; 
    Z=Xt'*(W*W')*Xt;
    M=lambda*Ys'*Lst-Nt*Ns'*Ys'*Xs'*(W*W')*Xt;
    % update Ft
    for i=1:1:nt
        bb=zeros(c,1);
        for j=1:1:nt
            fj=Ft(j,:)';
            bb=bb+lambda*fj*Stt(i,j)-Z(i,j)*N*fj;
        end
        b=bb-2*M(:,i);
        mm=b./(lambda*Dtt(i,i));
        [v,~] = EProjSimplex_new(mm);
        Ft(i,:)=v';
    end
    % update W
    F=[Ys;Ft];
    Nt=Nc(Ft);
    for tt=1:c
        if(Nt(tt,tt)==inf)
            Nt(tt,tt)=1e-10;
        end
    end
    Ntt=diag([1/nt,sum(Nt)]);
    Ftt=[ones(nt,1),Ft];
    A=Xs*Yss*Nss-Xt*Ftt*Ntt;
    [W, ~, ~]=eigs(A*A'+2*alpha*X*L*X'+0.1*eye(d),X*H*X', p, 'SM');
    temp1(iter,1)=norm(W'*A,'fro')^2;
    temp2(iter,1)=trace(F'*L*F);
    temp3(iter,1)=trace(W'*X*L*X'*W);
    temp4(iter,1)=trace(S'*S);
    obj=temp1+lambda*temp2+2*alpha*temp3+gamma*temp4;
    [~,predict_label] = max(Ft,[],2);
    acc = length(find(predict_label == Yt))./length(Yt);
    if iter==1
        acc_best=acc;
        Ft_pre=Ft;
        W_pre=W;
        S_pre=S;
    elseif iter>1 && acc_best<acc
        acc_best=acc;
        Ft_pre=Ft;
        W_pre=W;
        S_pre=S;
    end
    fprintf('dim=%d,lambda=%0.4f,gamma=%0.4f,alpha=%0.4f,The acc=%0.4f\n',p,lambda,gamma,options.alpha,acc);
end
end
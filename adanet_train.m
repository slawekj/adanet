function [index weight] = adanet_train(X,~,sPOSY,T,TA,Psi,Xi,s,delta)
%initialization of parameters
[N,P]=size(X);

% select number of points to sample
S=round(s*N);

weight=zeros(1,P);

for t=1:T
    allI=randperm(N);
    randI=allI(1:S);
    
    % select midpoint
    M=round((Psi+(Xi-Psi)*t/T)*N);
    
    % select buffer size
    B=round(delta*N);
    
    negA=1:M-1-B;
    negI=sPOSY(intersect(negA,randI));
    
    posA=(1:N-M-B)+M+B;
    posI=sPOSY(intersect(posA,randI));
    
    % now let's evaluate each of the features
    % we will use ADABOOST for that reason
    weight = weight + adaboost(X([negI;posI],:),[zeros(size(negI)); ones(size(posI))],TA);
end

weight=weight./T;
[weight, index]=sort(weight,'descend');
end

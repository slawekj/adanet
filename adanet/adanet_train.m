%   This is an implementation of ADANET algorithm for Gene Regulatory Network
%   inference from mRNA expression data, in form of a Matlab package.
%   Copyright (C) 2014  Janusz Slawek
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program, see LICENSE.
%

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
    % we will use ADABOOST for that purpose
    weight = weight + adaboost(X([negI;posI],:),[zeros(size(negI)); ones(size(posI))],TA);
end

weight=weight./T;
[weight, index]=sort(weight,'descend');
end

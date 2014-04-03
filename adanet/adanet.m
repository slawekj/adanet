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

function V = adanet(E,regulatorsI,T,t,Psi,Xi,s,delta)

[eN,eP]=size(E);
E=bsxfun(@minus,E,repmat(mean(E),eN,1));
E=bsxfun(@rdivide,E,std(E));

AllGenes=1:eP;
V=zeros(eP,eP);

for B=AllGenes
    As=setdiff(regulatorsI,B);
    
    T_all=E(As,B);
    F_all=E(As,As);
    [~,sPOS]=sort(T_all);
    
    [index weight]=adanet_train(F_all,T_all,sPOS,T,t,Psi,Xi,s,delta);
    weight(isnan(weight))=0;
    V(As(index),B)=weight;
    
    fprintf(['Neighbors of ' num2str(B) ' done\n']);
end

end

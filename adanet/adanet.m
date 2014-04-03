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

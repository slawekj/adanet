% set default parameters
C=30;
T=10;
Psi=0.25;
Xi=0.75;
S=0.67;
delta=0.05;

% create random data table

path(path,'src');

RandomData = normrnd(0,1,100,100);

V=adanet(RandomData,1:10,C,T,Psi,Xi,S,delta);

% you can apply post-processing steps here

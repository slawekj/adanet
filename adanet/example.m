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

% set default parameters
C=30;
T=10;
Psi=0.25;
Xi=0.75;
S=0.67;
delta=0.05;

path(path,'src');

% create random data table

RandomData = normrnd(0,1,100,100);
TFs        = 1:100;

V=adanet(RandomData,TFs,C,T,Psi,Xi,S,delta);

% you can apply post-processing steps here

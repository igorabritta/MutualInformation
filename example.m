% For example, generate two random distributions and calculate the mutual
%   information between them using the KDE method. For MATLAB R2015b

A=randn(10000,1);
B=randn(10000,1);
[I1,In1]=mi(A,B,'KDE');

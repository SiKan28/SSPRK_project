function [A,b]=shuosher2butcher(alpha,beta);
%function [A,b,c]=shuosher2butcher(alpha,beta);
%function [A,b,c]=shuosher2butcher(lambda,mu);
%
%By David Ketcheson
%
%Generate Butcher form of a Runge-Kutta method,
%given its Shu-Osher or modified Shu-Osher form
%
%For an m-stage method, alpha and beta (or lambda and mu) should be 
%matrices of dimension (m+1) x m 
%
%Note that MATLAB indexes from 1, while the Shu-Osher coefficients
%are usually indexed from zero.

m=size(alpha,2);
X=eye(m)-alpha(1:end-1,:);
A=X\beta(1:end-1,:);
b=beta(end,:)+alpha(end,:)*A;

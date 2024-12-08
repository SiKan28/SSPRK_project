function [alpha,beta]=butcher2shuosher(A,b,r);
%function [alpha,beta]=butcher2shuosher(A,b,r);
%
%By David Ketcheson
%
%Generate Shu-Osher form of an explicit Runge-Kutta method,
%given its Butcher form and radius of absolute monotonicity
%
%For an m-stage method, A should be an m x m matrix
%and b should be a column vector of length m.
%
%Note that MATLAB indexes from 1, while the Shu-Osher coefficients
%are usually indexed from zero.

if nargin<3 r = am_radius(A,b); end

s=size(A,1);
K=[A;b'];
G=eye(s)+r*A;
beta=K/G;
alpha=r*beta;
for i=2:s+1
  alpha(i,1)=1-sum(alpha(i,2:s));
end

%% 
function r = am_radius(A,b)
%function r = am_radius(A,b)
%
%By David Ketcheson
%
%Evaluates the Radius of absolute monotonicity
%of a Runge-Kutta method, given the Butcher array.
%
%For an m-stage method, A should be an m x m matrix
%and b should be a column vector of length m.
%
%Accuracy can be changed by modifying the value of eps.
%Methods with very large radii of a.m. (>50) will require
%rmax to be increased.

rmax=50; eps=1.e-12;

m=length(b); e=ones(m,1);
K=[A;b'];
rlo=0; rhi=rmax;

while rhi-rlo>eps  %use bisection
  r=0.5*(rhi+rlo);
  X=eye(m)+r*A; beta=K/X; ech=r*K*(X\e);
  if (min(beta(:))<-3.e-16 || max(ech(:))>1.+3.e-16)
    rhi=r;
  else
    rlo=r;
  end
end

if rhi==rmax % r>=rmax
  error('Error: increase value of rmax in am_radius.m');
else
  r=rlo;
end

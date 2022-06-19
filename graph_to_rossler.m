% Outputs a multivariate time-series (node-x-timesteps), as well as the
% binarized and symmetrized adjacency matrix used to generate that
% timeseries

function [data,A] = graph_to_rossler(A,T, a, b, c, d, freq)


% -----------------------------------------------------------------------
%   FUNCTION: phi_comp.m
%   PURPOSE:  Generates a multivariate time-series of coupled stochastic
%   Rössler oscillators, from an adjacency matrix. Also outputs a binarized
%   and symmetrized adjacency matrix, if binarization and symmetrization
%   were used.
%
%   INPUT:  A         -     An adjacency matrix. If not binary and symmetric,
%                           will be transformed to be
%
%           T         -     Number of time-points to simulate
%           a,b,c,d   -     Parameters of the model
%           freq      -     Intrinsic frequency of the oscillators


N=size(A,1);

% binarize and symmetrize the adjacency matrix
if ~issymmetric(A)
B=A+A';
B(1:N+1:end)=diag(A);
B(B>0)=1;
A=B;
end
Gmatrix = diag(degrees_und(A))-A; % Get the Laplacian matrix

% set a coupling strength according to the master stability function
gamma = eig(Gmatrix);
bounds = [.186/gamma(2) 4.614/gamma(end)];
lower_bound = min(bounds);
upper_bound = max(bounds);
e = (upper_bound-lower_bound)*.5+lower_bound;

dt = 0.001; % integration step

% run 50 times more time-steps than you set (will be downsampled later),
% and also add 80,000 time-points, which will be discarded (allows the
                                                            % oscillators to settle)
steps =T*50+80000;

w=.1*randn(N,1)+freq; % add some variance to the oscillators' intrinsic frequencies

x=zeros(steps,N);
y=x;
z=x;

% random initial conditions for the oscillators
xold = (30*rand(N,1)-15)/5;
yold = (30*rand(N,1)-15)/5;
zold = (40*rand(N,1)-5)/5;

for t= 1:steps

coupling = -Gmatrix*xold;
dxdt=-w.*yold-zold+e*coupling;
dydt=w.*xold+a*yold;
dzdt=b+zold.*(xold-c);

xnew=xold+dt*dxdt;
ynew=yold+dt*dydt+d*randn(N,1);
znew=zold+dt*dzdt;

xold=xnew;
yold=ynew;
zold=znew;

x(t,:)=xnew;
y(t,:)=ynew;
z(t,:)=znew;
end
y=y(80001:end,:);
y=downsample(y,50)'; % downsample
data=y;

%% Script for closed-loop performance of SafEDMD-based SOS controller
% Inputs: 
%   - none 
%
% Outputs: 
%   - none
%
%
% __author__ = "Robin Straesser"
% __contact__ = "robin.straesser@ist.uni-stuttgart.de"
% __date__ = "2025/01/17"

clear;clc;clearvars;close all;rng(0)
addpath("fcn\")
format long;

%% Define the system
% system dynamics
f   = @(x) [-2*x(1);x(2)-x(1)^2];
g    = @(x) [0;1];
sys.ode = @(x,u) f(x) + g(x)*u;
sys.n = 2;
sys.m = 1;

% parameters for the data collection and controller design
param.xmax =  1; 
param.xmin = -1;
param.d = 200;
degrees.alpha = 1;
degrees.beta = degrees.alpha;

% Lifting function
param.Phi = @(x) [1;x;x(2)-1/5*x(1)^2;x(1)*x(2)];
param.gradhPhi = @(x) [1,0;0,1;-2/5*x(1),1;x(2),x(1)]';

%% Run the data-driven design in continuous time with the SOS controller 
sys.timeVariant = 'continuous-time';

%% Data generation
[X0,X1] = generateData(sys,param);

%% Apply SafEDMD
[param,sys,X,Y] = SafEDMD(X0,X1,sys,param);

%% Controller design
eps.P = 1e-6;
eps.tau = 1e-7;
eps.rho = 1e-6;
eps.lambda = 1e-6;
eps.eta = 1e-6;

z = sdpvar(sys.N,1);
% controller denominator
ud = 1;
% specify performance input and output
param.performance = true;
sys.Bw = eye(sys.n);
sys.p = size(sys.Bw,2);
sys.C = [eye(sys.n),zeros(sys.n,sys.N-sys.n)];
sys.q = size(sys.C,1);
sys.D = zeros(sys.q,sys.m);
sys.tD = zeros(sys.q,sys.m*sys.N);
sys.Dw = zeros(sys.q,sys.p);
sys.Sw = zeros(sys.p,sys.q);
sys.Rw = eye(sys.q);
sys.Bwz = param.gradhPhi(z(1:sys.n))'*sys.Bw;

%% Run the controller design for different cx=cu values to compute L2 gain
cxs = [1e0,7e-1,6e-1,5e-1,4e-1,3e-1,2e-1,1e-1,5e-2,1e-2,5e-3,1e-3,5e-4,1e-4];
gammaMaxs = [150,150,150,80,80,20,20,20,20,10,5,1.5,1.5,1.5,1.5,1.5];
gamma = NaN(size(cxs));
for i = 1:length(cxs)
    param.cx = cxs(i); param.cu = param.cx;

    % bisection to find the smallest L2-gain bound
    fprintf('%i/%i: Optimize minimal gamma value...\n',i,length(cxs))
    gamma_max = gammaMaxs(i);
    gamma_min = 0;
    gamma(i) = 0;
    while gamma_max - gamma_min > 1e-5
        gamma_guess = (gamma_max + gamma_min)/2;
        sys.Qw = -gamma_guess^2*eye(sys.p);
        successful = true;
        fprintf('%i/%i: Current gamma test: %.4f\n',i,length(cxs),gamma_guess)
        try 
            evalc("controllerDesignSOS_continuous(sys,eps,param,ud,degrees,z,'trace(P)',true)");
        catch
            successful = false;
        end
        if successful
            gamma_max = gamma_guess;
            gamma(i) = gamma_guess;
        else
            gamma_min = gamma_guess;
        end
    end
    fprintf('%i/%i: Done. gamma=%.4f\n',i,length(cxs),gamma(i))
end
figure;
set(groot, 'defaultAxesColorOrder', get(gca,'colororder')); % Default color order
set(groot,'defaultAxesFontSize', 14); % Set font size
set(groot,'defaultLineLineWidth', 1.5); % Set line width
loglog(cxs,gamma,'ko')
grid on
xlim([min(cxs),max(cxs)])
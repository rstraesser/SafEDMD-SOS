%% Script for closed-loop stability of SafEDMD-based LMI and SOS controller
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

%% Run the data-driven design 
% for continuous and discrete time with both the LMI and SOS controllers 
timeVariants = {'continuous-time','discrete-time'};
colorsSOS = dictionary(timeVariants{1},'b',timeVariants{2},'c');
colorsLMI = dictionary(timeVariants{1},'r',timeVariants{2},'m');
figure;hold all;grid on;
set(groot, 'defaultAxesColorOrder', get(gca,'colororder')); % Default color order
set(groot,'defaultAxesFontSize', 14); % Set font size
set(groot,'defaultLineLineWidth', 1.5); % Set line width
for timeVariant = timeVariants
    sys.timeVariant = timeVariant{:};
    switch sys.timeVariant 
        case 'continuous-time'
            param.cx = 1e-1; 
            param.cu = 1e-1;
        case 'discrete-time'
            param.cx = 6e-3; 
            param.cu = 6e-3;
            param.DeltaT = 0.01;
    end

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
    switch sys.timeVariant 
        case 'continuous-time'
            % controller denominator
            ud = 1;
            [Kn,PinvSOS] = controllerDesignSOS_continuous(sys,eps,param,ud,degrees,z,'trace(P)',true);
    
        case 'discrete-time'
            % controller denominator
            [~,c,v] = polynomial(z,2*degrees.alpha,2*degrees.alpha);
            ud = 1 + ones(size(c))'*v;
            [Kn,PinvSOS] = controllerDesignSOS_discrete(sys,eps,param,ud,degrees,z,'trace(P)',false);
    end
    Knfunc = @(zvar) replace(Kn,z,zvar);
    if ~isequal(ud,1)
        udfunc = @(zvar) replace(ud,z,zvar);
        ufuncZ = @(z) 1/(udfunc(z))*Knfunc(z)*z;
    else
        udfunc = @(zvar) 1;
        ufuncZ = @(z) Knfunc(z)*z;
    end
    ufunc = @(x) ufuncZ(param.hPhi(x));

    %% Compare with LMI-based controller design
    sys.Pi.Rz = 1e-2;
    sys.Pi.Sz = zeros(sys.N,1);
    sys.Pi.Qz = -eye(sys.N);
    eps.F = 1e-6;
    eps.Lambda = 1e-7;
    eps.nu = 1e-7;
    %
    switch sys.timeVariant 
        case 'continuous-time'
            [KLMI,KwLMI,PinvLMI] = controllerDesignLMI_continuous(sys,eps,param);
    
        case 'discrete-time'
            [KLMI,KwLMI,PinvLMI] = controllerDesignLMI_discrete(sys,eps,param);
    end
    uLMIfunc = @(x) (eye(sys.m)-KwLMI*kron(eye(sys.m),param.hPhi(x))) \ (KLMI*param.hPhi(x));

    %% Plot region of attractions
    % generate grid over sampling region
    xx = param.xmin:0.001:param.xmax;
    yy = param.xmin:0.001:param.xmax;
    samples = [xx                           , xx                           , param.xmax*ones(1,length(yy)), param.xmin*ones(1,length(yy)); 
               param.xmax*ones(1,length(xx)), param.xmin*ones(1,length(xx)), yy                           , yy                           ];
    
    % calculate largest region of attraction for SOS controller via
    % bisection
    fprintf('%s: Calculate largest SOS Lyapunov sublevel set...',sys.timeVariant)
    cmax = 1e12; cmin = 0; cSOS = 0;
    while cmax - cmin > 1e-7
        cguess = (cmax + cmin)/2;
        touched = false;
        for i=1:size(samples,2)
            if param.hPhi(samples(:,i))'*PinvSOS*param.hPhi(samples(:,i)) <= cguess
                touched = true;
                break
            end
        end
        if touched
            cmax = cguess;
            cSOS = cguess;
        else
            cmin = cguess;
        end
    end
    fprintf('Done.\n')
    % calculate largest region of attraction for LMI controller via
    % bisection if LMI controller design was feasible
    if ~isnan(PinvLMI)
        fprintf('%s: Calculate largest LMI Lyapunov sublevel set...',sys.timeVariant)
        Rz_max = 1000; Rz_min = sys.Pi.Rz; Rz_guess = Rz_min;
        while Rz_max - Rz_min > 1e-7
            sys.Pi.Rz = (Rz_max + Rz_min)/2;
            touched = false;
            switch sys.timeVariant 
                case 'continuous-time'
                    [~,~,~,PinvLMI] = evalc('controllerDesignLMI_continuous(sys,eps,param)');
            
                case 'discrete-time'
                    [~,~,~,PinvLMI] = evalc('controllerDesignLMI_discrete(sys,eps,param)');
            end
            for i=1:size(samples,2)
                if param.hPhi(samples(:,i))'*PinvLMI*param.hPhi(samples(:,i)) <= 1
                    touched = true;
                    break
                end
            end
            if touched
                Rz_max = sys.Pi.Rz;
                Rz_guess = sys.Pi.Rz;
            else
                Rz_min = sys.Pi.Rz;
            end
        end
        sys.Pi.Rz = Rz_guess;
        switch sys.timeVariant 
            case 'continuous-time'
                [~,KLMI,KwLMI,PinvLMI] = evalc('controllerDesignLMI_continuous(sys,eps,param)');
        
            case 'discrete-time'
                [~,KLMI,KwLMI,PinvLMI] = evalc('controllerDesignLMI_discrete(sys,eps,param)');
        end
        uLMIfunc = @(x) (eye(sys.m)-KwLMI*kron(eye(sys.m),param.hPhi(x))) \ (KLMI*param.hPhi(x));
        fprintf('Done.\n')
    end

    % plot resulting region of attractions of SOS and LMI controllers by
    % checking which grid points are containted
    [XX,YY] = meshgrid(xx,yy);
    idxSOS = false(size(XX));
    idxLMI = false(size(XX));
    for i=1:size(XX,1)
        for j=1:size(XX,2)
            if param.hPhi([XX(i,j);YY(i,j)])'*PinvSOS*param.hPhi([XX(i,j);YY(i,j)]) <= cSOS
                idxSOS(i,j) = true;
            end
            if param.hPhi([XX(i,j);YY(i,j)])'*PinvLMI*param.hPhi([XX(i,j);YY(i,j)]) <= 1
                idxLMI(i,j) = true;
            end
        end
        fprintf('%s: Checked grid point: %i/%i.\n',sys.timeVariant,i,length(xx))
    end
    xRoASOS = [XX(idxSOS)';YY(idxSOS)'];
    xRoALMI = [XX(idxLMI)';YY(idxLMI)'];
    fprintf('%s: Get boundary points of SOS-RoA.\n',sys.timeVariant)
    idxSOS = boundary(xRoASOS(1,:)',xRoASOS(2,:)');
    fprintf('%s: Get boundary points of LMI-RoA.\n',sys.timeVariant)
    idxLMI = boundary(xRoALMI(1,:)',xRoALMI(2,:)');
    fprintf('%s: Plotting SOS RoA.\n',sys.timeVariant)
    xlim([param.xmin,param.xmax])
    ylim([param.xmin,param.xmax])
    plot(xRoASOS(1,idxSOS),xRoASOS(2,idxSOS),'.','color',colorsSOS(sys.timeVariant))
    fprintf('%s: Plotting LMI RoA.\n',sys.timeVariant)
    plot(xRoALMI(1,idxLMI),xRoALMI(2,idxLMI),'.','color',colorsLMI(sys.timeVariant))
    drawnow
end
function [Kn,Pinv] = controllerDesignSOS_discrete(sys,eps,param,ud,degrees,z,cost,rescale)
%% Function to design a controller for a continuous-time nonlinear system
% Inputs: 
%   - sys: system description
%   - eps: small margin to ensure positive definiteness
%   - param: parameters defining the data generation and area of interest
%   - ud: scalar polynomial denominator
%   - degrees: vector defining alpha and beta for the polynomial degrees
%   - z: substitute for the lifted state
%   - cost: optional, defined a objective for the optimization
%   - rescale: optional, boolean if the SOS program should be rescaled
%
% Outputs: 
%   - Kn: (m x N) polynomial matrix (numerator)
%   - Pinv: (N x N)-dimensional Lyapunov matrix 
%
% __author__ = "Robin Straesser"
% __contact__ = "robin.straesser@ist.uni-stuttgart.de"
% __date__ = "2025/01/17"
if nargin <= 6
    cost = '[]';
    rescale = true;
end
fprintf('Start optimization program...')
% Optimization: Decision variables
lmis = []; params = [];
    % (Inverse) Lyapunov matrix
    P = sdpvar(sys.N,sys.N,'symmetric');
    lmis = [lmis,P-eps.P*eye(size(P))>=0];
    % Controller matrix
    [Ln,Ln_coeff] = polymatrix(z,2*degrees.alpha-1,[sys.m,sys.N],'full');
    params = [params;vertcat(Ln_coeff{:})];

    % Decay factor rho
    rho = sdpvar(1);
    params = [params;rho];
    lmis = [lmis;rho>=eps.rho];

    % controller denominator needs to be SOS: check if ud is SOS
    [~,info_ud,~,Q_ud] = evalc("solvesos(sos(ud))");
    if info_ud.problem ~= 0 || min(eig(Q_ud{:})) < 0
        error("Chosen denominator is NOT SOS!")
    end
    
    if param.cx < 0 || param.cu < 0
        error('The matrix defining the proportional bound needs to be positive!')
    elseif param.cx == 0 && param.cu == 0
        tau = 0;
        F12 = []; F13 = []; F22 = []; F23 = []; F24 = []; F24 = []; F33 = []; F34 = []; F34 = [];
    else
        % multiplier residual
        [tau,tau_coeff] = polynomial(z,2*degrees.beta); % needs to be even degree
        params = [params;tau_coeff];
        lmis = [lmis,tau_coeff(1)>=eps.tau];

        F12 = zeros(sys.N,sys.N);
        F13 = zeros(sys.N,sys.m);
        if rescale
            F22 = tau/2*eye(sys.N);
            F24 = param.cx*P*ud;
            F33 = tau/2*eye(sys.m);
            F34 = param.cu*Ln;
        else
            F22 = tau/2/param.cx^2*eye(sys.N);
            F24 = P*ud;
            F33 = tau/2/param.cu^2*eye(sys.m);
            F34 = Ln;
        end
        F23 = zeros(sys.N,sys.m);

        if param.cx == 0
            F12 = []; F22 = []; F23 = []; F24 = [];
        end
        if param.cu == 0
            F13 = []; F23 = []; F33 = []; F34 = [];
        end
    end

    F11 = ud*P - rho*ud*eye(sys.N) - tau*eye(sys.N);
    F14 = ud*sys.A*P + sys.B0*Ln + sys.tB*kron(Ln,z);
    F44 = ud*P;
    
    F = [F11 ,F12 ,F13 ,F14;
         F12',F22 ,F23 ,F24;
         F13',F23',F33 ,F34;
         F14',F24',F34',F44];
    
    [~,cost] = evalc(cost);
    lmis = [lmis,sos(F)];
    opt = solvesos(lmis,cost,sdpsettings('solver','mosek','sos.clean',1e-12),params);

    switch opt.problem 
        case -1
            fprintf(2,'Optimization status: %s: ILL_POSED\n', opt.info)
            error('Controller design not successful!')
        case 0 
            fprintf('Optimization status: %s: PRIMAL_AND_DUAL_FEASIBLE\n', opt.info)
    
            % Check obtained decision variables
            Pinv = inv(double(P));
            Ln = replace(Ln,vertcat(Ln_coeff{:}),value(vertcat(Ln_coeff{:})));
            Kn = Ln / double(P);
        case 1
            fprintf(2,'Optimization status: %s: DUAL_INFEASIBLE\n', opt.info)
            error('Controller design not successful!')
        case 2
            fprintf(2,'Optimization status: %s: PRIMAL_INFEASIBLE\n', opt.info)
            error('Controller design not successful!')
        case 4
            fprintf(2,'Optimization status: %s: UNKNOWN\n', opt.info)
            error('Controller design not successful!')
        otherwise
            fprintf(2,'Optimization status: %s\n', opt.info)
            error('Controller design not successful!')
    end
end
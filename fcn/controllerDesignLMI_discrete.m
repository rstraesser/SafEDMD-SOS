function [K,Kw,Pinv,sys,compTime] = controllerDesignLMI_discrete(sys,eps,param)
%% Function to design a controller for a discrete-time nonlinear system
% Inputs: 
%   - sys: system description
%   - eps: small margin to ensure positive definiteness
%   - param: parameters defining the data generation and area of interest
%
% Outputs: 
%   - K: (m x N) controller gain (numerator)
%   - Kw: (m x mN) controller gain (denominator)
%   - Pinv: (N x N)-dimensional Lyapunov matrix 
%   - sys: system description
%   - compTime: computation time of the controller design
%
% __author__ = "Robin Straesser"
% __contact__ = "robin.straesser@ist.uni-stuttgart.de"
% __date__ = "2025/01/17"

% Optimization: Decision variables
lmis = [];
    % (Inverse) Lyapunov matrix
    P = sdpvar(sys.N,sys.N,'symmetric');
    lmis = [lmis,P >= eps.P*eye(sys.N)];
    % State-feedback gain
    L = sdpvar(sys.m,sys.N,'full');
    % Full-information feedback gain
    Lw = sdpvar(sys.m,sys.m*sys.N,'full');

    % multiplier bilinearity
    Piinv = inv([sys.Pi.Qz,sys.Pi.Sz;sys.Pi.Sz',sys.Pi.Rz]);
    tQz = Piinv(1:sys.N,1:sys.N);
    tSz = Piinv(1:sys.N,sys.N+1:end);
    tRz = Piinv(sys.N+1:end,sys.N+1:end);
    Lambda = sdpvar(sys.m,sys.m,'symmetric');
    lmis = [lmis,Lambda - eps.Lambda*eye(sys.m) >= 0];
    
    if param.cx < 0 || param.cu < 0
        error('The matrix defining the proportional bound needs to be positive!')
    elseif param.cx == 0 && param.cu == 0
        tau = 0;
        F13 = []; F23 = []; F33 = []; F34 = []; F35 = [];
    else
        % multiplier residual
        tau = sdpvar(1);
        lmis = [lmis,tau >= eps.tau];
        
        F13 = zeros(sys.N,sys.N+sys.m);
        F23 = -kron(eye(sys.m),tSz')*[zeros(sys.N*sys.m,sys.N),param.cu*Lw'];
        F33 = tau*blkdiag(1/2*eye(sys.N),1/2*eye(sys.m));
        F34 = [param.cx*P;param.cu*L];
        F35 = -[zeros(sys.N,sys.N*sys.m);param.cu*Lw];
    end
    
    F11 = P - tau*eye(sys.N);
    F12 = -sys.tB*kron(Lambda,tSz) - sys.B0*Lw*kron(eye(sys.m),tSz);
    F14 = sys.A*P + sys.B0*L;
    F15 = sys.tB*kron(Lambda,eye(sys.N)) + sys.B0*Lw;
    F22 = kron(Lambda,tRz) - Lw*kron(eye(sys.m),tSz) - (Lw*kron(eye(sys.m),tSz))';
    F24 = L;
    F25 = Lw;
    F44 = P;
    F45 = zeros(sys.N,sys.N*sys.m);
    F55 = -kron(Lambda,inv(tQz));

    F = [F11 ,F12 ,F13 ,F14 ,F15 ;
         F12',F22 ,F23 ,F24 ,F25 ;
         F13',F23',F33 ,F34 ,F35 ;
         F14',F24',F34',F44 ,F45 ;
         F15',F25',F35',F45',F55];
    
    lmis = [lmis,F - eps.F*eye(size(F)) >= 0];
    
    % invariance of region of attraction
    nu = sdpvar(1);
    lmis = [lmis,nu >= eps.nu];
    FI11 = nu*tRz - 1;
    FI12 = -nu*tSz';
    FI22 = nu*tQz + P;
    FI = [ FI11 ,FI12 ;
           FI12',FI22];
    lmis = [lmis,FI <= 0];
    
    % Solve optimization
    cost = -trace(P);
    [~,opt] = evalc("optimize(lmis,cost,sdpsettings('solver','mosek'))");

    compTime = opt.solvertime;
    K = NaN(size(L));
    Kw = NaN(size(Lw));
    Pinv = NaN(size(P));
    switch opt.problem 
        case -1
            fprintf(2,'Controller design was not successful: %s: ILL_POSED', opt.info)
        case 0 
            fprintf('Controller design completed. %s: PRIMAL_AND_DUAL_FEASIBLE\n', opt.info)
            % Store obtained decision variables
            P = double(P);
            Pinv = P \ eye(sys.N);
            Lambda = value(Lambda);
        
            K = double(L) / P;
            Kw = double(Lw)*kron(inv(Lambda),eye(sys.N));
        case 1
            fprintf(2,'Controller design was not successful: %s: DUAL_INFEASIBLE\n', opt.info)
        case 2
            fprintf(2,'Controller design was not successful: %s: PRIMAL_INFEASIBLE\n', opt.info)
        case 4
            fprintf(2,'Controller design was not successful: %s: UNKNOWN\n', opt.info)
        otherwise
            fprintf(2,'Controller design was not successful: %s\n', opt.info)
    end
end
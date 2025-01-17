function [param,sys,X,Y] = SafEDMD(X0,X1,sys,param)
%% Function for data generation 
% Inputs: 
%   - X0: data matrix (states)
%   - X1: data matrix (xp) for xp = f(x) + g(x)u
%   - sys: system description
%   - param: parameters defining the data generation and area of interest
%
% Outputs:  
%   - param: parameters defining the data generation and area of interest
%   - sys: system description
%   - X: lifted data matrix (lifted states)
%   - Y: lifted data matrix (lifted xp) for xp = f(x) + g(x)u
%
% __author__ = "Robin Straesser"
% __contact__ = "robin.straesser@ist.uni-stuttgart.de"
% __date__ = "2025/01/17"
%% Calculate data-driven approximation of the Koopman operator
fprintf('Apply SafEDMD to calculate the bilinear surrogate model...')
sys.N = size(param.Phi(zeros(sys.n,1)),1)-1;
param.hPhi = @(x) [zeros(sys.N,1),eye(sys.N)]*param.Phi(x);
    X=cell(sys.m+1,1);Y=cell(sys.m+1,1);

    %% Define the type of data tuples
    switch sys.timeVariant
        case 'continuous-time'
            if isfield(param,'gradhPhi')
                gradhPhi = param.gradhPhi;
            else
                % assumes Phi to be polynonimal
                xSym = sdpvar(sys.n,1);
                gradhPhiSym = full(jacobian(param.hPhi(xSym),xSym))';
                gradhPhi = @(x) value(replace(gradhPhiSym,xSym,x));
            end
            y = @(X0,X1) gradhPhi(X0)'*X1;

        case 'discrete-time'
            y = @(X0,X1) param.hPhi(X1);
    end
            
    %% Lift data to the Koopman space
    for k=0:sys.m % each input u=0, u=e_1, ..., u=e_m
        if k == 0
            X{k+1,1}=cell2mat(arrayfun(@(j) param.hPhi(X0{k+1,1}(:,j)), 1:param.d, 'UniformOutput', false));
        else
            X{k+1,1}=cell2mat(arrayfun(@(j) param.Phi(X0{k+1,1}(:,j)), 1:param.d, 'UniformOutput', false));
        end
        Y{k+1,1}=cell2mat(arrayfun(@(j) y(X0{k+1,1}(:,j),X1{k+1,1}(:,j)), 1:param.d, 'UniformOutput', false));
    end
        
    %% Koopman operator approximations
    % input u=0
    sys.A = Y{0+1,1}*pinv(X{0+1});
    sys.B0 = NaN(sys.N,sys.m);
    sys.tB = NaN(sys.N,sys.m*sys.N);
    % input u=e_1, ..., u=e_m
    for i=1:sys.m 
        [B0i_tBi] = Y{i+1,1}*pinv(X{i+1,1});   
        sys.B0(:,i) = B0i_tBi(:,1);
        sys.tB(:,(i-1)*sys.N+1:i*sys.N) = B0i_tBi(:,2:end) - sys.A;
    end
    sys.A(abs(sys.A)<1e-12) = 0;
    sys.B0(abs(sys.B0)<1e-12) = 0;
    sys.tB(abs(sys.tB)<1e-12) = 0;
    fprintf('Done.\n')
    
    %% Check stabilizability of (A,B0)
    if rank(ctrb(sys.A,sys.B0)) ~= sys.N
        eigA = eig(sys.A);
        eigA = eigA(abs(eigA)>=1);
        stabilizable = 1;
        for i = 1:length(eigA)
            if rank([eigA(i)*eye(sys.N) - sys.A,sys.B0]) < sys.N
                stabilizable = 0;
            end
        end
        if stabilizable
            fprintf('Linear system part (A,B0) is stabilizable.\n')
        else
            fprintf(2,'Linear system part (A,B0) is NOT stabilizable!\n')
        end
    else
        fprintf('Linear system part (A,B0) is controllable.\n')
    end
end
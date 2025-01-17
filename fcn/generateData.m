function [X0,X1,U] = generateData(sys,param)
%% Function for data generation 
% Inputs: 
%   - sys: system description
%   - param: parameters defining the data generation and area of interest
%
% Outputs:  
%   - X0: data matrix (states)
%   - X1: data matrix (xp) for xp = f(x) + g(x)u
%   - U: data matrix (inputs)
%
% __author__ = "Robin Straesser"
% __contact__ = "robin.straesser@ist.uni-stuttgart.de"
% __date__ = "2025/01/17"
%% Collect data of nonlinear system
randuni = @(dim) (2*rand(dim(1),dim(2))-1);
xwidth = (param.xmax - param.xmin)/2;
xcenter = param.xmax - xwidth;
time = tic;fprintf('Simulate measurement data and arrange it in matrices...')
    X0=cell(sys.m+1,1);U=cell(sys.m+1,1);X1=cell(sys.m+1,1);
    I_m = eye(sys.m);
    for k=0:sys.m % each input u=0, u=e_1, ..., u=e_m
        X0{k+1,1} = xcenter + xwidth.*randuni([sys.n,param.d]);
        if k == 0
            U{k+1,1} = zeros(sys.m,param.d);
        else
            U{k+1,1} = repmat(I_m(:,k),[1,param.d]);
        end
        X1{k+1,1}=NaN(sys.n,param.d);
        for j=1:param.d % uniformly sampled data
            switch sys.timeVariant
                case 'continuous-time'
                    X1{k+1,1}(:,j) = sys.ode(X0{k+1,1}(:,j),U{k+1,1}(:,j));
                case 'discrete-time'
                    [~,xnext] =  ode45(@(t,x) sys.ode(x,U{k+1,1}(:,j)),[0,param.DeltaT],X0{k+1,1}(:,j)); 
                    X1{k+1,1}(:,j) = xnext(end,:)';
                otherwise
                    error("Specify 'sys.timeVariant' as either 'discrete-time' or 'continuous-time'!")
            end
        end
    end
time = toc(time);fprintf('Done. Time: %fs\n',time)
end
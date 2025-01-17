function [Q_poly,Q_coeff] = polymatrix(x,d,dim,flag)    
%% Function for data generation 
% Inputs: 
%   - x: state of the polynomial matrix
%   - d: degree of the polynomial matrix
%   - dim: dimension of the polynomial matrix
%   - flag: symmetric property of the polynomial matrix
%
% Outputs:  
%   - Q_poly: polynomial matrix
%   - Q_coeff: coefficients of polynomial matrix
%
% __author__ = "Robin Straesser"
% __contact__ = "robin.straesser@ist.uni-stuttgart.de"
% __date__ = "2025/01/17"
%  Constructs Yalmip polynomial matrix  
    if dim(1) ~= dim(2) && isequal(flag,'symmetric')
        error('Symmetric matrix must be square!')
    end
    Q_poly = sdpvar(dim(1),dim(2));
    Q_coeff = cell(dim(1),dim(2));
    for j = 1:dim(2)
        switch flag 
            case 'symmetric'
                for i = 1:j
                    [Q_poly(i,j),Q_coeff{i,j}] = polynomial(x,d);
                    if i ~= j
                        Q_poly(j,i) = Q_poly(i,j);
                        Q_coeff{j,i} = Q_coeff{i,j};
                    end
                end
            case 'full'
                for i = 1:dim(1)
                    [Q_poly(i,j), Q_coeff{i,j}] = polynomial(x,d);
                end
        end
    end
end
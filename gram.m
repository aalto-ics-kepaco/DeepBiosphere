function [Kx] = gram( X1, X2, kernel_type, kernel_param)

% This function computes the kernel (gram) matrix between the matrices X1
% and X2.

% Input Arguments:
%   X1              Matrix with n1 rows and p columns
%   X2              Matrix with n2 rows and p columns
%   kernel_type     'linear' or 'gaussian'
%   kernel_param    If 'linear' was chosen, put to [] or leave out. If 
%                   'gaussian' was chosen, this parameter corresponds to 
%                   the width of the kernel.

% Output Arguments:
%   Kx              Kernel matrix with n1 rows and n2 columns

% © 08/06/2015 Viivi Uurtio, Aalto University
% viivi.uurtio@aalto.fi
%
% This code is for academic purposes only.
% Commercial use is not allowed.

    n1 = size(X1,2);
    n2 = size(X2,2);
    
    switch kernel_type
        
        case 'linear'
            Kx = X1'*X2;  
            % normalization
            Kx1 = X1'*X1;
            Kx2 = X2'*X2;
            Kx = Kx./ sqrt(repmat(diag(Kx1),1,n2) .* repmat(diag(Kx2)',n1,1));
            
        case 'gaussian'
            D = diag(X1'*X1)*ones(1,n2) + ones(n1,1)*diag(X2'*X2)' - 2*X1'*X2;
            Kx = exp(- kernel_param * D);
            
        otherwise
            warning('Unexpected kernel type.')
    end
    
end


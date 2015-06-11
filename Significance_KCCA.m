function [observed_corr, p_values] = Significance_KCCA(XN,YN,kernel_type,projections,perm_nr,kappa,kernel_param)

% Statistical significance tests on parameter values by permutation tests

% Input Arguments:
%   XN              Normalized and centered data in view X, , the rows
%                   correspond to the samples (examples) and the columns to
%                   the variables (or features).
%   YN              Normalized and centered data in view Y, , the rows
%                   correspond to the samples (examples) and the columns to
%                   the variables (or features).
%   projections     The number of projections.
%   perm_nr         Number of permutations
%   kernel_type     'linear' or 'gaussian'
%   kappa           The regularization parameter, 
%                   (if left empty, i.e. [], default kappa=1)
%   kernel_param    If 'linear' was chosen, put to [] or leave out. If 
%                   'gaussian' was chosen, this parameter corresponds to 
%                   the width of the kernel.

% Output Arguments:
%   observed_corr   Canonical correlation obtained from the data
%                   with the cross-validated parameters.
%   p_values        Statistical significance of the canonical correlation
%                   obtained from the data

% © 09/06/2015 Viivi Uurtio, Aalto University
% viivi.uurtio@aalto.fi
%
% This code is for academic purposes only.
% Commercial use is not allowed.


if ~exist('kappa','var') || isempty(kappa),
    kappa = 1;
end

if ~exist('kernel_param','var') || isempty(kernel_param),
    kernel_param = 1;
end

rng(1,'v5normal');
l = size(XN,1);

switch kernel_type
    case 'linear'
        KXN = XN*XN'; %Kernel for 1st view
        KYN = YN*YN'; %Kernel for 2nd view
        [~, ~, correlations] = kcanonca_reg_ver2(KXN,KYN,0.1,kappa);
        
        for i = 1:perm_nr %random permutations
            perm = randperm(l);
            YN_perm=YN(perm,:); %random permutation of the rows
            KYN = YN_perm*YN_perm';% Kernel for the permuted view
            [~, ~, permuted_can_corr] = kcanonca_reg_ver2(KXN,KYN,0.1,kappa);
            permuted_corr(i,:) = permuted_can_corr;
        end
        
        
    case 'gaussian'
        width=kernel_param;
        KXN = XN*XN'; %Kernel for 1st view
        KYN = gram( YN', YN', 'gaussian', width); % Kernel for the 2nd view
        [~, ~, correlations] = kcanonca_reg_ver2(KXN,KYN,0.1,kappa);
        
        for i = 1:perm_nr %random permutations
            perm = randperm(l);
            YN_perm=YN(perm,:); %random permutation of the rows
            KYN = gram( YN_perm', YN_perm', 'gaussian', width);
            [~, ~, permuted_can_corr] = kcanonca_reg_ver2(KXN,KYN,0.1,kappa);
            permuted_corr(i,:) = permuted_can_corr;
        end
end

% Computing p-values
observed_corr=correlations(1:projections,1)';
perm_corr=permuted_corr(:,1:projections);
p_values=zeros(1,projections);
for k=1:projections
    greater=sum(perm_corr(:,k) >= observed_corr(k));
    p_values(k)=greater/perm_nr;
end


end
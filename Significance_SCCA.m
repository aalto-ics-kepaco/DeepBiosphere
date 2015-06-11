function [observed_corr, p_values] = Significance_SCCA(XN,YN,mu,kernel_type,projections,perm_nr,kernel_param)
% Statistical significance tests on parameter values by permutation tests

% Input Arguments:
%   XN              Normalized and centered data in view X, , the rows
%                   correspond to the samples (examples) and the columns to
%                   the variables (or features).
%   YN              Normalized and centered data in view Y, , the rows
%                   correspond to the samples (examples) and the columns to
%                   the variables (or features).
%   mu              A vector for the range of parameter values that are
%                   tested in the cross-validation, e.g. 0.1:0.1:1
%   kernel_type     'linear' or 'gaussian'
%   projections     The number of projections.
%   perm_nr         Number of permutations
%   kernel_param    If 'linear' was chosen, put to [] or leave out. If 
%                   'gaussian' was chosen, this parameter corresponds to 
%                   a vector of values of the width of the kernel that
%                   are tested in cross validation, e.g. 10.^(-2:0.1:2).

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

%cvx_solver sedumi
rng(1,'v5normal');
permuted_corr=zeros(perm_nr,projections);
l = size(XN,1);

switch kernel_type
    case 'linear'
        [mu_star, ~] = Cross_Validation(XN,YN,mu);
        KYN = YN*YN';
        [seeds,~] = spectral_clustering(YN,projections);
        [~, ~, ~, observed_corr, ~, ~]...
            = SCCAwrapper2_corrected_corr_cvx(XN',KYN,seeds,1,0,mu_star,1);
        disp(['OBS CORR ' num2str(observed_corr)])
        
        for i = 1:perm_nr %random permutations
            disp(['PERMUTATION ' num2str(i)])
            perm = randperm(l);
            YN_perm=YN(perm,:); %random permutation of the rows
            KYN = YN_perm*YN_perm';% Kernel for the permuted view
            
            [seeds,~] = spectral_clustering(YN_perm,projections);
            
            [~, ~, ~, permuted_can_corr, ~, ~]...
                = SCCAwrapper2_corrected_corr_cvx(XN',KYN,seeds,1,0,mu_star,2);
            permuted_corr(i,:) = permuted_can_corr;
            disp(['PERM CORR ' num2str(permuted_can_corr)])
        end
            
    case 'gaussian'
        [mu_star, width_star,~] = CrossValidation_SparsityWidth(XN,YN,mu,kernel_param);
        KYN = gram( YN', YN', 'gaussian', width_star);
        [seeds,~] = spectral_clustering(YN,projections);
        [~, ~, ~, observed_corr, ~, ~]...
            = SCCAwrapper2_corrected_corr_cvx(XN',KYN,seeds,1,0,mu_star,1);
        
        for i = 1:perm_nr %random permutations
            disp(['PERMUTATION ' num2str(i)])
            perm = randperm(l);
            YN_perm=YN(perm,:); %random permutation of the rows
            KYN = gram( YN_perm', YN_perm', 'gaussian', width_star);
            [seeds,~] = spectral_clustering(YN_perm,projections);
            [~, ~, ~, permuted_can_corr, ~, ~]...
                = SCCAwrapper2_corrected_corr_cvx(XN',KYN,seeds,1,0,mu_star,2);
            permuted_corr(i,:) = permuted_can_corr;
            disp(['PERM CORR ' num2str(permuted_can_corr)])
        end
end

% Computing p-values
p_values=zeros(1,projections);
for k=1:projections
    greater=sum(permuted_corr(:,k) >= observed_corr(k));
    p_values(k)=greater/perm_nr;
end


end


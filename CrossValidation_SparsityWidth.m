function [mu_star, width_star, value] = CrossValidation_SparsityWidth(XN,YN,mu,width)

% This function computes the value of the sparsity parameter mu that
% maximizes the predictive canonical correlation coefficient when a
% Gaussian kernel is used for the view Y. 

% Input Arguments:
%   XN              Normalized and centered data in view X, the rows
%                   correspond to the samples (examples) and the columns to
%                   the variables (or features).
%   YN              Normalized and centered data in view Y, the rows
%                   correspond to the samples (examples) and the columns to
%                   the variables (or features).
%   mu              A vector for the range of parameter values that are
%                   tested in the cross-validation, e.g. 0.1:0.1:1
%   width           A vector for the range of kernel width values that are
%                   tested in the cross-validation, e.g. 10.^(-2:0.2:2)

% Output Arguments:
%   mu_star         The level of sparsity that yields a maximal predictive
%                   canonical correlation coefficient
%   width_star      The width of the kernel that yields a maximal predictive
%                   canonical correlation coefficient
%   value           The value of the maximal predictive canonical
%                   correlation coefficient

% © 09/06/2015 Viivi Uurtio, Aalto University
% viivi.uurtio@aalto.fi
%
% This code is for academic purposes only.
% Commercial use is not allowed.

disp('CROSS VALIDATION')

rng(1,'v5normal')
N=size(XN,1); %Total number of observations, i.e. rows
K=3;
Indices = crossvalind('Kfold', N, K);
projections=1;

% Preallocation for speed
q=zeros(size(mu,2),size(width,2));
q_cv=zeros(size(mu,2),size(width,2),K);
corrs=zeros(size(mu,2),size(width,2));

for k=1:size(mu,2)
    for j=1:size(width,2)
        KYN = gram(YN', YN', 'gaussian', width(j)); % Gaussian kernel
        [seeds,~] = spectral_clustering(YN,projections);
        
        [WsccaD1, ZsccaD1, ~, observed_corr1, ~, ~]...
            = SCCAwrapper2_corrected_corr_cvx(XN',KYN,...
            seeds,1,0,mu(k),2);
        corrs(k,j)=observed_corr1;
        
        score_a = XN*WsccaD1;
        score_b = KYN'*ZsccaD1;
        q(k,j) = score_a'*score_b/(norm(score_a)*norm(score_b));
        
        for i=1:K
            disp(['FOLD:  ' num2str(i)])
            
            testInd = (Indices == i);
            trainInd = ~testInd;
            
            trainXN=XN(trainInd,:);
            trainYN=YN(trainInd,:);
            
            testXN=XN(testInd,:);
            testYN=YN(testInd,:);
            
            [seeds1,~] = spectral_clustering(trainYN,projections);
            
            % TRAIN SET: computing the weights
            trainKYN = gram( trainYN', trainYN', 'gaussian', width(j)); % Gaussian kernel
            [WsccaD1_cv, ZsccaD1_cv, ~, ~, ~, ~]...
                = SCCAwrapper2_corrected_corr_cvx(trainXN',trainKYN,...
                seeds1,1,0,mu(k),2);%deflation is not required at this stage
            
            % TEST SET:
            testKYN = gram( trainYN', testYN', 'gaussian', width(j)); % Gaussian kernel
            % We compute the scores for the hold out set
            score_a = testXN*WsccaD1_cv;
            score_b = testKYN'*ZsccaD1_cv;
            if norm(score_a)*norm(score_b)==0
                q_cv(k,j,i)=0;
            else
                q_cv(k,j,i) = max(score_a'*score_b/(norm(score_a)*norm(score_b)));
            end
            
        end
    end
end

CV1=mean(q_cv,3);

[value, max_ind]=max(CV1(:));
[u,v] = ind2sub(size(CV1),max_ind);
mu_star=mu(u);
width_star=width(v);

disp('--------------------------')
disp(['CV mu:  ' num2str(mu_star)])
disp('--------------------------')
disp('--------------------------')
disp(['CV WIDTH:  ' num2str(width_star)])
disp('--------------------------')

end
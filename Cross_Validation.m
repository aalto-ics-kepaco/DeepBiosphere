function [mu_star, value] = Cross_Validation(XN,YN,mu)

% This function computes the value of the sparsity parameter mu that
% maximizes the predictive canonical correlation coefficient when a linear
% kernel is used for the view Y.

% Input Arguments:
%   XN              Normalized and centered data in view X, , the rows
%                   correspond to the samples (examples) and the columns to
%                   the variables (or features).
%   YN              Normalized and centered data in view Y, , the rows
%                   correspond to the samples (examples) and the columns to
%                   the variables (or features).
%   mu              A vector for the range of parameter values that are
%                   tested in the cross-validation, e.g. 0.1:0.1:1

% Output Arguments:
%   mu_star         The level of sparsity that yields a maximal predictive
%                   canonical correlation coefficient
%   value           The value of the maximal predictive canonical
%                   correlation coefficient

% © 09/06/2015 Viivi Uurtio, Aalto University
% viivi.uurtio@aalto.fi
%
% This code is for academic purposes only.
% Commercial use is not allowed.

disp('CROSS VALIDATION')
% Split the data randomly into k folds, each example in one fold.
rng(1,'v5normal')
N=size(XN,1); %Total number of observations, i.e. rows
K=3;
Indices = crossvalind('Kfold', N, K);

projections=1; 
density=zeros(K,size(mu,2));
q=zeros(K,size(mu,2));


for i=1:K 
    
    testInd = (Indices == i); 
    trainInd = ~testInd;
    
    trainXN=XN(trainInd,:);
    trainYN=YN(trainInd,:);
    trainKYN = trainYN*trainYN';
    
    [seeds,~] = spectral_clustering(trainYN,projections);
        
    testXN=XN(testInd,:);
    testYN=YN(testInd,:);
    testKYN = trainYN*testYN'; 
    
    disp(['FOLD: ' num2str(i)])
    disp(['SEED: ' num2str(seeds)])
       
    for j=1:length(mu)
        
        % TRAIN SET: computing the weights
        [WsccaD, ZsccaD, ~, ~, ~, ~]...
            = SCCAwrapper2_corrected_corr_cvx(trainXN',trainKYN,...
            seeds,1,0,mu(j),2);%deflation is not required at this stage
       
        % TEST SET:
        % We compute the scores for the hold out set
        score_a = testXN*WsccaD;
        score_b = testKYN'*ZsccaD;
        if norm(score_a)*norm(score_b)==0
            q(i,j)=0;
        else
            q(i,j) = max(score_a'*score_b/(norm(score_a)*norm(score_b)));
        end
        density(i,j)=nnz(WsccaD)/numel(WsccaD);

    end
end

% The optimal mu_star is the one that maximizes q
[value, I]=max(median(q));
mu_star=mu(I);

disp('--------------------------')
disp(['CV mu:  ' num2str(mu_star)])
disp('--------------------------')



end







function [mu_star, width_star, value]=PredictiveCC_GaussianKernel(XN,YN,mu,width)
% This function computes the value of the sparsity parameter mu that
% maximizes the predictive canonical correlation coefficient when a
% Gaussian kernel is used for the view Y. The predictive canonical correlation
% coefficient, canonical correlation coefficient and density of the data
% matrix are visualized.

% Input Arguments:
%   XN              Normalized and centered data in view X, , the rows
%                   correspond to the samples (examples) and the columns to
%                   the variables (or features).
%   YN              Normalized and centered data in view Y, , the rows
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

rng(1,'v5normal')
N=size(XN,1); %Total number of observations, i.e. rows
K=3; % 3-fold cross validation
Indices = crossvalind('Kfold', N, K); % divide the samples into 3 folds
projections=1;

% Preallocation for speed
density=zeros(1,length(width));
q=zeros(size(mu,2),size(width,2));
q_cv=zeros(size(mu,2),size(width,2),K);
corrs=zeros(size(mu,2),size(width,2));

for k=1:size(mu,2)
    for j=1:size(width,2)
        
        [seeds,~] = spectral_clustering(YN,projections);
        
        KYN = gram(YN', YN', 'gaussian', width(j)); % Gaussian kernel
        [WsccaD1, ZsccaD1, ~, observed_corr1, ~, ~]...
            = SCCAwrapper2_corrected_corr_cvx(XN',KYN,...
            seeds,1,0,mu(k),2);
        corrs(k,j)=observed_corr1;
        
        score_a = XN*WsccaD1;
        score_b = KYN'*ZsccaD1;
        q(k,j) = score_a'*score_b/(norm(score_a)*norm(score_b));
        density(j)=nnz(WsccaD1)/numel(WsccaD1);
        
        for i=1:K
            
            testInd = (Indices == i);
            trainInd = ~testInd;
            
            trainXN=XN(trainInd,:);
            trainYN=YN(trainInd,:);
            
            testXN=XN(testInd,:);
            testYN=YN(testInd,:);
            
            [seeds,~] = spectral_clustering(trainYN,projections);
            disp(['SEEDS ' num2str(seeds)])
            shuffle_seeds(i,:)=seeds;
            
            % TRAIN SET: computing the weights
            trainKYN = gram( trainYN', trainYN', 'gaussian', width(j)); % Gaussian kernel
            [WsccaD1_cv, ZsccaD1_cv, ~, ~, ~, ~]...
                = SCCAwrapper2_corrected_corr_cvx(trainXN',trainKYN,...
                seeds,1,0,mu(k),2); % deflation is not required at this stage
            
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
            fprintf('width: %.1f \n',width(j))
            fprintf('mu: %.1f \n',mu(k))
            
        end
    end
end

CV1=mean(q_cv,3);% mean over the cross-validated correlations

% figure on predictive canonical correlation
figure 
imagesc(CV1);
set(gca,'YTick', 1:size(mu,2),'YTickLabel',sprintf('%3.2f\n',mu));
set(gca,'XTick', 1:size(width,2),'XTickLabel',sprintf('%3.1f\n',width));
xlabel('Width')
ylabel('Scaling Factor')
h = colorbar;
ylabel(h, 'Predictive Correlation');
set(h,'fontsize',10);
set(gca, 'FontName', 'Arial','fontsize',10);
set(gcf, 'PaperPositionMode', 'auto');

[value, max_ind]=max(CV1(:));
[m,n] = ind2sub(size(CV1),max_ind);
mu_star=mu(m);
width_star=width(n);

% figure on canonical correlation obtained from the data
figure 
imagesc(q);
set(gca,'YTick', 1:size(mu,2),'YTickLabel',sprintf('%3.2f\n',mu));
set(gca,'XTick', 1:size(width,2),'XTickLabel',sprintf('%3.1f\n',width));
xlabel('Width')
ylabel('Scaling Factor')
h = colorbar;
ylabel(h, 'Canonical Correlation');
set(h,'fontsize',10);
set(gca, 'FontName', 'Arial','fontsize',10);


end


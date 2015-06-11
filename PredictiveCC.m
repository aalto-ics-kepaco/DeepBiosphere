function [mu_star, value]=PredictiveCC(XN,YN,mu)
% This function computes the value of the sparsity parameter mu that
% maximizes the predictive canonical correlation coefficient when a linear
% kernel is used for the view Y. The predictive canonical correlation
% coefficient, canonical correlation coefficient and density of the data
% matrix are visualized.

% Input Arguments:
%   XN              Normalized and centered data in view X, the rows
%                   correspond to the samples (examples) and the columns to
%                   the variables (or features).
%   YN              Normalized and centered data in view Y, the rows
%                   correspond to the samples (examples) and the columns to
%                   the variables (or features).
%   mu              A vector for the range of parameter values that are
%                   tested in the cross-validation, e.g. 0.1:0.1:1

% Output Arguments:
%   mu_star         The level of sparsity that yields a maximal predictive
%                   canonical correlation coefficient
%   value           The value of the maximal predictive canonical
%                   correlation coefficient

rng(1,'v5normal');
N=size(XN,1); %Total number of observations, i.e. rows
K=3;
Indices = crossvalind('Kfold', N, K);

projections=1;
density=zeros(1,size(mu,2));
q=zeros(1,size(mu,2));
KYN=YN*YN';

for i=1:K
    
    testInd = (Indices == i);
    trainInd = ~testInd;
    
    trainXN=XN(trainInd,:);
    trainYN=YN(trainInd,:);
    trainKYN = trainYN*trainYN';
    
    [seeds,~] = spectral_clustering(trainYN,projections);
    disp(['SEED: ' num2str(seeds)])
    
    testXN=XN(testInd,:);
    testYN=YN(testInd,:);
    testKYN = trainYN*testYN';
    
    disp(['FOLD: ' num2str(i)])
    for j=1:size(mu,2)
        [WsccaD1, ZsccaD1, ~, ~, ~, ~]...
            = SCCAwrapper2_corrected_corr_cvx(XN',KYN,...
            seeds,1,0,mu(j),2);
        
        score_a = XN*WsccaD1;
        score_b = KYN'*ZsccaD1;
        q(j) = score_a'*score_b/(norm(score_a)*norm(score_b));
        density(j)=nnz(WsccaD1)/numel(WsccaD1);
        
        % TRAIN SET: computing the weights
        [WsccaD1_cv, ZsccaD1_cv, ~, ~, ~, ~]...
            = SCCAwrapper2_corrected_corr_cvx(trainXN',trainKYN,...
            seeds,1,0,mu(j),2);% deflation is not required at this stage, 
        % the last parameter is 2.
        
        % TEST SET:
        % We compute the scores for the hold out set
        score_a = testXN*WsccaD1_cv;
        score_b = testKYN'*ZsccaD1_cv;
        if norm(score_a)*norm(score_b)==0
            q_cv(i,j)=0;
        else
            q_cv(i,j) = max(score_a'*score_b/(norm(score_a)*norm(score_b)));
        end
        fprintf('mu: %.3f \n',mu(j))
        %density_cv(i,j)=nnz(WsccaD1_cv)/numel(WsccaD1_cv);  
        
    end
end
CV=mean(q_cv);
[value,n]=max(CV);
mu_star=mu(n);

figure;
plot(mu,q);
hold on
plot(mu,CV)
plot(mu,density)
line([mu(n) mu(n)], [-1 1],'color','k');
line([0 max(mu)], [0 0],'color','k','LineStyle',':');
h=legend('Canonical Correlation','Predictive Correlation','Density of the Matrix','location','northeast');
p = get(h, 'pos');
p(1)=p(1)-0.08;
p(2)=p(2)+0.008;
set(h, 'pos', p);
set(h,'FontSize',7,'Interpreter','Latex')
legend('boxoff')
xlabel('$\mu$','fontsize',12,'Interpreter','Latex')
ylabel('Correlation and Density','Interpreter','Latex')
axis([min(mu) max(mu) -1 1])
set(gca,'XTick', mu(n), 'XTickLabel', num2str(mu(n)))
set(gca,'LooseInset',get(gca,'TightInset'))
axis square



end


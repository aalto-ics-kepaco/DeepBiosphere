function ScorePlot(analysis_type,kernel_type,XN,YN,projections,sample_name,analysis_param,kernel_param)

% This function runs either sparse canonical correlation analysis
% (scca) or kernel canonical correlation analysis (kcca) on a given
% two-view dataset and visualizes the projected data points on a score plot.

% Input Arguments:
%   analysis_type   'scca' for sparse canonical correlation analysis 
%                   'kcca' for kernel canonical correlation analysis
%   analysis_param  For 'scca': The level of sparsity, mu, in scca 
%                   (if left empty, i.e. [], default mu=1)
%                   For 'kcca': The regularization parameter, kappa,
%                   (if left empty, i.e. [], default kappa=1)
%   XN              Normalized and centered data in view X, the rows
%                   correspond to the samples (examples) and the columns to
%                   the variables (or features).
%   YN              Normalized and centered data in view Y, the rows
%                   correspond to the samples (examples) and the columns to
%                   the variables (or features).
%   sample_name     Names of the samples (Cell array of strings)
%                   If left empty, [], the samples are annotated with
%                   numbers
%   kernel_type     'linear' or 'gaussian'
%   kernel_param    If 'linear' was chosen, put to [] or leave out. If 
%                   'gaussian' was chosen, this parameter corresponds to 
%                   the width of the kernel. If left out, default is 1.


% © 08/06/2015 Viivi Uurtio, Aalto University
% viivi.uurtio@aalto.fi
%
% This code is for academic purposes only.
% Commercial use is not allowed.

                  
% Default parameters
if ~exist('analysis_param','var') || isempty(analysis_param)
    analysis_param = 1;
end

if ~exist('sample_name','var') || isempty(sample_name)
   a = (1:size(XN,1))'; b = num2str(a); c = cellstr(b);
   sample_name=c;
end

if ~exist('kernel_param','var') || isempty(kernel_param)
    kernel_param = 1;
end

switch kernel_type
    case 'linear'
        KXN = XN*XN'; %Kernel for 1st view
        KYN = YN*YN'; %Kernel for 2nd view
           
    case 'gaussian'
        KXN = XN*XN'; %Kernel for 1st view
        KYN = gram( YN', YN', 'gaussian', kernel_param); % Kernel for 2nd view               
end


% We select the two leading canonical projections
pro1=1;
pro2=2;

switch analysis_type
    % computing the scores (canonical variates)
    case 'scca'
        if projections < 2
            warning('Number of projections set to 2')
            projections=2;
        end
        scores=scca(XN,KYN,analysis_param);
    case 'kcca'
        scores=kcca(XN,KYN,analysis_param);
end

figure;
set(gca,'LooseInset',get(gca,'TightInset'))

hold on
h=zeros(size(sample_name,1),1);
for n=1:size(h,1)
    h(n)=plot(scores(n,pro1),scores(n,pro2),'bo');
end

[xmin,~]=min(min(scores(:,pro1:pro2)));[xmax,~]=max(max(scores(:,pro1:pro2)));
r=[xmin xmax];
rmax=max(abs(r));
r_new=[-rmax rmax -rmax rmax];
axis(r_new)
set(gca,'YTick',[-rmax 0 rmax],'YTickLabel',{'-','0','+'});
set(gca,'XTick',[-rmax 0 rmax],'XTickLabel',{'-','0','+'});
set(gca, 'xcolor', 'k', 'ycolor', 'k'); 
plot([r_new(1);r_new(2)],[0;0],'k-') 
plot([0;0],[r_new(3);r_new(4)],'k-')
box on

title('\bf Score Plot','FontSize',14,'Interpreter','latex')
xlabel('\bf 1$^{st}$ Projection','fontsize',14,'Interpreter','latex')
ylabel('\bf 2$^{nd}$ Projection','fontsize',14,'Interpreter','latex')
axis square

set(h(:),'Marker','None');
text(scores(:,1), scores(:,2), sample_name,'fontsize',14,'Interpreter','latex');
hold off


    function scores=scca(XN,KYN,analysis_param)
        rng(1,'v5normal')
        [seeds,~] = spectral_clustering(YN,projections);
        [W, beta, ~, ~, ~, ~]...
            = SCCAwrapper2_corrected_corr_cvx(XN',KYN,seeds,1,0,analysis_param,1);
        scores = XN*W;
        %scores = KYN*beta;
    end

    function scores=kcca(XN,KYN,analysis_param)
        [alpha, beta, ~] = kcanonca_reg_ver2(KXN,KYN,0.1,analysis_param);
        scores = KXN*alpha;
        %scores = KYN*beta;    
    end


end


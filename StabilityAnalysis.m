function StabilityAnalysis(XN,YN,mu,kernel_type,projections,kernel_param,mu_star)

% This function generates a plot in which the stability of the parameter mu
% around its optimum (the value that yields maximal predictive canonical
% correlation coefficient) can be assessed. The curves are interpolated
% through the tested values of mu.

% Input Arguments:
%   XN              Normalized and centered data in view X, , the rows
%                   correspond to the samples (examples) and the columns to
%                   the variables (or features).
%   YN              Normalized and centered data in view Y, , the rows
%                   correspond to the samples (examples) and the columns to
%                   the variables (or features).
%   mu              A vector for the range of parameter values
%                   e.g. 0.1:0.1:1
%   projections     The number of projections.
%   kernel_type     'linear' or 'gaussian'
%   kernel_param    If 'linear' was chosen, put to [] or leave out. If 
%                   'gaussian' was chosen, this parameter corresponds to 
%                   the width of the kernel (default is 1).
%   mu_star         The optimal value of mu obtained from cross validation.
%                   Put to [] if not known.

% © 09/06/2015 Viivi Uurtio, Aalto University
% viivi.uurtio@aalto.fi
%
% This code is for academic purposes only.
% Commercial use is not allowed.


if ~exist('kernel_param','var') || isempty(kernel_param)
    kernel_param = 1;
end

switch kernel_type
    case 'linear'
        KYN = YN*YN'; %Kernel for 2nd view
           
    case 'gaussian'
        KYN = gram( YN', YN', 'gaussian', kernel_param); % Kernel for 2nd view               
end
  
% Feature weights at various mu
rng(1,'v5normal');
[seeds,~] = spectral_clustering(YN,projections);

for j=1:length(mu)
    [WsccaD, ~, ~, ~, ~, ~]= SCCAwrapper2_corrected_corr_cvx(XN',KYN,...
        seeds,1,0,mu(j),1);
    p1(:,j)=WsccaD(:,1);
    p2(:,j)=WsccaD(:,2);
end

W=logical(p1);
Z=logical(p2);

for i=1:size(mu,2)
    rhow(i)=corr(W(:,4),W(:,i));
    rhoz(i)=corr(Z(:,4),Z(:,i));   
end

figure
hold on
fitx=linspace(0.1,2,100);
fity = interp1(mu,rhow,fitx,'pchip');
h(1)=scatter(mu,rhow,'DisplayName', '1^{st} projection');
line(fitx,fity,'HandleVisibility','off');
fitx2=linspace(0.1,2,100);
fity2 = interp1(mu,rhoz,fitx2,'v5cubic');
h(2)=scatter(mu,rhoz,'ms','DisplayName', '2^{nd} projection');
line(fitx2,fity2,'color','m','HandleVisibility','off');
[~,~,~,~] = legend('show');
legend('boxoff')
set(gca,'XTick',mu);
box on
axis([min(mu) max(mu) -1 1])
line([0 2], [0.0 0.0],'linestyle',':','color','k');
ylabel(['Linear Correlation Coefficient' 10 'Between the Weight Vectors'],'interpreter','latex')
xlabel('$\mu$','interpreter','latex')

if ~exist('mu_star','var') || isempty(mu_star)
    line([mu_star mu_star], [-1 1],'linestyle','-','color','k');
end


end


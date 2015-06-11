function CorrelationPlot(analysis_type,kernel_type,XN,YN,projections,Xvar_name,Yvar_name,analysis_param,kernel_param)
    
% This function runs either sparse canonical correlation analysis
% (scca) or kernel canonical correlation analysis (kcca) on a given
% two-view dataset and visualizes the correlations between the features on
% a correlation plot.

% Input Arguments:
%   analysis_type   'scca' for sparse canonical correlation analysis 
%                   'kcca' for kernel canonical correlation analysis
%   kernel_type     'linear' or 'gaussian'
%   XN              Normalized and centered data in view X, the rows
%                   correspond to the samples (examples) and the columns to
%                   the variables (or features).
%   YN              Normalized and centered data in view Y, the rows
%                   correspond to the samples (examples) and the columns to
%                   the variables (or features).
%   projections     The number of projections, minimum is 2.
%   Xvar_name       Variable names in view X, if [], variables are
%                   annotated a,b,c,...,aa,ab,ac,...,zx,zy,zz
%   Yvar_name       Variable names in view Y, if [] variables are annotated
%                   by numbers
%   analysis_param  The level of sparsity, mu, in scca 
%                   (if left empty, i.e. [], default mu=1)
%                   The regularization parameter, kappa, in kcca
%                   (if left empty, i.e. [], default kappa=1)
%   kernel_param    If 'linear' was chosen, put to [] or leave out. If 
%                   'gaussian' was chosen, this parameter corresponds to 
%                   the width of the kernel. If left out, default is 1.


% © 08/06/2015 Viivi Uurtio, Aalto University
% viivi.uurtio@aalto.fi
%
% This code is for academic purposes only.
% Commercial use is not allowed.
                  

% Setting the default parameters and type of analysis

if ~exist('analysis_param','var') || isempty(analysis_param),
    analysis_param = 1;
end

Alphabet=('a':'z').';
[I,J] = meshgrid(1:26,1:26);
Xlabels=[Alphabet(I(:)), Alphabet(J(:))];
Xlabels=strvcat(Alphabet,Xlabels);

if ~exist('Xvar_name','var') || isempty(Xvar_name)
    Xvar_name=cellstr(Xlabels(1:size(XN,2),:));
end

if ~exist('Yvar_name','var') || isempty(Yvar_name)
    a = (1:size(YN,2))'; b = num2str(a); c = cellstr(b);
    Yvar_name=c;
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
        [list_indices,X_name]=scca(XN,KYN,analysis_param);
    case 'kcca'
        [list_indices,X_name]=kcca(XN,KYN,analysis_param);
end

figure;
set(gca,'LooseInset',get(gca,'TightInset'))
hold on
h=zeros(size(list_indices));
for n=1:size(h,1)
    h(n)=plot(x(list_indices{n}),y(list_indices{n}),'bo');
end
set(h(:),'Marker','None');
e=cellstr(X_name); % Marking the variables of view X by their names
text(x, y, e,'fontsize',14,'Interpreter','latex'); 
g=plot(z,w,'rs');
set(g,'Marker','None');
f=cellstr(Yvar_name); % Marking the variables of view Y by their names
text(z, w, f,'fontsize',14,'Interpreter','latex'); 

xmin=-1;xmax=1;ymin=-1;ymax=1;
axis([xmin xmax ymin ymax]) 
set(gca,'YTick',[-1 -0.5 0 0.5 1]); 
set(gca,'XTick',[-1 -0.5 0 0.5 1]);
set(gca, 'xcolor', 'k', 'ycolor', 'k'); 
plot([xmin;xmax],[0;0],'k-') 
plot([0;0],[ymin;ymax],'k-') 
circle([0,0],0.5,1000,'k:');
circle([0,0],1,1000,'k:');
box on

title('\bf Correlation Plot','FontSize',14,'Interpreter','latex')
xlabel('\bf 1$^{st}$ Projection','fontsize',14,'Interpreter','latex')
ylabel('\bf 2$^{nd}$ Projection','fontsize',14,'Interpreter','latex')

hold off


    function [list_indices, X_name]=scca(XN,KYN,analysis_param)
        rng(1,'v5normal')
        [seeds,~] = spectral_clustering(YN,projections);
        [W, beta, ~, ~, ~, ~]...
            = SCCAwrapper2_corrected_corr_cvx(XN',KYN,seeds,1,0,analysis_param,1);
        scores = XN*W;
        %scores = KYN*beta;
        
        % indices for the sparse features
        vec1=logical(W(:,pro1));
        vec2=logical(W(:,pro2));
        
        % Scores vs original (Gonzalez et al. 2012, Mevik and Wehrens 2007)
        yycor = corr(YN, scores);
        yycor(isnan(yycor)) = 0;
        xxcor = corr(XN, scores);
        xxcor(isnan(xxcor)) = 0;
        
        % correlation coordinates multiplied
        x=xxcor(:,pro1).*vec1;
        y=xxcor(:,pro2).*vec2;
        
        z=yycor(:,pro1);
        w=yycor(:,pro2);
        
        indices=zeros(size(x,1),1);
        j=size(XN,2);
        for k=1:size(x,1)
            if (x(j,1)==0 && y(j,1)==0)
                x(j)=[];
                y(j)=[];
            else
                indices(k)=1;
            end
            j=j-1;
        end
        indices=flipud(indices);
        sparse=find(indices==1);% indices for variables in view X
        X_name=Xvar_name(sparse);
        
        % Mark each species differently
        C = unique(X_name);
        list_indices=cell(size(C,2),1);
        for i=1:size(C,2);
            IndexC = strfind(X_name, C{1,i});
            Index = find(not(cellfun('isempty', IndexC)));
            list_indices{i}=Index;
        end   
    end

    function [list_indices,X_name]=kcca(XN,KYN,analysis_param)
        [alpha, beta, ~] = kcanonca_reg_ver2(KXN,KYN,0.1,analysis_param);
        scores = KXN*alpha;
        %scores = KYN*beta;
        
        % Scores vs original (Gonzalez et al. 2012, Mevik and Wehrens 2007)
        yycor = corr(YN, scores);
        yycor(isnan(yycor)) = 0;
        xxcor = corr(XN, scores);
        xxcor(isnan(xxcor)) = 0;
        
        x=xxcor(:,pro1);
        y=xxcor(:,pro2);
        
        z=yycor(:,pro1);
        w=yycor(:,pro2);
        
        X_name=Xvar_name;
        % Mark each species differently
        C = unique(X_name);
        list_indices=cell(size(C,2),1);
        for i=1:size(C,2);
            IndexC = strfind(X_name, C{1,i});
            Index = find(not(cellfun('isempty', IndexC)));
            list_indices{i}=Index;
        end
    end


end


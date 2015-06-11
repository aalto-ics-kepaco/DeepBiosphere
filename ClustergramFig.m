function ClustergramFig(analysis_type,kernel_type,XN,YN,projections,Xvar_name,Yvar_name,analysis_param,kernel_param)

% This function visualizes the correlation coefficients between the 
% variables in view X and Y using the analytic geometry of the position
% vectors on the correlation plot.

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

% © 09/06/2015 Viivi Uurtio, Aalto University
% viivi.uurtio@aalto.fi
%
% This code is for academic purposes only.
% Commercial use is not allowed.

% Default parameters
if ~exist('analysis_param','var') || isempty(analysis_param),
    analysis_param = 1;
end

if ~exist('kernel_param','var') || isempty(kernel_param)
    kernel_param = 1;
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
        [X_name]=scca(XN,KYN,analysis_param);
    case 'kcca'
        [X_name]=kcca(XN,KYN,analysis_param);
end

points_view1=[x y];
points_view2=[z w];

angle=[]; dist=[]; dist_1=[]; dist_2=[];
% compute cosines of the angles between the vectors
for m=1:length(points_view1)
    for n=1:length(points_view2)
        angle(m,n)=acosd(dot(points_view1(m,:),points_view2(n,:))/...
            (norm(points_view1(m,:),2)*norm(points_view2(n,:),2)));
        dist(m,n)=norm(points_view1(m,:)-points_view2(n,:));
        dist_1(m)=norm([0 0]-points_view1(m,:));
        dist_2(n)=norm([0 0]-points_view2(n,:));
    end
end

corrs = arrayfun(@(x) cosd(x), angle);
corrs(isnan(corrs)==1)=0;

% Discriminating between positive and negative correlations
weights=0.01:0.001:1;
scale=linspace(0,max(max(dist)),size(weights,2));
dist_scaled=zeros(size(dist));

for row=1:size(corrs,1)
    for col=1:size(corrs,2)
        % High positive correlation, great weight
        if corrs(row,col)>=0
            weights=sort(weights,'descend');
            [~,ind] = histc(dist(row,col),scale);
            dist_scaled(row,col)=weights(ind)*dist_1(row)*dist_2(col);%
            W(row,col)=weights(ind);
            % High negative correlation, great weight
        else
            weights=sort(weights,'ascend');
            [~,ind] = histc(dist(row,col),scale);
            dist_scaled(row,col)=weights(ind)*dist_1(row)*dist_2(col);%f
        end
    end
end

matrix=corrs.*dist_scaled;
matrix(isnan(matrix)==1)=0;

cgo = clustergram(matrix,'Standardize','None','colormap','redbluecmap');
set(cgo,'ColumnLabels',Yvar_name,'RowLabels',X_name)


    function [X_name]=scca(XN,KYN,analysis_param)
        rng(1,'v5normal')
        [seeds,~] = spectral_clustering(YN,projections);
        [w, beta, ~, ~, ~, ~]...
            = SCCAwrapper2_corrected_corr_cvx(XN',KYN,seeds,1,0,analysis_param,1);
        scores = XN*w;
        %scores = KYN*beta;
        
        % indices for the sparse features
        vec1=logical(w(:,pro1));
        vec2=logical(w(:,pro2));
        
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
    end

    function [X_name]=kcca(XN,KYN,analysis_param)
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
    end


end


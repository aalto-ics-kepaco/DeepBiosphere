function [seeds,S] = spectral_clustering(X,k)
% Spectral clustering algorithm by Ng, Jordan and Weiss(2002)

%[opt_width] = opt_similarity(X,k);
% Similarity matrix from data X by Gaussian kernel
A = gram( X', X', 'gaussian',1); % Gaussian kernel
A(logical(eye(size(A)))) = 0; % diagonal elements = 0

% calculate degree matrix
degs = sum(A, 2);
D    = sparse(1:size(A, 1), 1:size(A, 2), degs);

% compute unnormalized Laplacian
L = D - A;
  
% avoid dividing by zero
degs(degs == 0) = eps;
% calculate D^(-1/2)
D = spdiags(1./(degs.^0.5), 0, size(D, 1), size(D, 2));

% calculate normalized Laplacian
L = D * L * D;

% compute the eigenvectors corresponding to the k smallest
% eigenvalues
diff   = eps;
[U, ~] = eigs(L, k, diff);

% normalize the eigenvectors row-wise
U = bsxfun(@rdivide, U, sqrt(sum(U.^2, 2)));

for i=1:100
    [~,~,~,~,s] = kmedoids(U,k);
    S(i,:)=sort(s); 
end

% for i=1:100
%     centroids = find_centroids(U,k);
%     S(i,:)=sort(centroids);
% end

seeds=mode(S);



end


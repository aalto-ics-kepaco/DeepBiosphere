function AN = normalize_data(A,unitlen)

if nargin < 2
    unitlen = 0;
end

%standarization of the VARIABLES
AN = A - ones(size(A,1),1)*mean(A);
nonzero = std(A) > 1E-6;

AN(:,nonzero) = AN(:,nonzero)./(ones(size(A(:,nonzero),1),1)*std(A(:,nonzero)));


%creating unit vectors
if unitlen % if the user wishes so
    
invrownorms = 1./sqrt(sum(AN.*AN,2));
AN = diag(invrownorms)*AN;

       
end

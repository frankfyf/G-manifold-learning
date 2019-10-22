function [ affinity ] = MFVDM_ps_SO3_s( Evec, Eval, t)
% Using power spectrum method to compute the affinity with single frequency
% k = 1 in SO3

% compute the filtered matrix
n = size(Evec, 1)/3;
Eval = Eval.^t;
evec_norm = sqrt(sum(abs(Evec).^2, 2));
evec = bsxfun(@times, Evec, 1./evec_norm);
W = bsxfun(@times, evec, Eval.')*evec'; %denoising the matrices

% compute the affinity for each pair
affinity = zeros(n,n);
for i = 1:n
    for j = i+1:n
        tmp = W((i-1)*3+1:i*3, (j-1)*3+1:j*3);
        affinity(i,j) = trace(tmp*tmp');
        disp(['Finish the power spectrum affinity for nodes ', num2str(i), ' and ', num2str(j)]);
    end
end

affinity = affinity + affinity';



end


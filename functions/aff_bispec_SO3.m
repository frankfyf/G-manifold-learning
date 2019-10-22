% function [ affinity ] = aff_bispec_SO3(Evec, Eval, t)
% 
% Author:   Yifeng Fan (yifengf2@illinois.edu)
% Date:     2019/10/17   
% 
% Description: Compute the bispectrum affinity between nodes, based on 
% the eigen-decomposition. The transformation is SO(3)
% 
% Parameters : Evec              -- eigenvectors, a 1 by k_max cell
%              Eval              -- eigencalues, a 1 by k_max cell
%              t                 -- diffusion time
% 
% Return     : affinity          -- a n by n matrix, the affinity between n nodes.

function [ affinity ] = aff_bispec_SO3( Evec, Eval, t)

% Weighted matrices filtering 
k_max = size(Evec, 2);
n = size(Evec{1}, 1);
W = cell(1,k_max);
for i = 0:k_max-1
    eval = Eval{i+1}.^(2*t);
    evec = Evec{i+1};
    evec_norm = sqrt(sum(abs(evec).^2, 2));
    evec = bsxfun(@times, evec, 1./evec_norm);
    W{i+1} = bsxfun(@times, evec, eval.')*evec'; 
end

% Compute the Clebsch-Gorden matrices
affinity = zeros(n,n, k_max, k_max);
cg_all = cell(k_max, k_max); 
for sigma = 0:k_max-1
    for delta = sigma:k_max-1
        if min(sigma+delta,k_max-1) - abs(sigma-delta) >= 0 
            cg_all{sigma+1,delta+1} = cg_matrix(sigma, delta, k_max-1);
        end
    end
end

% Compute the bispectrum affinity for each entry
for i = 1:n
    for j = i+1:n
        count = 0;
        for sigma = 0:k_max-1
            for delta = sigma:k_max-1
                % Compute the Kronecker product of F_{sigma} and F_{delta}
                if sigma + delta <= k_max-1
                    tmp_kron = kron(W{sigma+1}((i-1)*(2*sigma+1)+1:i*(2*sigma+1),(j-1)*(2*sigma+1)+1:j*(2*sigma+1)), ...
                        W{delta+1}((i-1)*(2*delta+1)+1:i*(2*delta+1),(j-1)*(2*delta+1)+1:j*(2*delta+1)));
                    tmp_a = abs(sigma-delta);
                    tmp_b = sigma+delta;
                    tmp_mid = zeros((tmp_a+tmp_b+1)*(tmp_b - tmp_a +1), (tmp_a+tmp_b+1)*(tmp_b - tmp_a +1));
                    count = 0;
                    for kk = tmp_a:tmp_b
                        tmp_mid(count+1:count+(2*kk+1),count+1:count+(2*kk+1)) = (W{kk+1}((i-1)*(2*kk+1)+1:i*(2*kk+1), (j-1)*(2*kk+1)+1:j*(2*kk+1)))';
                        count = count + (2*kk+1);
                    end
                    affinity(i,j,sigma+1, delta+1) = trace(tmp_kron * cg_all{sigma+1, delta+1} * tmp_mid * (cg_all{sigma+1, delta+1}'));
                    count = count + 1;
                end
            end
        end
        % Print the result
        disp(['Finish the bispectrum affinity for nodes ', num2str(i), ' and ', num2str(j)]);
    end
end
affinity = sum(sum(affinity,3),4)/count;
affinity = affinity + affinity';

end


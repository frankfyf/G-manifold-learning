% function [ Eval, Evec ] = get_eigen(angle, list, eigen_num, n)
% 
% Author:   Yifeng Fan (yifengf2@illinois.edu)
% Date:     2019/10/17   
% 
% Description: Build a series of affinity matrices based on the edge
% connections, then do eigen-decomposions. The transformation is SO(2)
% 
% Parameters : angle             -- rotation on edges 
%              list              -- edge connection list  
%              eigen_num         -- truncation of eigenvalues
%              n                 -- number of nodes
% 
% Return     : Eval              -- eigenvalues, a 1 by k_max cell
%              Evec              -- eigenvectors, a 1 by k_max cell

function [ Eval, Evec ] = get_eigen(angle, list, eigen_num, n)

I = list(:,1);
J = list(:,2);
k_max = size(eigen_num,2);

%Compute the eigenvalues and eigenvectors for each frequency k
parfor k = 1:k_max
    W = sparse(I, J, exp(sqrt(-1)*(k)*angle), n, n); % should be a sparse matrix
    D = sum(abs(W), 2);
    W = bsxfun(@times, 1./sqrt(D), W);
    W = bsxfun(@times, 1./sqrt(D).', W);
    [ u, d ] = eigs(W, eigen_num(k), 'lm');
    [ sorted_eigval, id ] = sort(real(diag(d)), 'descend');
    sorted_eigvec = u(:, id);
    sorted_eigval(isnan(sorted_eigval)) = 0;
    sorted_eigvec(:,isnan(sorted_eigval)) = 0; 
%     figure
%     bar(sorted_eigval)
%     xlim([0 50])
%     set(gca,'FontSize',12);
%     title(['k = ', num2str(k)])
%     saveas(gcf,['figures/eig_',num2str(k),'_','.fig']);
    Evec{k} = sorted_eigvec;
    Eval{k} = sorted_eigval;
end

end


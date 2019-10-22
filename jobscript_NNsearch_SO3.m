% Name: jobscript_NNsearch_SO3
% 
% Author:   Yifeng Fan (yifengf2@illinois.edu)
% Date:     2019/10/21   
% 
% Description: Jobscript example performing nearest neighbor search on a
% synthetic dataset, specifically: M = SO(4), G = SO(3), B = S^3
%   (1) Generate data points on a S^3 sphere with a G = SO(3) transformation 
%       and apply the random rewiring model for adding noise. 
%   (2) Compute affinities
%   (3) Nearest neighbor search is done by sorting the affinities. 
% We compare to the baseline: vector diffusion maps (VDM). Note that at 
% this SO(3) case we don't implement the optimal alignment affinity due to 
% the complexity.
%
% Note: need to install Robotics System Toolbox.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% You may add these at first %%
clear
close all
addpath(genpath('./'))
rng('default');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters %%
disp('Preprocessing...')
n = 1000; % number of data points 
k_max = 6; % maximum frequency
m_k = 10; % eigenvalue truncation for each frequency
id_nn = 50;

% Uniformly sample on S^3
q = randn(n,4); 
q = diag(1./sqrt(sum(q.^2,2)))*q; % normaliztion
corr = q*q';

% Preallocating 
dis_bispec= zeros(n, id_nn);
dis_ps= zeros(n, id_nn);
dis_VDM= zeros(n, id_nn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main %%

p_range = [0.2, 0.3, 0.5, 1]; % range of random rewiring probabilities
for num = 1:numel(p_range)
    rng('default')
    p = p_range(num);
    disp(['The random rewiring probability is ', num2str(p)]) 
    
    % Random graph by kappa-nearest neighbor criteria
    kappa = 20;
    disp('Building the random graph (kappa-nearest neighbor)...');
    [ list, wigner_D ] = random_graph_S3_knn(p, corr, kappa, k_max);

    % Eigen-decomposition
    disp('Eigen-decomposition...')
    [ Evec, Eval ] = get_eigen_SO3( wigner_D, list, m_k*ones(1,k_max+1), n );  % get the eigenfucntions
    disp(['The average nearest neighbors of each node is: ', num2str(floor(size(list,1)/n))]);
    
    % Compute the affinities 
    % Bispectrum affinity
    disp('Computing the bispectrum affinity...');
    t = 1;
    [affinity_bispec] = aff_bispec_SO3(Evec, Eval, t);
    % Power spectrum affinity
    disp('Computing the power spectrum affinity...');
    t = 1;
    [affinity_ps] = aff_ps_SO3(Evec, Eval, t);
    % VDM affinity (VDM is a special case of power spectrum affinity with k = 1.
    disp('Computing the VDM affinity...');
    t = 1;
    [affinity_VDM] = aff_ps_SO3_s(Evec{2}, Eval{2} ,t);
    
    % Nearest neighbor identification based on sorting the affinity
    disp('Sorting affinity...');
    [~, id_bispec] = sort(affinity_bispec, 2, 'descend');
    class_bispec = id_bispec(:,2:id_nn+1); 
    [~, id_ps] = sort(affinity_ps, 2, 'descend');
    class_ps = id_ps(:,2:id_nn+1);
    [~, id_VDM] = sort(affinity_VDM, 2, 'descend');
    class_VDM = id_VDM(:,2:id_nn+1);
    
    % Computing the geodesic distance between identified neighbors
    for i = 1:n
        for j = 1:knn
            dis_bispec(i,j) = acos(corr(i,class_bispec(i,j)))*180/pi;
            dis_ps(i,j) = acos(corr(i,class_ps(i,j)))*180/pi;
            dis_VDM(i,j) = acos(corr(i,class_VDM(i,j)))*180/pi;
        end
    end
    
    % Plot the histogram nearest neighbor search result
    [ps_x, ps_y] = hist(dis_ps(:), 0:3:180);
    [VDM_x, VDM_y] = hist(dis_VDM(:), 0:3:180);
    [bispec_x, bisp_y] = hist(dis_bispec(:), 0:3:180);
    figure
    plot(ps_y, ps_x/sum(ps_x), bispec_y, bispec_x/sum(bispec_x), VDM_y, VDM_x/sum(VDM_x));
    legend('Power spec. (ours)', 'Bispec. (ours)', 'VDM');
    set(gca, 'fontsize', 16);
    set(gca, 'XGrid', 'on');
    set(gca, 'YGrid', 'on');
    xlim([1,180])
    hLegend = findobj(gcf, 'Type', 'Legend');
    hLegend.FontSize = 20;
    hLegend.Box = 'off';
    xlabel('Geodesic distance between neighbors', 'Fontsize', 18, 'Interpreter', 'latex');
    ylabel('Proportion of neighbors', 'Fontsize', 18, 'Interpreter', 'latex');
    title(['p = ',num2str(p)]);
end

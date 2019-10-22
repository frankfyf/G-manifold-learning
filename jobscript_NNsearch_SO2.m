% Name: jobscript_NNsearch_SO2
% 
% Author:   Yifeng Fan (yifengf2@illinois.edu)
% Date:     2019/10/16   
% 
% Description: Jobscript example performing nearest neighbor search on a
% synthetic dataset, specifically: M = SO(3), G = SO(2), B = S^2
%   (1) Generate data points on a S^2 sphere with a G = SO(2) transformation 
%       and apply the random rewiring model for adding noise. 
%   (2) Compute affinities
%   (3) Nearest neighbor search is done by sorting the affinities. 
% We compare to the baseline: vector diffusion maps (VDM)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% You may add these at first %%
clear 
close all
addpath(genpath('./'));
rng('default') 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters %%
disp('Preprocessing...')
n = 10000; % number of data points 
q = qrand(n); % quaternions
rot = q_to_rot(q); % rotation matrices 
v = squeeze(rot(3, :, :)); % viewing directions of data points 
clear rot;
corr = v'*v; 

k_max = 10; % maximum frequency
m_k = 20; % eigenvalue truncation for each frequency
id_nn = 50;

% Euclidean and geodesic distance between data points 
e_dis = zeros(n,n);
for i = 1:n
    for j = i+1:n
        e_dis(i,j) = norm(v(:,i)-v(:,j),2);
    end
end
e_dis = e_dis + e_dis'; % Euclidean distance
g_dis = abs(acos(v'*v)); % Geodesic distance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main %%

p_range = [0.08, 0.09, 0.1, 0.5]; % range of random rewiring probabilities
for num = 1:numel(p_range)
    rng('default')
    p = p_range(num);
    disp(['The random rewiring probability is ', num2str(p)]) 
    
    % We have two standards for building the neighborhood graph
    % Random graph by kappa-nearest neighbor criteria
    disp('Building the random graph (kappa-nearest neighbor)...');
    kappa = 150; 
    [list, angle] = random_graph_knn(corr,kappa,p,q); 
    disp(['The average nearest neighbors of each node is: ', num2str(floor(size(list,1)/n))]);
    
%     % Random graph by epsilon-neighborhood criteria
%     disp(['Building the random graph (epsilon-neighborhood)...']);
%     epsilon = 0.97; 
%     [list, angle] = random_graph_epsilon(corr,epsilon,p,q); 
%     disp(['The average nearest neighbors of each node is: ', num2str(floor(size(list,1)/n))]);
    
    % Build the affinity matrices and eigen-decomposition
    disp('Eigen-decompostion...');
    [Eval, Evec] = get_eigen(angle,list,m_k*ones(1,k_max),n); 
    
    % Compute the affinities 
    % Optimal alignment affinity
    disp('Computing the optimal alignment affinity...');
    t = 1; % diffusion time
    [affinity_opt] = aff_opt(Evec, Eval, t);
    % Bispectrum affinity
    disp('Computing the bispectrum affinity...');
    t = 1;
    [affinity_bispec] = aff_bispec(Evec, Eval, t);
    % Power spectrum affinity
    disp('Computing the power spectrum affinity...');
    t = 1;
    [affinity_ps] = aff_ps(Evec, Eval, t);
    % VDM affinity (VDM is a special case of power spectrum affinity with k = 1.
    disp('Computing the VDM affinity...');
    t = 1;
    [affinity_VDM] = aff_ps(Evec(1), Eval(1), t);
    
    % Nearest neighbor identification based on sorting the affinity
    disp('Sorting affinity...');
    [~, id_bispec] = sort(affinity_bispec, 2, 'descend');
    class_bispec = id_bispec(:,2:id_nn+1); 
    [~, id_ps] = sort(affinity_ps, 2, 'descend');
    class_ps = id_ps(:,2:id_nn+1);
    [~, id_opt] = sort(affinity_opt, 2, 'descend');
    class_opt = id_opt(:,2:id_nn+1);
    [~, id_VDM] = sort(affinity_VDM, 2, 'descend');
    class_VDM = id_VDM(:,2:id_nn+1);
    
    % Compute the angles between viewing directions of identified neighbors
    disp('Checking nearest neighbor search result...');
    [ e_c_bispec, ~ ] = check_simulation_results(class_bispec, ones(size(class_bispec)), ones(size(class_bispec)), q);
    [ e_c_ps, ~ ] = check_simulation_results(class_ps, ones(size(class_ps)), ones(size(class_ps)), q);
    [ e_c_opt, ~ ] = check_simulation_results(class_opt, ones(size(class_opt)), ones(size(class_opt)), q);       
    [ e_c_VDM, ~ ] = check_simulation_results(class_VDM, ones(size(class_VDM)), ones(size(class_VDM)), q);
    
    % Plot the histogram nearest neighbor search result
    [y_c_bispec,x_c_bispec] = hist(acos(e_c_bispec(:))*180/pi,0:180);
    [y_c_ps,x_c_ps] = hist(acos(e_c_ps(:))*180/pi,0:180);
    [y_c_opt,x_c_opt] = hist(acos(e_c_opt(:))*180/pi,0:180);
    [y_c_VDM,x_c_VDM] = hist(acos(e_c_VDM(:))*180/pi,0:180);
    
    figure
    plot(x_c_ps, y_c_ps/sum(y_c_ps), x_c_bispec, y_c_bispec/sum(y_c_bispec), ...
         x_c_opt, y_c_opt/sum(y_c_opt), x_c_VDM, y_c_VDM/sum(y_c_VDM), 'linewidth', 3)
    legend('Power spec. (ours)', 'Bispec. (ours)', 'Opt. (ours)', 'VDM');
    set(gca, 'fontsize', 16);
    set(gca, 'XGrid', 'on');
    set(gca, 'YGrid', 'on');
    xlim([1,180])
    hLegend = findobj(gcf, 'Type', 'Legend');
    hLegend.FontSize = 20;
    hLegend.Box = 'off';
    xlabel('arccos$\left\langle v_i, v_j \right\rangle$', 'Fontsize', 18, 'Interpreter', 'latex');
    ylabel('Proportion of neighbors', 'Fontsize', 18, 'Interpreter', 'latex');
    title(['p = ',num2str(p)]);
end

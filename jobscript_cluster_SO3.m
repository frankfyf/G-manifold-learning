% Name: jobscript_cluster_SO3
% 
% Author:   Yifeng Fan (yifengf2@illinois.edu)
% Date:     2019/10/17   
% 
% Description: Jobscript example performing spectral clustering on
% synthetic dataset, transformation is SO(3). Specifically:
%   (1) Generate data points with k clusters and add noise based on the 
%       random rewiring model
%   (2) Compute affinities 
%   (3) Compute embedding based on the affinities 
%   (4) K-means for clustering
% We compare with vector diffusion maps (VDM) and the scaler spectral 
% clustering (Andrew Ng, 2002). Note that at this SO(3) case we don't
% implement the optimal alignment affinity due to the complexity.
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
K = 10; % number of clusters
k_max = 10; % maximum frequency 
m_k = K; % truncation of eigenvalues

n_points = 25*ones(1,K); % number of points in each cluster
n = sum(n_points); 
%tmp_index = randperm(n);
tmp_index = [1:n];
index = cell(1,K); % the index of points for each cluster
tmp = 0;
id_real = zeros(1,n);
for i = 1:K
    index{1,i} = tmp_index(tmp+1:tmp+n_points(i));
    id_real(1,tmp_index(tmp+1:tmp+n_points(i))) = i; % assign the cluster index
    tmp = tmp+n_points(i);
end

% The following is for visualizing the clustering. We will generate
% 2 dimension coordinates (x,y) for data points whose positions are close if in the same
% cluster. 
plot_flag = 0; % show plots (Y = 1/N = 0)?
coor = zeros(n,2);
shift = 2;
r_c = 3;
tmp_angle = 2*pi*rand(n,1);
tmp_radius = rand(n,1);
tmp_angle_c = 2*pi*[1:K]'/K;
tmp_radius_c = r_c*rand(K,1) + shift;

% We perform multiple trials
trial = 10; % number of trials
p_range = [0.16,0.20,0.25]; % range of random rewiring probabilities

% Preallocating
rand_index_ps = zeros(trial, numel(p_range));
rand_index_bispec = zeros(trial, numel(p_range));
rand_index_VDM = zeros(trial, numel(p_range));
rand_index_scalar = zeros(trial, numel(p_range));
id_ps = cell(trial, numel(p_range));
id_bispec = cell(trial, numel(p_range));
id_VDM = cell(trial, numel(p_range));
id_scalar = cell(trial, numel(p_range));
affinity_ps = cell(trial, numel(p_range));
affinity_bispec = cell(trial, numel(p_range));
affinity_VDM = cell(trial, numel(p_range));
affinity_scalar = cell(trial, numel(p_range));

for num = 1:numel(p_range)
    p = p_range(num);
    disp(['-------- The random rewiring probability is ', num2str(p), ' --------']); 
    rng('default')

    for time = 1:trial
        disp(['------------------ The trial ', num2str(time), ' starts ------------------']);
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main %%
        
        % Generate edge connection list and transformations
        disp('Building the random graph...');
        [ list, wigner_D ] = k_cluster_graph_SO3(p, index, k_max);

        % Eigen-decompostion
        disp('Eigen-decomposition...')
        [ Evec, Eval ] = get_eigen_SO3( wigner_D, list, m_k*ones(1,k_max+1), n ); 
        
        % Compute affinity
        disp('Computing power spectrum affinity...')
        t = 1;
        [ affinity_ps{time, num} ] = aff_ps_SO3( Evec, Eval, t );
        affinity_ps{time, num} = affinity_ps{time, num} - diag(diag(affinity_ps{time, num}));
        disp('Computing bispectrum affinity...')
        t = 1;
        [ affinity_bispec{time, num} ] = aff_bispec_SO3( Evec, Eval, t );
        affinity_bispec{time, num} = affinity_bispec{time, num} - diag(diag(affinity_bispec{time, num}));
        disp('Computing VDM affinity...')
        t = 1;
        [ affinity_VDM{time, num} ] = aff_ps_SO3_s( Evec{2}, Eval{2} ,t );
        affinity_VDM{time, num} = affinity_VDM{time, num} - diag(diag(affinity_VDM{time, num}));
        
        disp('Computing clustering embedding...')
        % Bispectrum affinity
        [Evec_bispec,~] = eigs(abs(affinity_bispec{time, num}),K);
        tmp_norm = sqrt(sum(abs(Evec_bispec.^2),2));
        Evec_bispec = bsxfun(@times, Evec_bispec, 1./tmp_norm);
        % Power spectrum affinity
        [Evec_ps,~] = eigs(abs(affinity_ps{time, num}),K);
        tmp_norm = sqrt(sum(abs(Evec_ps.^2),2));
        Evec_ps = bsxfun(@times, Evec_ps, 1./tmp_norm);
        % VDM affinity
        [Evec_VDM,~] = eigs(abs(affinity_VDM{time, num}),K);
        tmp_norm = sqrt(sum(abs(Evec_VDM.^2),2));
        Evec_VDM = bsxfun(@times, Evec_VDM, 1./tmp_norm);
        
        %% Clustering by the old method
        disp('Scalar sepctral clustering (Andrew Ng, 2002)...');
        affinity_scalar{time, num} = sparse(list(:,1), list(:,2), 1, n, n);
        [~,Evec_scalar] = get_eigen( zeros(size(list,1),1), list, K*ones(1,1), n);
        Evec_scalar = cell2mat(Evec_scalar);
        tmp_norm = sqrt(sum(abs(Evec_scalar.^2),2));
        Evec_scalar = bsxfun(@times, Evec_scalar, 1./tmp_norm);

        % Perform the k-means on the affinity for clustering
        disp('K-means for clustering...');
        id_bispec{time, num} = kmeans(Evec_bispec,K);
        id_ps{time, num} = kmeans(Evec_ps,K);
        id_VDM{time, num} = kmeans(Evec_VDM,K);
        id_scalar{time, num} = kmeans(Evec_scalar,K);

        % Evaluation: rand index
        rand_index_ps(time, num) = rand_index(id_ps{time, num}, id_real);
        rand_index_bispec(time, num) = rand_index(id_bispec{time, num}, id_real);
        rand_index_VDM(time, num) = rand_index(id_VDM{time, num}, id_real);
        rand_index_scalar(time, num) = rand_index(id_scalar{time, num}, id_real);
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visualization %%

    % Visualization for the clusters, the result comes from the first trial
    if plot_flag
        for i = 1:K
            % Generate the coordinates for data points
            coor(index{i},:) = [tmp_radius_c(i,1)*cos(tmp_angle_c(i,1)),tmp_radius_c(i,1)*sin(tmp_angle_c(i,1))] +...
                [tmp_radius(index{i}).*cos(tmp_angle(index{i})), tmp_radius(index{i}).*sin(tmp_angle(index{i}))];  
        end
        
        tmp_id = {id_ps{1,num}, id_bispec{1,num}, id_VDM{1,num}, id_scalar{1,num}};
        tmp_aff = {affinity_ps{1,num}, affinity_bispec{1,num}, affinity_VDM{1,num}, affinity_scalar{1,num}};
        tmp_name = {'Power spectrum', 'Bispectrum', 'VDM' , 'Scalar'};
        
        % Plot the clustering result
        for i = 1:numel(tmp_name)
            figure 
            scatter(coor(:,2), coor(:,1), 300, tmp_id{i}, '.')
            title(['p = ', num2str(p), ', ', tmp_name{i}], 'Fontsize', 10);
            axis equal
            xlim([-(r_c+shift+1),(r_c+shift+1)])
            ylim([-(r_c+shift+1),(r_c+shift+1)])
            xlabel('x', 'Fontsize', 16);
            ylabel('y', 'Fontsize', 16);
            set(gca, 'fontsize', 16);
            set(gca, 'XGrid', 'on');
            set(gca, 'YGrid', 'on')
        end
        
        % Plot the affinity matrix 
        for i = 1:numel(tmp_name)        
            figure
            imagesc(abs(tmp_aff{i}));
            title(['p = ', num2str(p), ', ', tmp_name{i}], 'Fontsize', 10);
            axis equal
            axis off
            colorbar
        end
    end    
end

% Compute the average rand index
p_range
rand_index_ps_mean = mean(rand_index_ps,1)
rand_index_bispec_mean = mean(rand_index_bispec,1)
rand_index_VDM_mean = mean(rand_index_VDM,1)
rand_index_scalar_mean = mean(rand_index_scalar,1)



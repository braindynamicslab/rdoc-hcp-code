% Calculating modularity for Mapper analysis of subjects in cohort 2. This code sets G1 to NVS correctly

% Set project root
projectRoot = '';

% Path to BDL code for Mapper
addpath(fullfile(projectRoot, 'bdl-mapper-v2-main'));

% Load data for cohort 2
load(fullfile(projectRoot, 'comm_centroid2.mat'));            % 5×347
load(fullfile(projectRoot, 'contrasts_all_2d_sub2.mat'));     % M×347
contrasts_all_2d_sub2 = data;
load(fullfile(projectRoot, 'row_labels_sub2.mat'));           % cell array of M
row_labels_sub2 = rowNames;

% Get number of subjects
n_maps_per_subject = 24;
n_subjects2 = size(contrasts_all_2d_sub2, 1) / n_maps_per_subject;

% Jackknife setup
mod_rdoc_sub2 = zeros(n_subjects2, 1);
mod_centroid_sub2 = zeros(n_subjects2, 1);

% Get domain indices for rdoc annotation
labels2 = {'NV-E1', 'CS-E2', 'NV-G1', 'PV-G2', 'CS-L1', 'CS-L2', 'CS-M1', 'SS-M2', 'SS-M3', 'SS-M4', 'SS-M5', 'SS-M6', 'CS-R1', 'CS-R2', 'CS-S1', 'SP-S2', 'CS-W1', 'CS-W2', 'CS-W3', 'CS-W4', 'CS-W5', 'CS-W6', 'CS-W7', 'CS-W8'}; % 24 task contrast labels
domains_rdoc = cellfun(@(x) x(1:2), labels2, 'UniformOutput', false);
[~, ~, domain_indices_rdoc] = unique(domains_rdoc);
metricType = 'euclidean'; ndim = 2; conn = .95;

% Parallel processing
parpool('local', 8);

parfor i = 1:n_subjects2
    leaveout = i;
    keep_idx = setdiff(1:n_subjects2, leaveout);

    start_idx = (leaveout - 1) * n_maps_per_subject + 1;
    end_idx = leaveout * n_maps_per_subject;
    selected_data2 = contrasts_all_2d_sub2;
    selected_data2(start_idx:end_idx, :) = [];

    selected_comm_rdoc2 = repmat(domain_indices_rdoc(:), length(keep_idx), 1);

    corrmat2 = corr(selected_data2', comm_centroid2');
    [~, comm_centroid2_labels] = max(corrmat2, [], 2);

    num_k = 88; res_val = 26; gain_val = 65; num_bin_clusters = 10;
    [nodeTpMat2, nodeBynode2, ~, ~, ~, ~] = runBDLMapper(selected_data2, metricType, res_val, gain_val, num_k, num_bin_clusters, ndim);

    [nodeComm_rdoc2, ~] = computeProportionCategorical(nodeTpMat2, selected_comm_rdoc2);
    mod_rdoc_sub2(i) = calMod(nodeBynode2, nodeComm_rdoc2);
    [nodeComm_comm2, ~] = computeProportionCategorical(nodeTpMat2, comm_centroid2_labels);
    mod_centroid_sub2(i) = calMod(nodeBynode2, nodeComm_comm2);
end

save(fullfile(projectRoot, 'modularity_jackknife_results_NVG1_sub2.mat'), ...
     'mod_rdoc_sub2', 'mod_centroid_sub2');
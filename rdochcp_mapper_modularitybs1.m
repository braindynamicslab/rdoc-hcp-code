% Calculating modularity for Mapper analysis of subjects in cohort 1. This code sets G1 to NVS correctly

% Set project root
projectRoot = '';

% Path to BDL code for Mapper
addpath(fullfile(projectRoot, 'bdl-mapper-v2-main'));

% Load data for cohort 1
load(fullfile(projectRoot, 'comm_centroid1.mat'));            % 4×347
load(fullfile(projectRoot, 'contrasts_all_2d_sub1.mat'));     % N×347
contrasts_all_2d_sub1 = data;
load(fullfile(projectRoot, 'row_labels_sub1.mat'));           % cell array of N
row_labels_sub1 = rowNames;

% Get number of subjects
n_maps_per_subject = 24;
n_subjects1 = size(contrasts_all_2d_sub1, 1) / n_maps_per_subject;

% Jackknife setup
mod_rdoc_sub1 = zeros(n_subjects1, 1);
mod_centroid_sub1 = zeros(n_subjects1, 1);

% Get domain indices for rdoc annotation
labels2 = {'NV-E1', 'CS-E2', 'NV-G1', 'PV-G2', 'CS-L1', 'CS-L2', 'CS-M1', 'SS-M2', 'SS-M3', 'SS-M4', 'SS-M5', 'SS-M6', 'CS-R1', 'CS-R2', 'CS-S1', 'SP-S2', 'CS-W1', 'CS-W2', 'CS-W3', 'CS-W4', 'CS-W5', 'CS-W6', 'CS-W7', 'CS-W8'}; % 24 task contrast labels
domains_rdoc = cellfun(@(x) x(1:2), labels2, 'UniformOutput', false);
[~, ~, domain_indices_rdoc] = unique(domains_rdoc);
metricType = 'euclidean'; ndim = 2; conn = .95;

% Parallel processing
parpool('local', 8);

parfor i = 1:n_subjects1
    leaveout = i;
    keep_idx = setdiff(1:n_subjects1, leaveout);

    start_idx = (leaveout - 1) * n_maps_per_subject + 1;
    end_idx = leaveout * n_maps_per_subject;
    selected_data1 = contrasts_all_2d_sub1;
    selected_data1(start_idx:end_idx, :) = [];

    selected_comm_rdoc1 = repmat(domain_indices_rdoc(:), length(keep_idx), 1);

    corrmat1 = corr(selected_data1', comm_centroid1');
    [~, comm_centroid1_labels] = max(corrmat1, [], 2);

    num_k = 99; res_val = 26; gain_val = 65; num_bin_clusters = 10;
    [nodeTpMat1, nodeBynode1, ~, ~, ~, ~] = runBDLMapper(selected_data1, metricType, res_val, gain_val, num_k, num_bin_clusters, ndim);

    [nodeComm_rdoc1, ~] = computeProportionCategorical(nodeTpMat1, selected_comm_rdoc1);
    mod_rdoc_sub1(i) = calMod(nodeBynode1, nodeComm_rdoc1);
    [nodeComm_comm1, ~] = computeProportionCategorical(nodeTpMat1, comm_centroid1_labels);
    mod_centroid_sub1(i) = calMod(nodeBynode1, nodeComm_comm1);
end

save(fullfile(projectRoot, 'modularity_jackknife_results_NVG1.mat'), ...
     'mod_rdoc_sub1', 'mod_centroid_sub1');
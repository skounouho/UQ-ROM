
addpath('./../matlab');

% % apply the zero Dirichlet boundary condition
% bcidx = [1:34:239, (1:34:239)+272];
% rob(bcidx, :, :) = 0;

% run in parallel
% workers = 16;
pool = parpool;

% turn off repetitive warning all all processors
parfor i = 1:pool.NumWorkers
    warning('off','MATLAB:logm:nonPosRealEig');
end

% define constants
N = 86370; % natoms * 3 DOF
n = 20; % order of samples
k = 5; % number of models (including global)
s = 500; % number of samples

% load rob
rob = zeros(N, n, k);
folders = ["adp" "bop" "eam" "meam" "global"];
parfor i = 1:size(rob, 3)
    filename = sprintf("%s/seg.rob\n", folders(i));
    opts = detectImportOptions(filename, 'FileType', 'text', 'Range', [1 1 N n], 'Delimiter', " ");
    M = readmatrix(filename, opts);
    [row,col] = size(M);
    temp = zeros(N, n);
    temp(:, 1:col) = M;

    rob(:, :, i) = temp;
end


%% project on the tangential space
U0 = rob(:, :, end);
% assume that the Log is well-defined on all the points
Vs = zeros(size(rob)); % projected samples in tangent space
tau = 1e-4;

parfor i = 1:k
    Vs(:,:,i) = real(stiefel_log(U0, rob(:,:,i), tau));
end

%% solve the quadratic programming problem
num_models = size(rob,3)-1;
X = reshape(Vs(:,:,1:num_models), [], num_models)';
H = X*X'; f = zeros(num_models,1); Aeq = ones(1,num_models); beq = 1; 
lb=zeros(1,num_models); ub=ones(1,num_models);
beta = quadprog(H, f, [], [], Aeq, beq, lb, ub);
beta = beta;
fprintf('beta = %s\n', join(string(beta), ','))

%% Dirichlet sampling 
rng(41);
N_samples = s;
w = gamrnd(repmat(beta',N_samples,1), 1, N_samples, length(beta));
w = w ./ sum(w,2);

% sorts so that BOP simulations come first
w = sortrows(w, 2, 'descend');

tangential_samples = w*X;
[~, I] = sort(w, 2, "ascend");
maxI = I(:,end);

%% project back to Stiefel manifold
tangential_samples = tangential_samples';
tangential_samples = reshape(tangential_samples, [N,n,N_samples]);
stiefel_samples = zeros(size(tangential_samples));
parfor i = 1:N_samples
    delta = tangential_samples(:,:,i);
    stiefel_samples(:,:,i) = stiefel_exp(U0, delta);
end

%% Frechet mean
% frechet_mean = calc_frechet_mean_mat(stiefel_samples, U0, 0.1);

%% save samples
p = round(floor(log10(s))) + 1;

parfor i = 1:s
    filename = sprintf("%0" + num2str(p) + "d/seg.rob", i);
    writematrix(stiefel_samples(:,:,i), filename, 'FileType', 'text', 'Delimiter', 'space');
end

%% save active selection
selection = string(zeros(size(maxI)));
for i = 1:s
    selection(i) = folders(maxI(i));
end
writematrix(selection, "selection.txt");

%% distance in canonical metric
% load data/stiefel_samples_2.5k_qp.mat
% U0 = stiefel_samples(:,:,k);
% frechet_mean = stiefel_samples(:,:,8);
% d = calc_dist_metric_cano(U0, frechet_mean);


% end parallel
delete(pool);

% turn on repetitive warning
warning('on','MATLAB:logm:nonPosRealEig')

exit;
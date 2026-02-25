% Clear workspace, close figures, close files, and clear command window
clear all; clc; fclose all;
global totalRunCount edata_with_Fpre_step ForwardCount last_time nMaterial prestress_time

% Ensure no parallel pool is open (good practice if using UseParallel later)
if ~isempty(gcp('nocreate'))
    delete(gcp('nocreate'));
end

%% --- Set up file paths and constants ---

path = struct();
% Get the current working directory 
path.parent = pwd; 

% Define the path to the data/src/results/VF folders
path.data = fullfile(path.parent, 'data');
path.src = fullfile(path.parent, 'src');
path.results = fullfile(path.parent, 'results');
path.VF = fullfile(path.parent, 'VF');

% FIle names
mymodel = 'Real_data_eye_fiber_choroid_with_pre_cyncl_HGO_4regions.feb';      % FE model file
myexpdata = 'Real_data_eye_fiber_choroid_with_pre_cyncl_HGO_4regions.log';    % Experimental data file
matFile = 'Real_data_eye_fiber_choroid_with_pre_cyncl_HGO_4regions.mat'; % Pre-Saved model data

%p_app = [0.15,0.02];                       % Reference and applied pressure
p_app=[-2.0;0.0];

% Last physical time at which prestress is calculated/frozen
prestress_time = 1.0;
%prestress_time = 0.0;

% Last Time of simulation
%last_time = 1.5;
last_time = 2;

% Penalty Factor for the contact
eps = 1000;

% Number of  different material models
nMaterial = 11;      

%--------------------------------------------------------------------------
% Gauss Order for Numerical Integration
%--------------------------------------------------------------------------
% Gauss order
gauss_order = 2;       % Selects 2x2x2 Gauss quadrature for elements

% Bounds for optimization variables [c1, c2]

lb = [0.6, 0.6,0.6, 0.6,0.6, 0.6,0.6, 0.6];   % lower bounds
ub = [1.4, 1.4,1.4, 1.4,1.4, 1.4,1.4, 1.4];   % upper bounds
Normalizer = [50,100,3266,8.62,172.4,308,5638,11.4];
corresponding = [1,2,2,3,3,4,4,4];
parameter = {'c1','c1','k','c1','k','c','k1','k2'};

%lb = [0.6, 0.6];   % lower bounds
%ub = [1.4, 1.4];   % upper bounds
%Normalizer = [0.1,0.5];
%corresponding = [1,2];
%parameter = {'c1','c1'};

count_corresponding = zeros(size(corresponding));
is_unique = zeros(size(corresponding)); 

for i = 1:length(corresponding)
    count_corresponding(i) = sum(corresponding(1:i) ==corresponding(i));
    is_unique(i) = sum(corresponding == corresponding(i)) == 1;
end

changing_matrix = [num2cell(corresponding); num2cell(count_corresponding);...
    num2cell(is_unique);parameter];


nvars = numel(lb);  % number of variables

%% --- Operations for copying and combining parameters in the parameter matrix ---
target_rows = {[5], [6], [7], [8]};

source_rows = {4, 4, 4, 4};
  
% each matches source_rows index 
source_cols = { [1 2 3 4], [1 2 3 4] ,[1 2 3 4], [1 2]}; 

target_cols = {[1 2 3 4], [1 2 3 4] ,[1 2 3 4], [1 1]};

% direct copies for 1:4
weights = {[1 1 1 1], [1 1 1 1] ,[1 1 1 1],[1 0.0031]};   

% Generating the ops struct
mat_size=[nMaterial,7];
ops_matrix_struct = create_ops(target_rows, target_cols, source_rows, source_cols, weights,mat_size);

%% --- Set fmincon optimization options ---
options = optimoptions('fmincon', ...
    'Display', 'iter', ...
    'Algorithm', 'interior-point', ...
    'FiniteDifferenceType', 'central', ...
    'UseParallel', false, ... % Parallel disabled
    'MaxFunctionEvaluations', 250, ...
    'StepTolerance', 1e-3);  

%% --- Generate starting points using Latin Hypercube Sampling ---
% randomize based on the ID
task_id = str2double(getenv('SLURM_ARRAY_TASK_ID'));
if isnan(task_id)
    % This isn't a job array, fallback to unique seed from JOBID or time:
    jobid = str2double(getenv('SLURM_JOB_ID'));
    if isnan(jobid)
        % Last resort: use time-based seed
        rng('shuffle');
    else
        rng(jobid);
    end
else
    rng(task_id);
end


N = 7; % number of random additional points
lhs_points = lhsdesign(N, length(Normalizer));  % generates points in [0,1]

% Scale lhs_points to actual bounds
start_points = bsxfun(@plus, lb, bsxfun(@times, lhs_points, (ub - lb)));

% Randomize the row order
randomized_order = randperm(size(start_points, 1));
start_points = start_points(randomized_order, :);

%start_points(1,:) = [1,1]; 

% Prepend custom start point
totalRunCount = 1; % Starts at 2 to skip calculation of the virtual field
ForwardCount = 1;

%% --- Prepare arrays to hold results ---
n_start = size(start_points, 1);       % number of start points
run_times = zeros(n_start, 1);         % time taken for each run
x_opts = zeros(n_start, nvars);        % each run's optimum x values
fvals = zeros(n_start, 1);             % each run's minimum cost
exitflags = zeros(n_start, 1);         % fmincon exit flags
outputs = cell(n_start, 1);            % fmincon output structs
model = struct();                      % model data (nodes, elements, surfaces, etc.)
edata = struct();                      % simulation data (experimental data or results)

%% --- Load mesh and experimental data --------------------------------------
%           Read Model and Node/Element Data, Compute Cumulative Fpre
%  --------------------------------------------------------------------------

% Loads mesh and log file, including node coordinates, element connectivity, results at timesteps
fullMatPath = fullfile(path.data, matFile);
if isfile(fullMatPath)
    % If the results file exists, load the data
    loaded = load(fullMatPath);           % Loads to a struct 'loaded'
    model = loaded.model;
    edata = loaded.edata;
    disp('Loaded model and edata from saved .mat file.');
else
    [model, edata] = preliminary_reading(path.data, mymodel, myexpdata,edata,model,prestress_time);
    save(fullMatPath, 'model', 'edata');
    disp('Read experiment dataparsed and saved results.');
end;

%% --- Define the cost function (anonymous wrapper) ---
cost_function = @(x) get_cost2regions_calc_Fpre(...
    path, mymodel, model, edata, x, p_app, gauss_order, prestress_time,eps,...
    changing_matrix,Normalizer,ops_matrix_struct);

%% --- Perform optimization separately at each starting point ---
start_points(1,:)
t=tic; % Start timing
for i = 1:n_start
    fprintf('Starting optimization %d/%d...\n', i, n_start);

    %Define optimization problem for this starting point
    problem = createOptimProblem('fmincon', ...
        'objective', cost_function, ...
        'x0', start_points(i,:), ... % use this start point
        'lb', lb, ...
        'ub', ub, ...
        'options', options);

    % Run optimization
    [x_opt, fval, exitflag, output] = fmincon(problem.objective, problem.x0, [], [],...
        [], [], lb, ub, [], options); 
    
    % Store results in corresponding arrays
    x_opts(i,:) = x_opt;
    fvals(i) = fval;

    %Refreshing counters
    totalRunCount = 2;
    ForwardCount = 1;

    x_run = x_opt .* Normalizer    
    start_points(i+1,:) = x_opt;
end

t_elapsed = toc(t); % End timing

%% --- Print optimized parameters and cost ---
disp('Optimized material parameters:');
disp(start_points);
fprintf('Final cost: %.6e\n', fvals(n_start));

%% --- Print convergence time for each starting point with its cost ---
disp('Convergence time for each starting point (seconds):');
fprintf('%d sec\n',t_elapsed)
for i = 1:n_start
    fprintf('Start %d: x0=[%.6f, %.6f, %.6f, %.6f, %.6f] --> cost=%.6e; x=[%.6f, %.6f, %.6f, %.6f, %.6f]\n', ...
        i, start_points(i,1), start_points(i,2),start_points(i,3),start_points(i,4),start_points(i,5),...
        fvals(i),x_opts(i,1), x_opts(i,2), x_opts(i,3), x_opts(i,4), x_opts(i,5) );
end


% Creating out file
save_optimization_results(path.results, t_elapsed, start_points, fvals, x_opts);
% Set working directory from command line
cd(data_path);

% Load data and extract interesting blocks
load('./data/EEM_NMR_LCMS.mat');
X = Y.('data');
Z = Z.('data');

% Normalize data
X = X / sqrt(sum(X(:).^2));
Z = normalize(Z, 'scale');
Z = Z / sqrt(sum(Z(:).^2));

% Call CMTF
T = struct;
T.object = { tensor(X), tensor(Z) };
T.size = [28 13324 8 168];
T.modes = { [1, 2, 3], [1, 4] };
R = 6;
n_repeat = 100;

options              = ncg('defaults');
% options.Display      ='final';
options.MaxFuncEvals = 100000;
options.MaxIters     = 10000;
options.StopTol      = 1e-12;

beta      = 1e-3;

%%% Save all models %%%
% Save CMTF
for i = 1:n_repeat
    tic
    [Zhat, G, output] = acmtf_opt(T, R, 'alg_options', options, 'beta_cp', 0, 'beta_pca', 0, 'init', 'random');
    time = toc;
    save(strcat('./models/CMTF/', string(i), '.mat'), 'Zhat');
    save(strcat('./models/CMTF/times/', string(i), '.mat'), 'time');
end

% Save ACMTF
for i = 1:n_repeat
    tic
    [Zhat, G, output] = acmtf_opt(T, R, 'alg_options', options, 'beta_cp', beta, 'beta_pca', beta, 'init', 'random');
    time = toc;
    save(strcat('./models/ACMTF/', string(i), '.mat'), 'Zhat');
    save(strcat('./models/ACMTF/times/', string(i), '.mat'), 'time');
end

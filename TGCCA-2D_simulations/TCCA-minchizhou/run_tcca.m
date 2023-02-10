try
    % Load parameters
    load(strcat(data_path, "/params.mat"));
    
    rng(seed);

    % Call TCCA with parameters
    [betax, betay, rho] = tcca(X, Y, rx, ry, "tau", tau, "covstr", covstr, "replicates", replicates);
    
    scorex = NaN;
    scorey = NaN;
    
    if compute_score
        ax = arrayfun(@(x) kron(betax{2}(:, x), betax{1}(:, x)), 1:rx, 'UniformOutput', false);
        ax = cellfun(@(x) x(:) / norm(x(:), 2), ax, 'UniformOutput', false);
        ax = cat(2, ax{:});
        prodx = ax' * ax;
        prodx = prodx - diag(diag(prodx));
        scorex = sum(abs(prodx(:))) / (rx * 2);
        
        ay = arrayfun(@(x) kron(betay{2}(:, x), betay{1}(:, x)), 1:ry, 'UniformOutput', false);
        ay = cellfun(@(x) x(:) / norm(x(:), 2), ay, 'UniformOutput', false);
        ay = cat(2, ay{:});
        prody = ay' * ay;
        prody = prody - diag(diag(prody));
        scorey = sum(abs(prody(:))) / (ry * 2);
    end
    
    betax = double(betax);
    betay = double(betay);
    betax = betax(:);
    betay = betay(:);

    % Save results
    save(strcat(data_path, "/results.mat"), "betax", "betay", "rho", "scorex", "scorey");

    exit;

catch ME
    disp(ME.identifier);
    disp(ME);
    disp("Matlab exited due to an error");
    exit;
end

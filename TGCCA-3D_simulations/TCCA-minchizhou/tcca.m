function [betax, betay, rho] = tcca(X, Y, rx, ry, varargin)
% TCCA Tensor canonical correlation analysis
%
% INPUT
%   X: p1-by-p2-by-...-by-px-by-n tensor data
%   Y: q1-by-q2-by-...-by-qy-by-n tensor data
%   rx: rank of X signal
%   ry: rank of Y signal
%
% OPTIONAL INPUT
%   tau: ridge shrinkage parameter. Default is 1e-3.
%   betax0: p1-by-p2-by-...-by-px tensor as an initial canonical tensor for
%          tensor data X
%   betay0: q1-by-q2-by-...-by-qy tensor as an initial canonical tensor for
%          tensor data Y
%   covstr: option for separable covariance TCCA model, can have 'none' for
%          TCCA or 'sep' for sep TCCA
%   xl0maxprop: proportion of nonzero element selected for betax
%   yl0maxprop: proportion of nonzero element selected for betay
%
% OUTPUT
%   betax: rank rx (p1,...,px) tensor signal
%   betay: rank ry (p1,...,py) tensor signal
%
% Reference: Eun Jeong Min, Eric Chi, and Hua Zhou. Tensor Canonical
% Correlation Analysis. Stat. 2020
%
% Maintainer: Eun Jeong Min <ej.min@catholic.ac.kr>


%% parse inputs
argin = inputParser;
argin.addRequired('X', @(x) isa(x,'tensor') || isnumeric(x));
argin.addRequired('Y', @(x) isa(x,'tensor') || isnumeric(x));
argin.addRequired('rx', @(x) x>=1);
argin.addRequired('ry', @(x) x>=1);
argin.addParamValue('tau', 1e-3, @(x) isnumeric(x) && x >= 0 && x <= 1);
argin.addParamValue('betax0', [], @isnumeric);
argin.addParamValue('betay0', [], @isnumeric);
argin.addParamValue('covstr', 'none', @(x) strcmp(x,'none')||strcmp(x,'sep'));
argin.addParamValue('display', 'off', @(x) strcmp(x,'off')||strcmp(x,'iter')...
  ||strcmp(x,'rep'));
argin.addParamValue('maxiter', 100, @(x) isnumeric(x) && x>0);
argin.addParamValue('replicates', 5, @(x) isnumeric(x) && x>0);
argin.addParamValue('tolfun', 1e-8, @(x) isnumeric(x) && x>0);
argin.addParamValue('xl0maxprop', [], @(x) isnumeric(x) && x>0);
argin.addParamValue('yl0maxprop', [], @(x) isnumeric(x) && x>0);
argin.parse(X, Y, rx, ry, varargin{:});

tau = argin.Results.tau;
betax0 = argin.Results.betax0;
betay0 = argin.Results.betay0;
covstr = argin.Results.covstr;
display = argin.Results.display;
maxiter = argin.Results.maxiter;
tolfun = argin.Results.tolfun;
xl0maxprop = argin.Results.xl0maxprop;
yl0maxprop = argin.Results.yl0maxprop;
if isempty(betax0)
  replicates = argin.Results.replicates;
else
  replicates = 1; % no replicates if supplied starting point
end

% convert data (X, Y) into tensors
X = tensor(X);
Y = tensor(Y);

% retrieve dimensions
dx = ndims(X) - 1;
px = size(X);
n = px(end);
dy = ndims(Y) - 1;
py = size(Y);

% cardinality constraint
if isempty(xl0maxprop)
  xl0maxprop = max(px) * rx;
end
if isempty(yl0maxprop)
  yl0maxprop = max(py) * ry;
end

% center data
Xc = double(tenmat(X, dx+1, 't')); % n columns
Xc = tensor(bsxfun(@minus, Xc, mean(Xc, 2)), size(X));
Yc = double(tenmat(Y, dy+1, 't')); % n columns
Yc = tensor(bsxfun(@minus, Yc, mean(Yc, 2)), size(Y));
clear X Y;

% pre-compute sample covariances
varxy = cell(dx, dy);
varxx = cell(dx, 1);
varyy = cell(dy, 1);
for i = 1:dx
  Xi = double(tenmat(Xc, [i 1:i-1 i+1:dx], dx+1));
  % estimate Cov(Xi)
  if strcmpi(covstr,'none')
    varxx{i} = (Xi * Xi') / n;
  elseif strcmpi(covstr,'sep')
    varxx{i} = zeros(px(i), px(i));
    for k = 1:n
      Xislice = reshape(Xi(:, k), px(i), prod(px(1:end-1)) / px(i));
      varxx{i} = varxx{i} + Xislice * Xislice';
    end
    varxx{i} = varxx{i} / n;
  end
  for j = 1:dy
    Yj = double(tenmat(Yc, [j 1:j-1 j+1:dy], dy+1));
    if i == 1
      % estimate Cov(Yj)
      if strcmpi(covstr,'none')
        varyy{j} = (Yj * Yj') / n;
      elseif strcmpi(covstr,'sep')
        varyy{j} = zeros(py(j), py(j));
        for k = 1:n
          Yjslice = reshape(Yj(:, k), py(j), prod(py(1:end-1)) / py(j));
          varyy{j} = varyy{j} + Yjslice * Yjslice';
        end
        varyy{j} = varyy{j} / n;
      end
    end
    % estimate Cov(Xi, Yj)
    varxy{i, j} = (Xi * Yj') / n;
  end
end
% scaling factor matrices
if strcmpi(covstr,'sep')
  normx = norm(Xi, 'fro')^2 / n;
  for i = 1:dx
    varxx{i} = varxx{i} / normx^(1 - 1/dx);
  end
  normy = norm(Yj, 'fro')^2 / n;
  for j = 1:dy
    varyy{j} = varyy{j} / normy^(1 - 1/dy);
  end
end
clear Xi Yj;

% main loop
opts.issym = 1;
rho_best = -inf;
for rep = 1:replicates

  % starting point
  if isempty(betax0)
    betax = ktensor(arrayfun(@(j) randn(px(j),rx), 1:dx, ...
      'UniformOutput',false));
  else
    betax = betax0;
  end
  if isempty(betay0)
    betay = ktensor(arrayfun(@(j) rand(py(j),ry), 1:dy, ...
      'UniformOutput',false));
  else
    betay = betay0;
  end
  % cumulative Hadamard products
  if strcmpi(covstr,'sep')
    cumhdx = ones(rx, rx);
    for i = 1:dx
      cumhdx = cumhdx .* (betax{i}' * varxx{i} * betax{i});
    end
    cumhdy = ones(ry, ry);
    for j = 1:dy
      cumhdy = cumhdy .* (betay{j}' * varyy{j} * betay{j});
    end
  end

  objval = -inf;
  for iter = 1:maxiter

    % mode of X
    i = mod(iter-1, dx) + 1;
    % mode of Y
    j = mod(iter-1, dy) + 1;

    % initialize cumulative Katri-Rao product
    if i == 1
      cumkrx = ones(1, rx);
    end
    if j == 1
      cumkry = ones(1, ry);
    end

    % compute the covariance matrices in subproblem
    if i == dx
      matlt = kron(cumkrx', eye(px(i)));
    else
      matlt = kron(khatrirao([betax.U(dx:-1:i+1); cumkrx])', eye(px(i)));
    end
    if j == dy
      matrt = kron(cumkry, eye(py(j)));
    else
      matrt = kron(khatrirao([betay.U(dy:-1:j+1); cumkry]), eye(py(j)));
    end
    Cxy = matlt * varxy{i, j} * matrt;
    if strcmpi(covstr,'none')
      Cxx = matlt * varxx{i} * matlt';
      Cyy = matrt' * varyy{j} * matrt;

      % add ridge to Cxx and Cyy
      Cxx = Cxx + tau * speye(size(Cxx, 1));
      Cyy = Cyy + tau * speye(size(Cyy, 1));

    elseif strcmpi(covstr,'sep')
      cumhdx = cumhdx ./ (betax{i}' * varxx{i} * betax{i});
      Cxx = kron(cumhdx, varxx{i});
      cumhdy = cumhdy ./ (betay{j}' * varyy{j} * betay{j});
      Cyy = kron(cumhdy, varyy{j});

      % add ridge to Cxx and Cyy
      Cxx = Cxx + tau ^ (1 / dx) * speye(size(Cxx, 1));
      Cyy = Cyy + tau ^ (1 / dy) * speye(size(Cyy, 1));
    end
    % add ridge to Cxx and Cyy
    % Cxx = Cxx + 1e-6 * speye(size(Cxx, 1));
    % Cyy = Cyy + 1e-6 * speye(size(Cyy, 1));
    % Cxx = Cxx + tau * speye(size(Cxx, 1));
    % Cyy = Cyy + tau * speye(size(Cyy, 1));

    if sum(sum(isnan(Cyy))) > 0
        break;
    end
    % solve the subproblem by generalized eigenvalue problem
    A = sparse([zeros(size(Cxx)) Cxy; Cxy' zeros(size(Cyy))]);
    B = sparse([Cxx zeros(size(Cxy)); zeros(fliplr(size(Cxy))), Cyy]);
    objvalold = objval;
    [evec, objval] = eigs(A, B, 1, 'lm', opts);
    betax{i} = sign(objval) * reshape(evec(1:(px(i)*rx)), px(i), rx);
    betay{j} = reshape(evec((px(i)*rx+1):end), py(j), ry);
    objval = abs(objval); % resolve sign indeterminancy

    % enforce cardinality constraint
    xisort = sort(abs(vec(betax{i})), 'descend');
    betax{i}(abs(betax{i}) < xisort(min(floor(xl0maxprop*px(i)), px(i) * rx))) = 0;
    yjsort = sort(abs(vec(betay{j})), 'descend');
    betay{j}(abs(betay{j}) < yjsort(min(floor(yl0maxprop*py(j)), py(j) * ry))) = 0;

    % accumulate Katri-Rao products
    cumkrx = khatrirao(betax{i}, cumkrx);
    cumkry = khatrirao(betay{j}, cumkry);

    % accumulate Hadamard products for separable covariance case
    if strcmpi(covstr,'sep')
      cumhdx = cumhdx .* (betax{i}' * varxx{i} * betax{i});
      cumhdy = cumhdy .* (betay{j}' * varyy{j} * betay{j});
    end

    % check stopping criterion
    if strcmpi(display,'iter')
      disp(objval);
    end
    if abs(objval - objvalold) < tolfun * (abs(objvalold) + 1) && ...
        iter > max(dx, dy)
      break;
    end

  end

  % record if we have a better correlation
  canvarx = double(ttt(Xc, tensor(betax), 1:dx));
  canvary = double(ttt(Yc, tensor(betay), 1:dy));
  rho_rep = corr(canvarx, canvary);
  if strcmpi(display, 'rep')
    disp(['replicate ' num2str(rep)]);
    disp(rho_rep);
  end
  if rho_rep > rho_best
    rho_best = rho_rep;
    betax_best = betax;
    betay_best = betay;
  end

end
rho = rho_best;
betax = betax_best;
betay = betay_best;

%%%triBeyn test
rng(1);  

ex.A = cell(3,1);
ex.mA = cell(3,1);
ex.B = cell(3,1);
ex.C = cell(3,1);
ex.D = cell(3,1);

n = 2;

for i = 1:3
    ex.A{i} = randn(n,n);
    ex.mA{i} = -ex.A{i};
    ex.B{i} = 2*randn(n,n);
    ex.C{i} = 2*randn(n,n);
    ex.D{i} = 2*randn(n,n);
end

M = [ex.A, ex.B, ex.C, ex.D];
M_multipareig = [ex.mA, ex.B, ex.C, ex.D];

% Solve with multipareig
[lambda, ~, ~, ~] = multipareig(M_multipareig);
disp('Eigenvalues (each row is [lambda1, lambda2, lambda3]):');
disp(lambda);

% Solve with triBeyn
[A0,A1] = triBeyn(M);
e3 = eig(A1,A0);

e3 = e3(abs(e3) < 1);

% plug in and solve 2D subproblems
n3 = numel(e3);
e2 = zeros(n3,1);

for k = 1:n3
    lam3 = e3(k);
    
    Mk = cell(3,3);  % (const, λ1, λ2)
    
    for i = 1:3
        A = M{i,1};
        B = M{i,2};
        C = M{i,3};
        D = M{i,4};
        
        Mk{i,1} = A + lam3 * D;  % new constant term
        Mk{i,2} = B;             % λ1 coefficient
        Mk{i,3} = C;             % λ2 coefficient
    end

    [A0,A1] = biBeyn(Mk(1:2,:));
    tmp = eig(A1,A0);
    [A0,A1] = biBeyn(Mk(2:3,:));
    tmp2 = eig(A1,A0);
    [~,ind] = min(tmp-tmp2.',[],'all');
    mat = tmp+tmp2.';
    e2(k) = mat(ind)/2;
end

% plug in and solve 1D subproblems
e1 = zeros(n3,1);

for k = 1:n3
    lam3 = e3(k);
    lam2 = e2(k);
    
    % Build three 1-parameter problems (in λ1)
    % Each is: (A_i + lam2*C_i + lam3*D_i) + λ1 * B_i
    
    eigs1 = cell(3,1);
    
    for i = 1:3
        A = M{i,1};
        B = M{i,2};
        C = M{i,3};
        D = M{i,4};
        
        A0 = A + lam2*C + lam3*D;  % constant term
        A1 = B;                    % λ1 coefficient
        
        % Solve 1-parameter problem
        [A0,A1] = beyn({A0, A1});
        eigs1{i} = eig(A1,A0);
    end
    
    % Now match eigenvalues across the three problems
    % Start with pairwise matching (1 vs 2), then include 3
    
    E12 = abs(eigs1{1} - eigs1{2}.');   % pairwise distances
    [~, idx12] = min(E12(:));
    [i1, i2] = ind2sub(size(E12), idx12);
    
    candidate = (eigs1{1}(i1) + eigs1{2}(i2)) / 2;
    
    % Match with third set
    [~, i3] = min(abs(candidate - eigs1{3}));
    
    % Final averaged value
    e1(k) = (candidate + eigs1{3}(i3)) / 2;
end


eigs = [e1 e2 e3];
eigs = eigs(abs(e1) < 1 & abs(e2) < 1 & abs(e3) < 1,:);

tol = 5e-2;

% Filter true solutions inside unit polydisc
mask = all(abs(lambda) <= 1 + tol, 2);
lambda_in = lambda(mask,:);

% Sort
lambda_sorted = sortrows(lambda_in);
eigs_sorted   = sortrows(eigs);

% One-sided matching
n_true = size(lambda_sorted,1);
matched = false(n_true,1);

for i = 1:n_true
    lam_true = lambda_sorted(i,:);
    dists = vecnorm(eigs_sorted - lam_true, 2, 2);    
    if any(dists < tol)
        matched(i) = true;
    end
end

% Report
n_matched = sum(matched);

fprintf('Matched %d / %d true solutions\n', n_matched, n_true);

if n_matched == n_true
    disp('All true solutions were found (extras allowed).');
else
    missing = find(~matched);
    fprintf('Missing %d solutions.\n', numel(missing));
end

% Residual Check
res_tol = 5e-2;
valid_eigs = [];

for k = 1:size(eigs, 1)
    lam = eigs(k, :);
    svs = zeros(3,1);
    for i = 1:3
        W_i = M{i,1} + lam(1)*M{i,2} + lam(2)*M{i,3} + lam(3)*M{i,4};
        svs(i) = min(svd(W_i));
    end
    
    if max(svs) < res_tol
        valid_eigs = [valid_eigs; lam];
    end
end

fprintf('\nPassed residual check: %d / %d solutions\n', size(valid_eigs,1), size(eigs,1));

% 3D Visualization
figure('Color', 'w', 'Name', '3D Solution Map');
hold on; grid on;

% ================= SETTINGS & OPTIONS ================= %
cubeColor   = [0.2 0.5 0.8]; % Blue
markerColor = [0.85 0.3 0.2];% Red
cubeAlpha   = 0.1;           % Transparency
lineAlpha   = 0.8;           % Drop line visibility
limit       = -1;            % Plane boundaries for drop lines
% ====================================================== %

% Draw Transparent Blue Unit Cube
v = [-1 -1 -1; 1 -1 -1; 1 1 -1; -1 1 -1; -1 -1 1; 1 -1 1; 1 1 1; -1 1 1];
f = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
patch('Vertices', v, 'Faces', f, 'FaceColor', cubeColor, ...
      'FaceAlpha', cubeAlpha, 'EdgeColor', cubeColor*0.7);

if ~isempty(valid_eigs)
    X = real(valid_eigs(:,1));
    Y = real(valid_eigs(:,2));
    Z = real(valid_eigs(:,3));
    
    % Drop Lines
    for i = 1:length(X)
        % Line to Z-plane
        line([X(i) X(i)], [Y(i) Y(i)], [limit Z(i)], 'Color', [0.4 0.4 .4 lineAlpha], 'LineStyle', ':','LineWidth', 5);
    end
    
    % Solution Markers
    scatter3(X, Y, Z, 600, markerColor, 'filled', 'MarkerEdgeColor', 'k');
end

% Labels
xlabel('Re(\lambda_1)'); ylabel('Re(\lambda_2)'); zlabel('Re(\lambda_3)');
xticks([-1 -.5 0 .5 1]);
yticks([-1 -.5 0 .5 1]);
zticks([-1 -.5 0 .5 1]);
title('3D Solution Projections');
view(-35, 30);
axis equal;
set(gca, 'Projection', 'orthographic');
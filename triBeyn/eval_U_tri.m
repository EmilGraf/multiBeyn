function Ut = eval_U_tri(t, z, Qcoeff)
    % eval_U Evaluates the component U_t at point z in C^3.
    %
    % INPUTS:
    %   t        : Integer index (1, 2, or 3)
    %   z        : 3x1 complex vector representing the point z
    %   Q,Qh,dQH
   
    [Q, QH, dQH] = eval_Q_all(z, Qcoeff);

    % Determine the derivative indices r and s based on input t
    if t == 1
        r1 = 2; s1 = 3;  r2 = 3; s2 = 2;
    elseif t == 2
        r1 = 1; s1 = 3;  r2 = 3; s2 = 1;
    elseif t == 3
        r1 = 1; s1 = 2;  r2 = 2; s2 = 1;
    else
        error('Input t must be 1, 2, or 3.');
    end

    % Evaluate the R terms
    R_r1_s1 = compute_R(r1, s1, Q, QH, dQH);
    R_r2_s2 = compute_R(r2, s2, Q, QH, dQH);
    
    Ut = R_r1_s1 - R_r2_s2;
end

function R = compute_R(r, s, Q, QH, dQH)
    % Evaluates R_{rs} assuming i=1, j=2, k=3
    
    % Term 1A: lead index = 1, set = {1,2,3} -> I = {2,3}
    term1A = compute_drhoH(1, [2, 3], r, Q, QH, dQH);
    
    % Term 1B: lead index = 2, set = {2,3} -> I = {3}
    term1B = compute_drhoH(2, 3, s, Q, QH, dQH);
    
    % Term 2A: lead index = 2, set = {1,2,3} -> I = {1,3}
    term2A = compute_drhoH(2, [1, 3], r, Q, QH, dQH);
    
    % Term 2B: lead index = 1, set = {1,3} -> I = {3}
    term2B = compute_drhoH(1, 3, s, Q, QH, dQH);

    R = term1A * term1B - term2A * term2B;
end

function drhoH = compute_drhoH(idx, I, ell, Q, QH, dQH)
    % compute_drhoH Calculates \overline{\partial}\tilde{\rho}_{idx, I \cup \{idx\}, \ell}^H
    
    % Compute D_{(I)}
    D_I = zeros(size(Q{1}));
    for m = 1:3
        if ismember(m, I)
            D_I = D_I + QH{m} * Q{m};
        else
            D_I = D_I + Q{m} * QH{m};
        end
    end
    DI_inv = inv(D_I);

    % Compute S_{(I), \ell}
    S_I_ell = zeros(size(Q{1}));
    for m = 1:3
        if ismember(m, I)
            S_I_ell = S_I_ell + dQH{m, ell} * Q{m};
        else
            S_I_ell = S_I_ell + Q{m} * dQH{m, ell};
        end
    end

    % Evaluate based on the size of I
    if length(I) > 1
        % First term in the wedge product (e.g., I = {j, k})
        drhoH = (-1)^(idx-1) * (dQH{idx, ell} - QH{idx} * DI_inv * S_I_ell) * DI_inv;

    elseif length(I) == 1
        % Second term in the wedge product (e.g., I = {k})
        Q_k_inv = inv(Q{I(1)});
        drhoH = (-1)^(idx-1) * (dQH{idx, ell} - QH{idx} * DI_inv * S_I_ell) * DI_inv * Q_k_inv;
        
    else
        error('Subset I cannot be empty.');
    end
end

function [Qz, QHz, dQHz] = eval_Q_all(z, Qcoeff)
    % Outputs:
    % Qz   : 1x3 cell, Q_m(z)
    % QHz  : 1x3 cell, Q_m(z)^H
    % dQHz : 3x3 cell, (dQ_m/dz_ell)^H
    
    Qz   = cell(1,3);
    QHz  = cell(1,3);
    dQHz = cell(3,3);
    
    for m = 1:3
        % Evaluate Q_m(z)
        Qz{m} = Qcoeff{m,1} + ...
                z(1)*Qcoeff{m,2} + ...
                z(2)*Qcoeff{m,3} + ...
                z(3)*Qcoeff{m,4};
        
        % Hermitian
        QHz{m} = Qz{m}';
        
        % Derivatives (constant in z)
        for ell = 1:3
            dQHz{m,ell} = Qcoeff{m,ell+1}';
        end
    end
end


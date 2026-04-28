function Ut = eval_U_bi(t, z, Qcoeff)
    % eval_U Evaluates the component U_t at point z in C^2.
    
    % Retrieve matrices evaluated at z
    [Q, QH, dQH] = eval_Q_all(z, Qcoeff); 
    
    % Fixed indices for the 2D formulation
    idx = 1; % Corresponds to i = 1
    I = 2;   % Corresponds to subset {j} where j = 2

    % Map component index t to derivative index ell
    if t == 1
        ell = 2;
    elseif t == 2
        ell = 1;
    else
        error('Input t must be 1 or 2.');
    end
    
    % Evaluate the 1-form directly
    Ut = compute_drhoH(idx, I, ell, Q, QH, dQH);
end

function drhoH = compute_drhoH(idx, I, ell, Q, QH, dQH)
    % compute_drhoH Calculates \overline{\partial}\tilde{\rho}_{idx, I \cup \{idx\}, \ell}^H
    
    % Compute D_{(I)}
    D_I = zeros(size(Q{1}));
    for m = 1:2
        if ismember(m, I)
            D_I = D_I + QH{m} * Q{m};
        else
            D_I = D_I + Q{m} * QH{m};
        end
    end
    DI_inv = inv(D_I);

    % Compute S_{(I), \ell} 
    S_I_ell = zeros(size(Q{1}));
    for m = 1:2
        if ismember(m, I)
            S_I_ell = S_I_ell + dQH{m, ell} * Q{m};
        else
            S_I_ell = S_I_ell + Q{m} * dQH{m, ell};
        end
    end

    % Evaluate the differential form
    % Since I is always {2}, |I| = 1.
    Q_I_inv = inv(Q{I});
    
    drhoH = (-1)^(idx-1) * (dQH{idx, ell} - QH{idx} * DI_inv * S_I_ell) * DI_inv * Q_I_inv;
end

function [Qz, QHz, dQHz] = eval_Q_all(z, Qcoeff)
    % Outputs:
    % Qz   : 1x2 cell, Q_m(z)
    % QHz  : 1x2 cell, Q_m(z)^H
    % dQHz : 2x2 cell, (dQ_m/dz_ell)^H
    
    Qz   = cell(1,2);
    QHz  = cell(1,2);
    dQHz = cell(2,2);
    
    for m = 1:2
        % Evaluate Q_m(z)
        Qz{m} = Qcoeff{m,1} + ...
                z(1)*Qcoeff{m,2} + ...
                z(2)*Qcoeff{m,3};
        
        % Hermitian
        QHz{m} = Qz{m}';
        
        % Derivatives (constant in z)
        for ell = 1:2
            dQHz{m,ell} = Qcoeff{m,ell+1}';
        end
    end
end

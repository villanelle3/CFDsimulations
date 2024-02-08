function [source_term_u, source_term_v] = wall_source_term_spectral(u_hat, v_hat, y, Re)
    [~, ~, ~, nu, ~, ~, ~, ~, ~, ~, ~, kappa, ~, ~] = constantes(Re);
    % kappa = Von Karman constant
    yplus_crit = 11.06; % Critical wall distance in viscous sublayer
    % Calculate wall distance in wall units (y+)
    u_tau = sqrt(nu * abs(u_hat(2))/abs(y(2))); % Friction velocity
    yplus = y * u_tau / nu;

    % Compute source term based on the logarithmic law for u_hat
    source_term_u = zeros(size(u_hat));
    for i = 1:length(y)
        if yplus(i) <= yplus_crit
            % In the viscous sublayer
            source_term_u(:, i) = 0; % No source term needed
        else
            % In the log-law layer
            source_term_u(:, i) = kappa / log(yplus(i)); % Logarithmic law source term for u_hat
        end
    end

    % Compute source term based on the logarithmic law for v_hat
    source_term_v = zeros(size(v_hat));
    for i = 1:length(y)
        if yplus(i) <= yplus_crit
            % In the viscous sublayer
            source_term_v(:, i) = 0; % No source term needed
        else
            % In the log-law layer
            source_term_v(:, i) = kappa / log(yplus(i)); % Logarithmic law source term for v_hat
        end
    end

    % Repeat the source terms for each Fourier mode
    % source_term_u = repmat(source_term_u, size(u_hat, 1), 1);
    % source_term_v = repmat(source_term_v, size(v_hat, 1), 1);
end

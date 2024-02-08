function [viscx, viscy] = calculate_viscous_term(u_hat, v_hat, k_hat, eps_hat, kx, ky, k2, Re)
    [~, ~, ~, nu, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = constantes(Re);
    % Calcula as derivadas espectrais
    du_dx_hat = 1i * kx' .* u_hat;
    dv_dy_hat = 1i * ky .* v_hat;

    % Calcula a viscosidade turbulenta
    mut_hat = calculate_turbulent_viscosity(k_hat, eps_hat, nu);

    % Calcula os termos viscosos espectrais
    viscx = mut_hat .* (-k2 .* u_hat - du_dx_hat);
    viscy = mut_hat .* (-k2 .* v_hat - dv_dy_hat);
end

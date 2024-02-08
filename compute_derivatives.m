function [rhs_u, rhs_k] = compute_derivatives(u_hat, v_hat, k_hat, eps_hat, mut, kx, ky, k2, nx, ny, Re)

    [viscx, viscy] = visc(Re, nx, ny, kx, ky, mut, u_hat, v_hat, k_hat, eps_hat);
    % [~, ~, ~, mudynamic, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = constantes(Re);
    % nief = calculate_turbulent_viscosity(k_hat, eps_hat, mudynamic);
    % viscx=nief.*k2.*u_hat;
	% viscy=nief.*k2.*v_hat;


    [tnlxt, tnlyt] = tnlinear(kx, ky, nx, ny, u_hat, v_hat);

    % rhs_u = - tnlxt + viscx;
    % rhs_k = - tnlyt + viscy;
    rhs_u = -tnlxt - viscx;
    rhs_k = -tnlyt - viscy;

    [rhs_u, rhs_k] = projecao(rhs_u, rhs_k, nx, ny, kx, ky, k2);
end
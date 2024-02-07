function [rhs_u, rhs_k] = compute_derivatives(u_hat, v_hat, mut, kx, ky, k2, nx, ny, Re)

    [viscx, viscy] = visc(Re, nx, ny, kx, ky, mut, u_hat, v_hat);
    [tnlxt, tnlyt] = tnlinear(kx, ky, nx, ny, u_hat, v_hat);

    s = zeros(nx, ny);
    % rhs_u = - tnlxt + viscx + s;
    % rhs_k = - tnlyt + viscy;
    rhs_u = - tnlxt - viscx + s;
    rhs_k = - tnlyt - viscy;

    [rhs_u, rhs_k] = projecao(rhs_u, rhs_k, nx, ny, kx, ky, k2);
end
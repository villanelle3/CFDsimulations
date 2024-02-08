function [viscx, viscy] = visc(Re, nx, ny, kx, ky, mut, ut, vt, k_hat, eps_hat)
    [~, ~, ~, mudynamic, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = constantes(Re);
    dudx = zeros(nx, ny); 
    dudy = zeros(nx, ny); 
    dvdx = zeros(nx, ny); 
    dvdy = zeros(nx, ny);
    for i = 1:nx
        for j = 1:ny
            dudx(i,j) = 1i * kx(i) * ut(i,j);
            dudy(i,j) = 1i * ky(j) * ut(i,j);
            dvdx(i,j) = 1i * kx(i) * vt(i,j);
            dvdy(i,j) = 1i * ky(j) * vt(i,j);
        end
    end
    dudx = real(ifft2(dudx));
    dudy = real(ifft2(dudy));
    dvdx = real(ifft2(dvdx));
    dvdy = real(ifft2(dvdy));

    % Compute effective viscosity
    % nief = mut + mudynamic;
    nief = calculate_turbulent_viscosity(k_hat, eps_hat, mudynamic);

    
    % Compute stress tensor components in Fourier space
    tauxx = 2.0 * nief .* dudx;
    tauyy = 2.0 * nief .* dvdy;
    tauxy = nief .* (dudy + dvdx);
    
    % Apply FFT to stress tensor components
    tauxx = fft2(tauxx);
    tauyy = fft2(tauyy);
    tauxy = fft2(tauxy);
    
    % Initialize viscous terms
    viscx = zeros(size(dudx));
    viscy = zeros(size(dudy));
    
    % Compute viscous terms in physical space
    for i = 1:size(dudx, 1)
        for j = 1:size(dudx, 2)
            viscx(i, j) = 1i * (kx(i) * tauxx(i, j) + ky(j) * tauxy(i, j));
            viscy(i, j) = 1i * (kx(i) * tauxy(i, j) + ky(j) * tauyy(i, j));
        end
    end
    
end

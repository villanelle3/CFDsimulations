function [kx, ky, k2] = calcular_wave_numbers2(Nx, Ny, Lx, Ly)

    kx = (2*pi/Lx) * [0:Nx/2-1, 0, -Nx/2+1:-1]; % Wave numbers in x-direction
    ky = (2*pi/Ly) * [0:Ny/2-1, 0, -Ny/2+1:-1]; % Wave numbers in y-direction

    
    % Calculando k^2 (o quadrado do número de onda)
    [kx_grid, ky_grid] = meshgrid(kx, ky);
    k2 = kx_grid.^2 + ky_grid.^2 + eps/2;
end

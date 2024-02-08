function [kx, ky, k2] = define_wave_numbers2(nx, ny, Lx, Ly)
    % Define streamwise wave numbers (kx)
    kx = (0:(nx/2-1)) * (2 * pi / Lx); % Normalized to span from 0 to 2*pi/Lx
    kx = [kx, (-nx/2:-1) * (2 * pi / Lx)];

    % Define wall-normal wave numbers (ky)
    ky = (0:(ny/2-1)) * (pi / Ly); % Normalized to span from 0 to pi/Ly
    ky = [ky, (-ny/2:-1) * (pi / Ly)];

    % Compute the square of the wave numbers (k2)
    [KX, KY] = meshgrid(kx, ky);
    k2 = KX.^2 + KY.^2 + eps/2;
end

function [kx, ky, k2] = define_wave_numbers(nx, ny)
    % Define streamwise wave numbers (kx)
    kx = (0:nx-1)' * 2 * pi / nx; % Normalized to span from 0 to 2*pi

    % Define wall-normal wave numbers (ky)
    ky = (0:ny/2) * pi / ny; % Normalized to span from 0 to pi

    % Append negative wall-normal wave numbers
    if mod(ny, 2) == 0 % If ny is even
        ky = [ky, -fliplr(ky(2:end-1))]; % Exclude ky = 0 and duplicate negative values
    else
        ky = [ky, -fliplr(ky(2:end))]; % Exclude ky = 0 and duplicate negative values
    end

    % Compute the square of the wave numbers
    [KX, KY] = meshgrid(kx, ky);
    k2 = KX.^2 + KY.^2;
end

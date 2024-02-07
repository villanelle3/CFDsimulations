% function [ax, ay] = projecao(ax, ay, kx, ky, k2)
%     % Inputs:
%     % ax, ay: complex arrays representing velocity components
%     % kx, ky: wavenumbers
%     % k2: square of wavenumbers
% 
%     % Compute projection
%     tmp = (ax .* kx + ay .* ky) ./ k2;
%     ax = ax - tmp .* kx;
%     ay = ay - tmp .* ky;
% end

function [ax, ay] = projecao(ax, ay, nx, ny, kx, ky, k2)
    tmp = zeros(nx, ny);
    for i = 1:nx
        for j = 1:ny
            tmp(i, j) = (ax(i, j) * kx(i) + ay(i, j) * ky(j)) / k2(i, j);
            ax(i, j) = ax(i, j) - tmp(i, j) * kx(i);
            ay(i, j) = ay(i, j) - tmp(i, j) * ky(j);
        end
    end
end

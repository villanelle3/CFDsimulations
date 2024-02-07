function [kx, ky, k2] = numero_de_onda(Nx, Ny)
    kx = [0:Nx/2-1 0 -Nx/2+1:-1];        
    ky = [0:Ny/2-1 0 -Ny/2+1:-1]';     
    [k2xm, k2ym] = meshgrid(kx.^2, ky.^2);
    k2 = k2xm + k2ym + eps/2;
    k2 = k2';
end
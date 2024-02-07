function [kx, ky, k2] = calcular_wave_numbers(nx, ny, dx, dy)
    % Pré-alocação dos vetores de número de onda na direção "x" e "y"
    kx = zeros(1, nx);
    ky = zeros(1, ny);
    k2 = zeros(nx, ny);

    % Definição dos vetores de número de onda na direção "x"
    for i = 1:nx/2+1
        kx(i) = 2.0*pi*(i-1)/(nx*dx);
    end

    for i = nx/2+2:nx
        kx(i) = 2.0*pi*(i-1-nx)/(nx*dx);
    end

    % Definição dos vetores de número de onda na direção "y"
    for j = 1:ny/2+1
        ky(j) = 2.0*pi*(j-1)/(ny*dy);
    end

    for j = ny/2+2:ny
        ky(j) = 2.0*pi*(j-1-ny)/(ny*dy);
    end

    % Constante pequena para estabilidade numérica
    small = 1.0e-30;

    % Cálculo do quadrado do número de onda k2
    for j = 1:ny
        for i = 1:nx
            k2(i,j) = kx(i)^2 + ky(j)^2 + small;
        end
    end
end

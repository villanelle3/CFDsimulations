s = WallConditions.ComputeSourceTerm(u, k);


% ------------------------------------------------------------------------
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
            k2 = KX.^2 + KY.^2 + eps/2;
        end
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
        function [kx, ky, k2] = numero_de_onda(Nx, Ny)
            kx = [0:Nx/2-1 0 -Nx/2+1:-1];        
            ky = [0:Ny/2-1 0 -Ny/2+1:-1]';     
            [k2xm, k2ym] = meshgrid(kx.^2, ky.^2);
            k2 = k2xm + k2ym + eps/2;
            k2 = k2';
        end
        % ------------------------------------------------------------------------
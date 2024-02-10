classdef WallConditions
    methods (Static)
        function s = ComputeSourceTerm(u, k)
            % Calcula o termo fonte s nas paredes superior e inferior
            
            % Parâmetros
            rho = Canal.rho;
            cmu = Canal.c_mu;
            kappa = Canal.kappa;
            B = Canal.B;
            dx = Canal.dx;
            dy = Canal.dy;
            nx = Canal.nx;
            ny = Canal.ny;
            E = exp(kappa * B);
            
            % Inicializar o termo fonte s
            s = zeros(nx, ny);
            
            % Calcular para a parede superior (ny)
            for i = 2:nx-1
                ystar = rho * cmu^0.25 * sqrt(k(i, ny-1)) * dy / Canal.nu;
                tauwall = rho * kappa * cmu^0.25 * k(i, ny-1)^0.5 * u(i, ny-1)/log(E * ystar);
                s(i, ny-1) = rho * cmu^0.25 * sqrt(k(i, ny-1)) * dx / (log(E * ystar) / kappa);
            end
            
            % Calcular para a parede inferior (1)
            for i = 2:nx-1
                ystar = rho * cmu^0.25 * sqrt(k(i, 2)) * dy / Canal.nu;
                tauwall = rho * kappa * cmu^0.25 * k(i,ny-1)^0.5 * u(i, 2)/log(E*ystar);
                s(i, 2) = rho * cmu^0.25 * sqrt(k(i, 2)) * dx / (log(E * ystar) / kappa);
            end
        end

        function waall_star = CalculateWallStar(k)
            % Verificar se waall_star >= 30 para que a função de parede seja válida.
            waall_star = Canal.rho * Canal.c_mu^0.25 * sqrt(k(Canal.nx - 1, 2)) * (Canal.dy)/Canal.nu;
        end
    end
end

classdef Canal
    properties (Constant)
        % Parâmetros do canal
        Lx = 10;              % Comprimento do canal
        Ly = 1;               % Largura do canal
        nx = 8;               % Número de pontos de grade na direção do comprimento
        ny = 16;              % Número de pontos de grade na direção da parede
        Re = 40000;           % Número de Reynolds
        nu = 1/Canal.Re;      % Viscosidade cinemática
        uinlet = 1.0;

        % Constantes do modelo
        rho = 1; 
        c_mu = 0.09;
        sigma_k = 1;
        sigma_eps = 1.22;
        c_1eps = 1.44;
        c_2eps = 1.92;
        omegaeps = 0.3;
        omegak = 0.3;   
        kappa = 0.41;   
        B = 5.5;
        ti = 0.1;       % Intensidade turbulenta
    end
    
    properties (Constant, Access = public)
        x = linspace(0, Canal.Lx, Canal.nx);   
        y = linspace(0, Canal.Ly, Canal.ny);
        dx = Canal.Lx / (Canal.nx - 1);         % Espaçamento de grade na direção do comprimento (x)
        dy = Canal.Ly / (Canal.ny - 1);         % Espaçamento de grade na direção da largura (y)
    end
    
    methods (Static)
        function [u, v, k, epsilon, mu_t] = AplicarCondicaoInicial()
            % Método para aplicar as condições iniciais
            
            % Inicializar matrizes para as variáveis do escoamento
            u = Canal.uinlet * ones(Canal.nx, Canal.ny); % Componente da velocidade na direção do escoamento (constante em um canal)
            v = zeros(Canal.nx, Canal.ny); % Componente da velocidade na direção da largura (zero nas paredes)
            
            % Garantir que as primeiras e últimas linhas de u e v representem as paredes
            u(:, 1) = 0; % Velocidade na parede inferior
            u(:, end) = 0; % Velocidade na parede superior
            v(:, 1) = 0; % Velocidade na parede inferior na direção da largura
            v(:, end) = 0; % Velocidade na parede superior na direção da largura
            
            % Calcular energia cinética turbulenta (k)
            k = 2/3 * (u .* Canal.ti).^2;
            
            % Calcular taxa de dissipação turbulenta (epsilon)
            Ls = 0.07 * Canal.Ly;
            epsilon = Canal.c_mu^(3/4) * k.^(3/2) ./ Ls;
            
            % Calcular viscosidade turbulenta (mu_t)
            mu_t = Canal.rho * Canal.c_mu .* k.^2 ./ epsilon;
            
            % Substituir valores NaN e Inf em mu_t por 0
            mu_t(isnan(mu_t) | isinf(mu_t)) = 0;
        end
        
        function [kx, ky, k2] = CalcularNumerosOnda()
            % Método para calcular os vetores número de onda
                        
            % Definir os vetores de números de onda na direção do comprimento (x) e da largura (y)
            kx = 2 * pi / Canal.Lx * [0:1:(Canal.nx/2 - 1), 0, -Canal.nx/2+1:1:-1]; % Vetor de números de onda em x
            ky = 2 * pi / Canal.Ly * [0:1:(Canal.ny/2 - 1), 0, -Canal.ny/2+1:1:-1]; % Vetor de números de onda em y
          
            % Calcular o quadrado do vetor número de onda
            [KX, KY] = meshgrid(kx, ky);
            k2 = KX.^2 + KY.^2 + eps/2;
        end
    end
end

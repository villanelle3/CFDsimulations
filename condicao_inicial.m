function [x, y, u, v, k, eps, mut, rows, columns, dx, dy] = condicao_inicial(Re, nx, ny, xmax, xmin, ymax, ymin)

    [rho, ~, ~, ~, ~, ~, ~, ~, ~, cmu, ~, ~, ~, ~] = constantes(Re);
    
    Lx = xmax - xmin;           % Largura do canal
    Ly = ymax - ymin;           % Altura do canal
    dx = Lx/(nx-1);         
    dy = Ly/(ny-1);
    x = xmin:dx:xmax;       
    y = ymin:dy:ymax; 

    % --------------------------------- u ---------------------------------

    u = ones(nx, ny)';          
    rows = height(u);       
    columns = width(u); 
    u(1, :) = 0;                % Parede de cima
    u(rows, :) = 0;             % Parede de baixo

    % --------------------------------- v ---------------------------------

    v = zeros(rows, columns);

    % --------------------------------- k ---------------------------------

    I = 0.1;                    % Intensidade turbulenta
    k = 2/3*(u.*I).^2;

    % -------------------------------- eps --------------------------------
    
    Ls = 0.07*(ymax - ymin);
    eps = cmu^(3/4).*k.^(3/2)./Ls;
    
    % -------------------------------- mut --------------------------------

    mut = rho*cmu.*k.^2./eps;
    indD = isnan(mut)|isinf(mut);
    mut(indD) = 0;

end
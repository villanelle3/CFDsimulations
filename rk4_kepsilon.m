function [u_hat, v_hat] = rk4_kepsilon(u_hat, v_hat, mut, kx, ky, k2, nx, ny, Re, dt)
    % RK4 time integration

    % Step 1: Compute k1
    [u_k1, v_k1] = compute_derivatives(u_hat, v_hat, mut, kx, ky, k2, nx, ny, Re);
    % [u_k1, v_k1, k_k1, epsilon_k1] = compute_derivatives(u_hat, v_hat, k_hat, epsilon_hat);

    % Step 2: Compute k2
    [u_k2, v_k2] = compute_derivatives(u_hat + 0.5*dt*u_k1, v_hat + 0.5*dt*v_k1, mut, kx, ky, k2, nx, ny, Re);
    % [u_k2, v_k2, k_k2, epsilon_k2] = compute_derivatives(u_hat + 0.5*dt*u_k1, v_hat + 0.5*dt*v_k1, ...
    %                                                      k_hat + 0.5*dt*k_k1, epsilon_hat + 0.5*dt*epsilon_k1);

    % Step 3: Compute k3
    [u_k3, v_k3] = compute_derivatives(u_hat + 0.5*dt*u_k2, v_hat + 0.5*dt*v_k2, mut, kx, ky, k2, nx, ny, Re);
    % [u_k3, v_k3, k_k3, epsilon_k3] = compute_derivatives(u_hat + 0.5*dt*u_k2, v_hat + 0.5*dt*v_k2, ...
    %                                                      k_hat + 0.5*dt*k_k2, epsilon_hat + 0.5*dt*epsilon_k2);

    % Step 4: Compute k4
    [u_k4, v_k4] = compute_derivatives(u_hat + dt*u_k3, v_hat + dt*v_k3, mut, kx, ky, k2, nx, ny, Re);
    % [u_k4, v_k4, k_k4, epsilon_k4] = compute_derivatives(u_hat + dt*u_k3, v_hat + dt*v_k3, ...
    %                                                      k_hat + dt*k_k3, epsilon_hat + dt*epsilon_k3);

    % Update variables using RK4 formula
    u_hat = u_hat + (1/6)*(u_k1 + 2*u_k2 + 2*u_k3 + u_k4)*dt;
    v_hat = v_hat + (1/6)*(v_k1 + 2*v_k2 + 2*v_k3 + v_k4)*dt;
    % k_hat = k_hat + (1/6)*(k_k1 + 2*k_k2 + 2*k_k3 + k_k4)*dt;
    % epsilon_hat = epsilon_hat + (1/6)*(epsilon_k1 + 2*epsilon_k2 + 2*epsilon_k3 + epsilon_k4)*dt;


    % Regularize k e ep. Não faz sentido ter valores negativos.
    % k(k<0) = 0;
   	% ep(ep<0) = 0;

    % Apply boundary conditions to ensure zero velocity at the top and bottom walls
    u_hat(1,:) = 0;      % Velocidade na parede superior
    u_hat(end,:) = 0;    % Velocidade na parede inferior
    v_hat(1,:) = 0;      % Velocidade na parede superior
    v_hat(end,:) = 0;    % Velocidade na parede inferior
    
end









% -------------------------------------------------------------------------
function EPS = fix_nan(EPS)
    indD = isnan(EPS)|isinf(EPS);
    EPS(indD) = 0;
end

function var = wallBC(var)
    rows = height(var);
    var(1, :) = 0;                % Parede de cima
    var(rows, :) = 0;             % Parede de baixo
end

function [sx] = WallBC_NS(Nx, Ny, K, rho, cmu, dx, dy, E, kappa, mudynamic)

    % Calcular a tensão de cisalhamento da parede como um termo fonte
    sx = zeros(Nx, Ny);  
    
    K_upperWall = K(2, 2:end-1);      % Parede de cima (penultima linha da matriz) 
    K_lowerWall = K(Nx-1, 2:end-1);   % Parede de baixo (penultima linha) 
    
    for i = 1:length(K_upperWall)
        ystar = rho * cmu^0.25 * sqrt(K_upperWall(i)) * (dy)/mudynamic;
        sx(2, i+1) = rho*cmu^0.25*K_upperWall(i)^0.5*(dx)/(log(E*ystar)/kappa);
        ystar = rho * cmu^0.25 * sqrt(K_lowerWall(i)) * (dy)/mudynamic;
        sx(Nx-1, i+1) = rho*cmu^0.25*K_lowerWall(i)^0.5*(dx)/(log(E*ystar)/kappa);
    end

end



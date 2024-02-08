% Limpa a janela de comando e a área de trabalho
clear;
clc;

% Parâmetros do problema
Re = 40000;     % Número de Reynolds
nx = 16;        % Número de pontos na direção x
ny = 8;         % Número de pontos na direção y
xmin = 0;       % Limite inferior em x
xmax = 10;      % Limite superior em x
ymin = 0;       % Limite inferior em y
ymax = 1;       % Limite superior em y

%---------------------------- Condição inicial ----------------------------

% Gera a condição inicial do problema
[x, y, u, v, K, epsilon, mut, rows, columns, dx, dy] = ...
    condicao_inicial(Re, nx, ny, xmax, xmin, ymax, ymin);

% Transformada de Fourier das componentes de velocidade
u_hat = fft2(u);   
v_hat = fft2(v);
k_hat = fft2(K);
eps_hat = fft2(epsilon);

% Números de onda
% [kx, ky, k2] = calcular_wave_numbers(rows, columns, dx, dy);    % kx=1x8, ky=1x16, k2 = 8x16.  -0.5366,  0.3384, 0.1616 

[kx, ky, k2] = define_wave_numbers2(rows, columns, 10, 1);
k2 = k2';

% Parâmetros de integração temporal
tf = 10;            % Tempo final da simulação
dt = 10e-4;         % Passo de tempo

% -------------------------------------------------------------------------

% Loop principal de integração temporal
for i = 1:ceil(tf/dt)
    [u_hat, v_hat] = rk4_kepsilon(u_hat, v_hat, k_hat, eps_hat, mut, kx, ky, k2, rows, columns, Re, dt);
    % Atualiza o tempo
    tempo = i*dt;
    % disp(tempo)
end

u = real(ifft2(u_hat));
v = real(ifft2(v_hat));

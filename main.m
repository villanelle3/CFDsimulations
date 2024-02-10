clear; clc;

% Condição Inicial. Os dados de entrada são definidos na classe Canal
[kx, ky, k2] = Canal.CalcularNumerosOnda();
[u, v, k, epsilon, mu_t] = Canal.AplicarCondicaoInicial();

waall_star = WallConditions.CalculateWallStar(k);
s = WallConditions.ComputeSourceTerm(u, k);

% contourf(x, y', u');
% colorbar;
% xlabel('Comprimento do Canal (x)');
% ylabel('Largura do Canal (y)');
% title('Perfil de Velocidade na Direção do Escoamento (u)');




function mut = calculate_turbulent_viscosity(k_hat, eps_hat, nu)
    % Calcula a viscosidade turbulenta usando o modelo k-epsilon

    % Defina o coeficiente Cmu
    Cmu = 0.09;

    % Transforma as variáveis para o espaço físico
    k = real(ifft2(k_hat));
    eps = real(ifft2(eps_hat));

    % Evita divisão por zero
    eps(eps < 1e-12) = 1e-12;

    % Calcule a escala de comprimento turbulento
    l = (Cmu * k.^2 ./ eps).^0.5;

    % Calcule a viscosidade turbulenta
    mut = nu + nu * l.^2 .* sqrt(k) ./ eps;
end
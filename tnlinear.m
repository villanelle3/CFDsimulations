function [tnlx, tnly] = tnlinear(kx, ky, nx, ny, ut, vt)
    utfis = real(ifft2(ut));
    vtfis = real(ifft2(vt));

    uu = utfis .* utfis;
    vv = vtfis .* vtfis;
    uv = utfis .* vtfis;

    uu = fft2(uu);
    vv = fft2(vv);
    uv = fft2(uv);

    tnlxa = zeros(nx, ny); % Preallocate tnlxa matrix
    tnlya = zeros(nx, ny); % Preallocate tnlya matrix

    for i = 1:nx
        for j = 1:ny
            tnlxa(i,j) = 1i * kx(i) * uu(i,j) + 1i * ky(j) * uv(i,j);
            tnlya(i,j) = 1i * kx(i) * uv(i,j) + 1i * ky(j) * vv(i,j);
        end
    end

    dudx = zeros(nx, ny); % Preallocate dudx matrix
    dudy = zeros(nx, ny); % Preallocate dudy matrix
    dvdx = zeros(nx, ny); % Preallocate dvdx matrix
    dvdy = zeros(nx, ny); % Preallocate dvdy matrix
    
    for i = 1:nx
        for j = 1:ny
            dudx(i,j) = 1i * kx(i) * ut(i,j);
            dudy(i,j) = 1i * ky(j) * ut(i,j);
            dvdx(i,j) = 1i * kx(i) * vt(i,j);
            dvdy(i,j) = 1i * ky(j) * vt(i,j);
        end
    end

    dudx = real(ifft2(dudx));
    dudy = real(ifft2(dudy));
    dvdx = real(ifft2(dvdx));
    dvdy = real(ifft2(dvdy));

    tnlxc = utfis .* dudx + vtfis .* dudy;
    tnlyc = utfis .* dvdx + vtfis .* dvdy;

    tnlxc = fft2(tnlxc);
    tnlyc = fft2(tnlyc);

    tnlx = 0.5 * (tnlxa + tnlxc);
    tnly = 0.5 * (tnlya + tnlyc);

end

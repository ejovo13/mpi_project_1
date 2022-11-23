function [x, y, z, col] = plot_med_data()

    path_med = "../csv/ETOPO1_med.csv";
    load(path_med)
    tmp = ETOPO1_med;

    sph = zeros(size(tmp));

    sph(:,1) = phi_to_the(tmp(:,2));
    sph(:,2) = lam_to_phi(tmp(:,1));
    sph(:,3) = tmp(:,3);

    nt = sqrt(length(sph) / 2);
    np = 2 * nt;
    dt = pi / nt;
    dp = 2 * pi / np;

    r = reshape(sph(:,3), [nt np]);

    col = r;
    r = r + 100000; % arbitrary height to show differences in altitutde

    azimuth = linspace(-pi, pi - dp, np)';
    elevation = linspace(pi/2, (-pi/2) + dt, nt)';

    [th, ph] = meshgrid(azimuth, elevation);

    [x, y, z] = sph2cart(th, ph, r);

    s = surface(x, y, z, col);

    xlabel('th')
    ylabel('ph')
    % colorbar
    % colormap winter
    colormap jet

    axis equal
    shading interp

    CAZ = 15.350734615445207;
    CEL = 18.233845566149686;

    view(CAZ, CEL)
    axis off

    saveas(gcf, strcat("med.jpg"))

end

function phi = lam_to_phi(lam)

    phi = lam;

end

function theta = phi_to_the(phi)

    theta = (pi / 2) - phi;

end
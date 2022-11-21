% Let's use the medium data set and the information stored in the 
% diff files to plot degrees of different models

% plot_wide_range()
plot_small_range()

function plot_small_range()

    small = [1, 2, 3, 4, 5, 10, 20, 30, 40, 50];
    for i = small
        plot_med_model(i)
        figure
    end

    plot_med_data()
end

function plot_wide_range()

    plot_med_data()

    wide = [0, 20, 50, 100, 150, 200, 300];
    for i = wide
        figure
        plot_med_model(i)
    end

end

function plot_med_data()

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

    surface(x, y, z, col);

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

% Plot the model of a given degree
function plot_med_model(degree)

    diff_csv = strcat("../csv/diff/diff_med_", num2str(degree), ".csv");
    path_med = "../csv/ETOPO1_med.csv";

    diff_table = readtable(diff_csv);
    diff_data = diff_table{:,:};
    
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
    
    r_model = reshape(sph(:,3) + diff_data, [nt np]);

    col = r_model;
    r_model_displaced = r_model + 100000;
    
    azimuth = linspace(-pi, pi - dp, np)';
    elevation = linspace(pi/2, (-pi/2) + dt, nt)';
    
    [th, ph] = meshgrid(azimuth, elevation);
    
    [x, y, z] = sph2cart(th, ph, r_model_displaced);
    
    surface(x, y, z, col);
    
    % surface(th, ph, r)
    xlabel('th')
    ylabel('ph')
    title(strcat("Degree: ", num2str(degree)))
    % colorbar
    % colormap winter
    colormap jet
    
    axis equal
    shading interp

    CAZ = 15.350734615445207;
    CEL = 18.233845566149686;

    view(CAZ, CEL)

    axis off

    saveas(gcf, strcat("med_", num2str(degree), ".jpg"))

end

function phi = lam_to_phi(lam)

    phi = lam;

end

function theta = phi_to_the(phi)

    theta = (pi / 2) - phi;

end
% Let's use the medium data set and the information stored in the 
% diff files to plot degrees of different models

% plot_wide_range()
% plot_small_range()

small = [1, 2, 3, 4, 5, 10, 20, 30, 40, 50];
wide = [0, 20, 50, 100, 150, 200, 300];

varied = [1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 100, 150, 200, 300, 400, 500, 600];

[X, Y, Z, C] = get_XYZC_med(varied);

% with X, Y, Z, and C already computed,
% plot the animation

function [X, Y, Z, C] = plot_range(range)

    [X, Y, Z, C] = get_XYZC_med(range);

    n_degree = length(range);
    % number of frames in between two adjacent models
    n_transition = 5;
    n_frames = n_transition * n_degree;

    s = surface(X(:,:,1), Y(:,:,1), Z(:,:,1), C(:,:,1));

    CAZ = 15.350734615445207;
    CEL = 18.233845566149686;

    view(CAZ, CEL)
    shading interp
    colormap jet

    scale = [linspace(0, 1, n_transition)];

    
    for i = 1:(n_degree - 1)
        for j = 1:n_transition
            s.XData = convex_comb(X(:,:,i), X(:,:,i + 1), scale(j));
            s.YData = convex_comb(Y(:,:,i), Y(:,:,i + 1), scale(j));
            s.ZData = convex_comb(Z(:,:,i), Z(:,:,i + 1), scale(j));
            s.CData = convex_comb(C(:,:,i), C(:,:,i + 1), scale(j));
            pause(0.01)
        end 
    end

    % plot_med_data()
end

function z = convex_comb(x, y, lambda) 
    z = (1 - lambda) * x + lambda * y;
end

function plot_wide_range()

    plot_med_data();

    wide = [0, 20, 50, 100, 150, 200, 300];
    for i = wide
        figure
        plot_med_model(i);
    end

end

function [X, Y, Z, C] = get_XYZC_med(range) 

    nt = 540;
    np = 1080;
    N = nt * np;
    % small = [1, 2, 3, 4, 5, 10, 20, 30, 40, 50];

    X = zeros(nt, np, length(range));
    Y = zeros(nt, np, length(range));
    Z = zeros(nt, np, length(range));
    C = zeros(nt, np, length(range));

    for i = 1:length(range) 
        [x, y, z, c] = get_xyzc_med(range(i));
        X(:,:,i) = x; 
        Y(:,:,i) = y;
        Z(:,:,i) = z;
        C(:,:,i) = c;
    end

end 

function [x, y, z, col] = get_xyzc_med(degree) 

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
    r_model_displaced = r_model + 2e5;
    
    azimuth = linspace(-pi, pi - dp, np)';
    elevation = linspace(pi/2, (-pi/2) + dt, nt)';
    
    [th, ph] = meshgrid(azimuth, elevation);
    
    [x, y, z] = sph2cart(th, ph, r_model_displaced);

end
% Plot the model of a given degree
function s = plot_med_model(degree)

    [x, y, z, col] = get_xyzc_med(degree)
    
    s = surface(x, y, z, col);
    
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
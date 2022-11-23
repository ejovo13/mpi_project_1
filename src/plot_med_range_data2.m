function plot_med_range_data2(X, Y, Z, C, range, filename)

    nt = 540;
    np = 1080;
    N = nt * np;
    n_degree = length(range);
    % number of frames in between two adjacent models
    n_transition = 20;
    n_frames = n_transition * n_degree;

    s = surface(X(:,:,1), Y(:,:,1), Z(:,:,1), C(:,:,1));
    shading interp
    axis equal
    axis off 
    
    CAZ = -1.351331452423071e+02;
    CEL = -37.606821401124506;

    
    view(CAZ, CEL)
    colormap jet

    % Animation

    scale = [linspace(0, 1, n_transition)];
    title(strcat("L", num2str(range(1))))

    gif(filename)
    
    for i = 1:(n_degree - 1)
        for j = 1:n_transition
            s.XData = convex_comb(X(:,:,i), X(:,:,i + 1), scale(j));
            s.YData = convex_comb(Y(:,:,i), Y(:,:,i + 1), scale(j));
            s.ZData = convex_comb(Z(:,:,i), Z(:,:,i + 1), scale(j));
            s.CData = convex_comb(C(:,:,i), C(:,:,i + 1), scale(j));
            gif
        end
        title(strcat("L", num2str(range(i + 1))))
    end

end

function z = convex_comb(x, y, lambda) 
    z = (1 - lambda) * x + lambda * y;
end
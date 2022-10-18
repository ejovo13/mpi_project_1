path_med = '/home/ejovo/MAIN/S7/PPAR/project_1/csv/ETOPO1_med.csv';
path_tiny = '/home/ejovo/MAIN/S7/PPAR/project_1/csv/ETOPO1_tiny.csv';
path_unit = '/home/ejovo/MAIN/S7/PPAR/project_1/csv/ETOPO1_unit.csv';
path_small = '/home/ejovo/MAIN/S7/PPAR/project_1/csv/ETOPO1_small.csv';
path_small_prediction = '/home/ejovo/MAIN/S7/PPAR/project_1/csv/prediction.csv';

path_pred = '/home/ejovo/MAIN/S7/PPAR/project_1/build/src/prediction.csv';


% load(path_unit);
% tmp = ETOPO1_unit;

% load(path_tiny);
% tmp = ETOPO1_tiny;

% load(path_med);
% tmp = ETOPO1_med;

% load(path_small);
% tmp = ETOPO1_small;

% load(path_small_prediction)
% tmp = prediction;

load(path_pred);
tmp = prediction;

sph = zeros(size(tmp));

sph(:,1) = phi_to_the(tmp(:,2));
sph(:,2) = lam_to_phi(tmp(:,1));
sph(:,3) = tmp(:,3);

% Now that I have the coordinates in spherical coordinates, let's convert to cartesian and then plot the surface
% we actually have to intelligently arange the array. So here we go

% first let's get the number of phi and theta
nt = sqrt(length(sph) / 2);
np = 2 * nt;
dt = pi / nt;
dp = 2 * pi / np;

r = reshape(sph(:,3), [nt np]);

col = r;
r = r + 100000; % arbitrary height to show differences in altitutde

% r = r + 1e9;


azimuth = linspace(-pi, pi - dp, np)';

% disp("azimuth:")
% disp(azimuth);
elevation = linspace(pi/2, (-pi/2) + dt, nt)';
% disp('elevation: ')
% disp(elevation)

[th, ph] = meshgrid(azimuth, elevation);

[x, y, z] = sph2cart(th, ph, r);

surface(x, y, z, col);

% surface(th, ph, r)
xlabel('th')
ylabel('ph')
colorbar
% colormap winter
colormap jet

axis equal
shading interp

function phi = lam_to_phi(lam)

    phi = lam;

end

function theta = phi_to_the(phi)

    theta = (pi / 2) - phi;

end

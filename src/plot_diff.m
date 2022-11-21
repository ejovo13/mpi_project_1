% Let's go ahead and plot the differences for a small dataset

path_diff_5 = '/home/ejovo/MAIN/S7/PPAR/project_1/csv/diff/diff_5_small.csv';
path_diff_10 = '/home/ejovo/MAIN/S7/PPAR/project_1/csv/diff/diff_10_small.csv';
path_diff_20 = '/home/ejovo/MAIN/S7/PPAR/project_1/csv/diff/diff_20_small.csv';
path_diff_50 = '/home/ejovo/MAIN/S7/PPAR/project_1/csv/diff/diff_50_small.csv';
path_diff_100 = '/home/ejovo/MAIN/S7/PPAR/project_1/csv/diff/diff_100_small.csv';
path_diff_200 = '/home/ejovo/MAIN/S7/PPAR/project_1/csv/diff/diff_200_small.csv';
path_diff_300 = '/home/ejovo/MAIN/S7/PPAR/project_1/csv/diff/diff_300_small.csv';
path_diff_400 = '/home/ejovo/MAIN/S7/PPAR/project_1/csv/diff/diff_400_small.csv';
% path_diff_500 = '/home/ejovo/MAIN/S7/PPAR/project_1/csv/diff/diff_500_small.csv';



path_small = '/home/ejovo/MAIN/S7/PPAR/project_1/csv/ETOPO1_small.csv';

% data_diff = readtable(path_diff_5);
% data_diff = readtable(path_diff_10);
% data_diff = readtable(path_diff_20);
% data_diff = readtable(path_diff_50);
% data_diff = readtable(path_diff_100);
% data_diff = readtable(path_diff_200);
% data_diff = readtable(path_diff_300);
% data_diff = readtable(path_diff_400);
% data_diff = readtable(path_diff_500);
data_small = readtable(path_small);
sph = zeros(size(data_small));

sph(:,1) = phi_to_the(data_small{:,2});
sph(:,2) = lam_to_phi(data_small{:,1});
sph(:,3) = data_small{:,3};

diff = data_diff{:,1};

% first let's get the number of phi and theta
nt = sqrt(length(sph) / 2);
np = 2 * nt;
dt = pi / nt;
dp = 2 * pi / np;

r = reshape(sph(:,3), [nt np]);

% col = r;
col = reshape(diff, [nt np]);


r = r + 100000; % arbitrary height to show differences in altitutde
azimuth = linspace(-pi, pi - dp, np)';
elevation = linspace(pi/2, (-pi/2) + dt, nt)';

[th, ph] = meshgrid(azimuth, elevation);

% [x, y, z] = sph2cart(th, ph, r);
[x, y, z] = sph2cart(th, ph, ones(size(r)));

surface(x, y, z, col);

xlabel('th')
ylabel('ph')
% clim([-1000, 1000])
% caxis([-1000, 1000])
colorbar
colormap bone

axis equal
shading interp





function phi = lam_to_phi(lam)

    phi = lam;

end

function theta = phi_to_the(phi)

    theta = (pi / 2) - phi;

end
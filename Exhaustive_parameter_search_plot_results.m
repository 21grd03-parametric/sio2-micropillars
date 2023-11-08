%% Load the data from the scan

clear
close all

disp(['------------MATLAB code developed by----------',newline, ...
      '---------Jérémy Werlé LENS 16/12/2022---------',newline, ...
      '-------------------showing--------------------',newline, ...
      '--Results of the exhaustive parameter search--'])

%*************************************************************
%       Load simulation results and Input parameters
%*************************************************************

% Load simulations results
g.Res.Rup = readmatrix('Rup.csv');
g.Res.Tup = readmatrix('Tup.csv');
g.Res.Abs = readmatrix('Abs.csv');

% Load input parameters
Input_param = readmatrix('Input_param.csv');

g.my_spacing = linspace(Input_param(1,1), Input_param(1,2), Input_param(1,3));     % [um]
g.my_height = linspace(Input_param(2,1), Input_param(2,2), Input_param(2,3));      % [um]
g.my_diameter = linspace(Input_param(3,1), Input_param(3,2), Input_param(3,3));    % [um]

% All possible combinations of structures
[g.My_Spacing, g.My_height, g.My_Diameter] = ndgrid(g.my_spacing, g.my_height, g.my_diameter);

%%
%*************************************************************
%                   Find the interpolant
%*************************************************************

for i = size(g.Res.Rup,2):-1:1 % Loop over the number of tested structures
    % computation of the average emissivity over all the wavelength
    g.Res.Abs_mean(i) = mean(g.Res.Abs(:,i));
    Rup_mean(i) = mean(g.Res.Rup(:,i));
    Tup_mean(i) = mean(g.Res.Tup(:,i));
end

% Reshape matrix of average emissivity
Abs_mean_reshape = reshape(g.Res.Abs_mean, size(g.My_Spacing));

% Interpolate Abs_mean on a finer grid
F = griddedInterpolant(g.My_Spacing, g.My_height, g.My_Diameter, Abs_mean_reshape, 'spline');

% higher number of points to have a more refined plot
nb = 60;
my_spacing = linspace(g.my_spacing(1), g.my_spacing(end), nb);
my_height = linspace(g.my_height(1), g.my_height(end), nb+1);
my_diameter = linspace(g.my_diameter(1), g.my_diameter(end), nb+2);


% Creation of an ndgrid to use the gridded interpolant
[My_Spacing, My_height, My_Diameter] = ndgrid(my_spacing, my_height, my_diameter);

% Creation of a meshgrid to use the isosurface function
[My_height2, My_Spacing2, My_Diameter2] = meshgrid(my_height, my_spacing, my_diameter);

Abs_mean_fine = F(My_Spacing, My_height, My_Diameter);

%%
%*************************************************************
%              Experimental structure (Sample 2)
%*************************************************************

% Experimental structure
H_Z = 1.6;    % in  um
D_Z = 0.4;    % in  um
S_Z = 0.2;    % in  um

Abs_mean_reshape_fine_Z = F(S_Z, D_Z, H_Z); % Spacing, height, diameter

%%
%***********************************************************************
%  Single 3D plot of the Average Emissivity isosurface contour surface
%***********************************************************************

isoemiss = 0.955; % set an emissivity level for the isosurface plot

figure;
hold on;
% Isosurface structure
fv = isosurface(My_height2, My_Spacing2, My_Diameter2, Abs_mean_fine, isoemiss);
p = patch(fv, 'FaceColor', 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.4);
isonormals(My_height2, My_Spacing2, My_Diameter2, Abs_mean_fine, p)

% Experimental structure point
scatter3(H_Z, S_Z, D_Z, 'o', 'filled', ...
                        'MarkerEdgeColor', 'red', ...
                        'MarkerFaceColor', 'g', ...
                        'SizeData', 50)

% Lines
plot3([H_Z H_Z], [my_spacing(1) my_spacing(end)], [D_Z D_Z], '--', 'Color', 'k', 'LineWidth', 2)
plot3([my_height(1) my_height(end)], [S_Z S_Z],[D_Z D_Z], '--', 'Color', 'k', 'LineWidth' ,2)
plot3([H_Z H_Z], [S_Z S_Z],[my_diameter(1) my_diameter(end)], '--', 'Color', 'k', 'LineWidth', 2)

% Light parameters
camlight
camlight(-80,-10)
lighting gouraud
view(3)

% Labels
xlabel('height [µm]')
ylabel('spacing [µm]')
zlabel('diameter [µm]')
title('average emissivity optimization')

% Legend
my_legend = legend(['isosurface avg. emiss. = ', num2str(isoemiss)], 'experimental sample');
set(my_legend, 'Position',[0.7 0.8 0.1 0.1]);

% Graph parameters
grid
axis([g.my_height(1), g.my_height(end), ...
      g.my_spacing(1), g.my_spacing(end), ...
      g.my_diameter(1), g.my_diameter(end)])

saveas(gcf, 'emissivity_isosurface.png')
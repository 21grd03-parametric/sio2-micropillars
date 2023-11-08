%% Exhaustive parameter search

clear
close all
clc

disp(['---------------MATLAB code developed by---------------',newline, ...
      '-------------Jérémy Werlé LENS 16/12/2022-------------',newline, ...
      '--------------Exhaustive parameter search-------------',newline, ...
      '--using RETICOLO developed by J.P.Hugonin, P.Lalanne--'])


%*************************************************************
%                      Initialization
%*************************************************************

% /!\ Declaration of the useful path to be changed by the user /!\
path.simu = pwd;
path.reticolo = 'insert\your\path\to\RETICOLO V9\reticolo_allege_v9';
cd(path.simu);

addpath(genpath(path.reticolo));        % add the Reticolo functions
retio;                                  % clear old results

%*************************************************************
%                 Defintion of the wavelength
%*************************************************************

My_Param.wave.L_min = 8;                  % [µm]
My_Param.wave.L_max = 13;                 % [µm]
My_Param.wave.D_L = 0.1;                  % [µm]
My_Param.wave.Lam = My_Param.wave.L_min:My_Param.wave.D_L:My_Param.wave.L_max;
Lam = My_Param.wave.Lam;                  % [µm]

%*************************************************************
%                 Import Material properties
%*************************************************************

% Refractive index coming fromrfractive index info:https://refractiveindex.info/?shelf=main&book=SiO2&page=Kischkat
Name = 'SiO2_J. Kischkat, S. Peters.txt';

data = load(fullfile(path.simu,Name));
xq = My_Param.wave.L_min : My_Param.wave.D_L : My_Param.wave.L_max; % interpolation vector
nSiO2 = interp1(data(:,1), data(:,2), xq, 'linear').';
kSiO2 = interp1(data(:,1), data(:,3), xq, 'linear').';

My_Param.material.n.SiO2 = nSiO2;
My_Param.material.k.SiO2 = kSiO2;
n.air = 1;

My_Param.material.Cpx = My_Param.material.n.SiO2 + 1j * My_Param.material.k.SiO2;

My_Param.simu_type = 'loop_all';

%%
% *************************************************************
%            Geometric parameters of the structure
% *************************************************************

nb = 10;

spacing_min = 0;
spacing_max = 0.4;
nb_spacing = nb-2;

height_min = 0.1;
height_max = 4;
nb_height = nb+2;

diameter_min = 0.2;
diameter_max = 1;
nb_diameter = nb+2;

my_spacing = linspace(spacing_min,spacing_max,nb_spacing);        % [µm]
my_height = linspace(height_min,height_max,nb_height);            % [µm]
my_diameter = linspace(diameter_min,diameter_max,nb_diameter);    % [µm]


% All possible combinations of structures
[My_Spacing, My_height, My_Diameter] = ndgrid(my_spacing, my_height, my_diameter);

My_Param.struct.d = My_Diameter;
My_Param.struct.spacing = My_Spacing;
My_Param.struct.h_pillar = My_height;

My_Param.struct.thick = 2*1e3;

%*************************************************************
%            Reticolo simulation parameters
%*************************************************************

N_H = numel(My_Spacing);      % number of structures
N = length(Lam);              % number of wavelength
Tot_length = N*N_H;           % total number of simulations

My_Param.simu.h = [10 10];    % number of Fourier harmonics in 2D
My_Param.angle_theta = 0;     % incident angle of the light

angle_delta = 0;              % angle of polarization
parm = res0;                  % default parameter for the resolution
parm.sym.x = 0;               % structure symmetry under x direction
parm.sym.y = 0;               % structure symmetry under y direction
parm.sym.pol = 0;             % Polarization -1 TM 1 TE 0 for both

parm.res1.nx = 256;           % number of points to plot the structure under x direction
parm.res1.ny = 256;           % number of points to plot the structure under y direction

parm.res1.champ = 0;          % do not compute the electric and magnetic fields
parm.res1.trace = 0;          % do not plot 2D image of the layers

%%
%*************************************************************
%               Construction of the textures
%*************************************************************

% Loop over the number of structures
for j = 1:N_H
    [row,col,prof] = ind2sub(size(My_Spacing), j);

    % Definition of the geometrical parameters

    My_Param.struct.l(row,col,prof) = My_Param.struct.d(row,col,prof)+My_Param.struct.spacing(row,col,prof);
    % triangular lattice
    My_Param.struct.L(row,col,prof) = sqrt(3)*(My_Param.struct.l(row,col,prof));
    My_Param.struct.P(j,:) = [My_Param.struct.L(row, col, prof), My_Param.struct.l(row, col, prof)];

    for i = 1:N
        % Definition of the texture

        % 1st layer: air (infinite)
        My_Param.textures{1,i,j} = n.air;

        % 2nd layer: SiO2 pillar
        My_Param.textures{2,i,j}= {n.air [0,0,My_Param.struct.d(row,col,prof), My_Param.struct.d(row,col,prof), My_Param.material.Cpx(i), 20],...
                                         [-My_Param.struct.L(row,col,prof)/2, My_Param.struct.l(row,col,prof)/2, My_Param.struct.d(row,col,prof), My_Param.struct.d(row,col,prof), My_Param.material.Cpx(i), 20]};

        % 3rd layer: residual layer of glass
        My_Param.textures{3,i,j} = My_Param.material.Cpx(i);

        % 4th layer: air (infinite)
        My_Param.textures{4,i,j} = n.air;

        % Definition of the profile of thicknesses
        My_Param.profile{j} = {[0, My_Param.struct.h_pillar(row,col,prof), My_Param.struct.thick,0], 1:size(My_Param.textures,1)};
    end
    % Definition of the profile
end


%%
%*************************************************************
%               Construction of the textures
%*************************************************************

tic;
m = 1;

reverseStr = '';
fprintf('\nprogress: ');

for j = 1:N_H                               % Loop over the declared parameter

    for i = 1:N                             % Loop over the wavelength

        wavelength = Lam(i);                % Lam

        for k = size(My_Param.textures,1):-1:1
            textures{k} = My_Param.textures{k,i,j}; % extraction of the good texture
        end

        k_parallel = sin(My_Param.angle_theta(1)*pi/180); % parallel component of wave vector

        aa = res1(wavelength, My_Param.struct.P(j,:), textures, My_Param.simu.h, k_parallel(1), angle_delta, parm);
        result1 = res2(aa, My_Param.profile{j}, parm);

        % Transmission TM polarisation
        tp1(i,j) = result1.TMinc_top_transmitted;
        Tpp(i,j) = sum((tp1(i,j).efficiency(tp1(i,j).theta~=90))); % avoid pathological theta=90 value

        % Reflection TM polarisation
        rp1(i,j) = result1.TMinc_top_reflected;
        Rpp(i,j) = sum((rp1(i,j).efficiency(rp1(i,j).theta~=90)));

        % Transmission TE polarisation
        ts1(i,j) = result1.TEinc_top_transmitted;
        Tss(i,j) = sum((ts1(i,j).efficiency(ts1(i,j).theta~=90)));

        % Reflection TE polarisation
        rs1(i,j) = result1.TEinc_top_reflected;
        Rss(i,j) = sum((rs1(i,j).efficiency(rs1(i,j).theta~=90)));

        % Averageing of TE and TM for the unpolarized light
        Res.Tup(i,j) = (Tss(i,j) + Tpp(i,j))/2;
        Res.Rup(i,j) = (Rss(i,j) + Rpp(i,j))/2;
        Res.Abs(i,j) = 1-(Res.Tup(i,j) + Res.Rup(i,j));
        retio;      % clear the old intermediate results
        m = m+1;
        msg = sprintf('%d / %d', m, Tot_length);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        pause(0.001)
        clear textures
    end
    clear i
end

fprintf(' ... done\n')

toc;

%%
%*************************************************************
%                       Save Results
%*************************************************************
writematrix(g.Res.Rup, 'Rup.csv')
writematrix(g.Res.Tup, 'Tup.csv')
writematrix(g.Res.Abs, 'Abs.csv')

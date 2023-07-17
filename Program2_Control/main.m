clear all; clc;
addpath('./CosseratRods', './BigMag', './BigMag/Classes', './MatData', './Plot')

% Decide whether to plot the actuation system after solving for the
% configuration
plot_figures = true;

%% Load materials (pdms, nitinol, stainlessSteel, smpRubber, smpGlass)

% Load material states (function at the end of this script)
material = load_materials;

% Load base position (global (A) Frame) of the AGC (scripts in the Program1 folder)
load('MatData/basePos.mat');

% Load target positions in the global frame (scripts in the Program1 folder)
load('MatData/target1.mat');
load('MatData/target2.mat');
load('MatData/target3.mat');

% Select target
target_position_AFrame = target3; %target1, target2, target3

%% Rod parameters

% Initialize struct for the rubber-phase AGC
C1.f.radius_inner = 0.6e-3;
C1.f.radius_outer = 2e-3;
C1.f.length = 45e-3;
C1.f.elastic_modulus = material.smpRubber.E; % Elastic modulus
C1.f.shear_modulus = material.smpRubber.G; % Shear modulus
C1.f.density = material.smpGlass.rho; % Density


% Initialize struct for the glass-phase AGC
C1.s.radius_inner = C1.f.radius_inner;
C1.s.radius_outer = C1.f.radius_outer;
C1.s.length = C1.f.length;
C1.s.elastic_modulus = material.smpGlass.E; % Elastic modulus
C1.s.shear_modulus = material.smpGlass.G; % Shear modulus
C1.s.density = material.smpGlass.rho; % Density

% Initialize struct for the needle
C2.radius_inner = 0e-3;
C2.radius_outer = 0.5e-3;
C2.length = C1.f.length;
C2.elastic_modulus = material.stainlessSteel.E; % Elastic modulus
C2.shear_modulus = material.stainlessSteel.G; % Shear modulus
C2.density = material.stainlessSteel.rho; % Density

% Initialize struct for the magnet at the AGC tip
magnet.radius_inner = 0.75e-3;
magnet.radius_outer = 2e-3;
magnet.length = 7e-3;
magnet.density = 7000;
magnet.weight = pi * (magnet.radius_outer^2 - magnet.radius_inner^2) * magnet.length * magnet.density * [0; 0; -9.81];
magnet.moment_norm = pi * (magnet.radius_outer^2 - magnet.radius_inner^2) * magnet.length * 1.1 / (pi*4e-7);


%% Initialize Cosserat rods
% Rubber-phase AGC
C1.f.object = CosseratRod(C1.f.radius_inner, C1.f.radius_outer, ...
    C1.f.length, C1.f.elastic_modulus, C1.f.shear_modulus, C1.f.density);

% Glass-phase AGC
C1.s.object = CosseratRod(C1.s.radius_inner, C1.s.radius_outer, ...
    C1.s.length, C1.s.elastic_modulus, C1.s.shear_modulus, C1.s.density);

% Needle
C2.object = CosseratRod(C2.radius_inner, C2.radius_outer, ...
    C2.length, C2.elastic_modulus, C2.shear_modulus, C2.density);


%% Compute torque for shape configuration

% Define homogeneous transformation from the AGC body frame to actuation frame
H.p.p_ab = basePos;  % AGC body position in global frame
H.R.R_ab = roty(90); % AGC body orientation relative to global frame
H.H_ab = getTransform(H.R.R_ab, H.p.p_ab); % Body-to-Actuation frame transformation

% Transform target position to the AGC body frame
target_position_BFrame = pinv(H.H_ab) * [target_position_AFrame; 1];
target_position_BFrame = target_position_BFrame(1:3);

% Define boundary conditions for the flexible AGC
C1.f.bc.prox.p = [0;0;0];     C1.f.bc.prox.q = [1;0;0;0];  %  Proximal
C1.f.bc.dist.n = H.H_ab(1:3,1:3)' * magnet.weight;  C1.f.bc.dist.m = [0;0;0];  %  Distal

% Simulation framework (Fig.(3B)). Compute torque to deflect the flexible
% AGC and align the CTR with the target
[torque_BFrame, Y_C1f_star_BFrame, Y_C1f_BFrame] = ...
    getTorque(C1.f.object, C1.s.object, C2.object, C1.f.bc.prox.p, C1.f.bc.prox.q, C1.f.bc.dist.n, C1.f.bc.dist.m, target_position_BFrame);

%% Transform torque to global frame and obtain magnet position & dipole moment

% Transform body-frame torque to global actuation frame
torque_AFrame = H.R.R_ab * torque_BFrame;

% Extract body-frame magnet position and dipole moment
magnet.p.BFrame = Y_C1f_BFrame(1:3,end);
magnet.mu.BFrame = quat2rotm(Y_C1f_BFrame(4:7,end)') * [0; 0; -1]; % The north-pole of the tip-magnet is directed down the length of the AGC
magnet.mu.BFrame = magnet.mu.BFrame / norm(magnet.mu.BFrame);

% Transform to global actuation frame
magnet.p.AFrame = H.H_ab * [magnet.p.BFrame; 1];   magnet.p.AFrame = magnet.p.AFrame(1:3);
magnet.mu.AFrame = H.R.R_ab * magnet.mu.BFrame;
magnet.mu.AFrame = magnet.mu.AFrame * magnet.moment_norm;


%% Magnetic actuation system

% Initialize object (System name: BigMag)
BigMag = BigMag_Object;
BigMag_neutral = BigMag_Object; % Unconfigured system

% Compute system configuration based on magnet position and dipole moment,
% as well as desired torque expressed in the actuation frame
[BigMag_positions_final, BigMag_currents_final] = getBigMagConfiguration(...
    BigMag, magnet.p.AFrame, magnet.mu.AFrame, ...
    torque_AFrame);

% Print required configuration
disp(['BigMag positions final: ' num2str(BigMag_positions_final')])
disp(['BigMag currents final: ' num2str(BigMag_currents_final')])

% Update the actuation system object
BigMag = BigMag.update(BigMag_currents_final, BigMag_positions_final);


%% Compute centerline polynomials of the AGC and CTR, write to file

% Initialize rubber-phase AGC model with boundary conditions
C1.f.model = ConcentricTubeRobot([C1.f.object]); % A CTR comprising one rod (rubber SMP)
C1.f.bc.prox.p = [0;0;0];     C1.f.bc.prox.q = [1;0;0;0];     C1.f.bc.y0 = [C1.f.bc.prox.p; C1.f.bc.prox.q; zeros(6,1)];
C1.f.bc.dist.n = H.H_ab(1:3,1:3)' * magnet.weight;      C1.f.bc.dist.m = torque_BFrame;     C1.f.bc.yL = [C1.f.bc.dist.n; C1.f.bc.dist.m];

% Solve BVP to obtain the shape solution of the AGC
Y_C1f_BFrame = C1.f.model.solve(C1.f.bc.y0, C1.f.bc.yL, 1e-8, 30);  %  Local frame
Y_C1f_AFrame = H.H_ab * [Y_C1f_BFrame(1:3,:); ones(1, size(Y_C1f_BFrame,2))];  % Global frame

% Fit centerline polynomial to the shape solution of the AGC
p1s_BFrame = getPolynomialFromShape(Y_C1f_BFrame, C1.f.length, 4);  % Polynomial of order 4

% Write the AGC shape polynomial to file
p1s_AFrame = H.H_ab * [p1s_BFrame; sym(1)];
p1s_AFrame = p1s_AFrame(1:3);
matlabFunction(p1s_AFrame,'File','polynomial_C1s_global_simulation');

% Compute pre-curvature of the AGC shape polynomial and assign to the glass-phase AGC
[frenetSerret_R_BFrame, frenetSerret_u_BFrame] = getCurvature(p1s_BFrame, 0);
C1.s.object = C1.s.object.setPrecurvatureHandle(frenetSerret_u_BFrame);

% Initialize CTR model
C12.model = ConcentricTubeRobot([C1.s.object, C2.object]);

% Proximal CTR boundary conditions
C12.bc.prox.p = [0;0;0];    C12.bc.prox.q = rotm2quat(frenetSerret_R_BFrame(0))';  C12.bc.prox.theta = 0; % Proximal (known)
% Proximal CTR optimization vector (unknown material state values)
C12.guess.n0 = [0;0;0];    C12.guess.mbxy0 = [0;0];   C12.guess.mb1z0 = 0;   C12.guess.mb2z0 = 0; % Proximal (guess)
% Proximal CTR state vector
C12.y0 = [C12.bc.prox.p; C12.bc.prox.q; C12.guess.n0; C12.guess.mbxy0; C12.guess.mb1z0; C12.guess.mb2z0; C12.bc.prox.theta]; % Proximal
% Distal CTR boundary conditions
C12.bc.dist.nL = C1.f.bc.dist.n;   C12.bc.dist.mbxy = torque_BFrame(1:2);   C12.bc.dist.mb1z = torque_BFrame(3); C12.bc.dist.mb2z = 0; % Distal

% Solve BVP to obtain the shape solution of the CTR
Y_C12s_BFrame = C12.model.solve(C12.y0, [C12.bc.dist.nL; C12.bc.dist.mbxy; C12.bc.dist.mb1z; C12.bc.dist.mb2z], 1e-8, 60);
Y_C12s_AFrame = H.H_ab * [Y_C12s_BFrame(1:3,:); ones(1, size(Y_C12s_BFrame,2))];

% Write the CTR shape polynomial to file
P12s_BFrame = getPolynomialFromShape(Y_C12s_BFrame, C1.f.length, 4);  % Polynomial of order 4
P12s_AFrame = H.H_ab * [P12s_BFrame; sym(1)];
P12s_AFrame = P12s_AFrame(1:3);
matlabFunction(P12s_AFrame,'File','polynomial_C12s_global_simulation');



%% Draw
if (plot_figures)
    
    targetSize = 40;
    cameraPosition = [0.1, 0.2, 0.1];
    
    % Neutral setup
    BigMag_neutral.draw(1.2); hold on;
    BigMag_neutral.draw_field(0:10e-3:30e-3, 0:30:360, 0:30:360,3);
    plot3(target1(1), target1(2), target1(3), 'k.', 'MarkerSize', targetSize);
    plot3(target2(1), target2(2), target2(3), 'k.', 'MarkerSize', targetSize);
    plot3(target3(1), target3(2), target3(3), 'k.', 'MarkerSize', targetSize);
    plot_rod(H.H_ab, Y_C1f_star_BFrame, 2e-3,0.5, [1 0 0])
    plot3(Y_C1f_AFrame(1,1), Y_C1f_AFrame(2,1), Y_C1f_AFrame(3,1), 'k.', 'MarkerSize', 35);
    campos(cameraPosition)
    
    % Reconfigured BigMag
    BigMag.draw(1.2); hold on;
    BigMag_neutral.draw_field(0:10e-3:30e-3, 0:30:360, 0:30:360,3);
    plot3(target1(1), target1(2), target1(3), 'k.', 'MarkerSize', targetSize);
    plot3(target2(1), target2(2), target2(3), 'k.', 'MarkerSize', targetSize);
    plot3(target3(1), target3(2), target3(3), 'k.', 'MarkerSize', targetSize);
    plot_rod(H.H_ab, Y_C1f_star_BFrame, 2e-3,0.5, [1 0 0])
    plot3(Y_C1f_AFrame(1,1), Y_C1f_AFrame(2,1), Y_C1f_AFrame(3,1), 'k.', 'MarkerSize', 35);
    campos(cameraPosition)
    
    % Actuating BigMag
    BigMag.draw(1.2); hold on;
    BigMag.draw_field(0:10e-3:30e-3, 0:30:360, 0:30:360,3);
    plot3(target1(1), target1(2), target1(3), 'k.', 'MarkerSize', targetSize);
    plot3(target2(1), target2(2), target2(3), 'k.', 'MarkerSize', targetSize);
    plot3(target3(1), target3(2), target3(3), 'k.', 'MarkerSize', targetSize);
    plot_rod(H.H_ab, Y_C1f_BFrame, 2e-3,0.5, [1 0 0])
    plot3(Y_C1f_AFrame(1,1), Y_C1f_AFrame(2,1), Y_C1f_AFrame(3,1), 'k.', 'MarkerSize', 35);
    campos(cameraPosition)
    
    % Cooling sheath
    figure
    BigMag.draw_field(0:10e-3:30e-3, 0:30:360, 0:30:360,3);
    plot3(target1(1), target1(2), target1(3), 'k.', 'MarkerSize', targetSize);
    plot3(target2(1), target2(2), target2(3), 'k.', 'MarkerSize', targetSize);
    plot3(target3(1), target3(2), target3(3), 'k.', 'MarkerSize', targetSize);
    plot_rod(H.H_ab, Y_C1f_BFrame, 2e-3,0.5, [0 0 1])
    plot3(Y_C1f_AFrame(1,1), Y_C1f_AFrame(2,1), Y_C1f_AFrame(3,1), 'k.', 'MarkerSize', 35);
    xlabel('$x$ (m)', 'Interpreter', 'Latex');  ylabel('$y$ (m)', 'Interpreter', 'Latex');  zlabel('$z$ (m)', 'Interpreter', 'Latex');
    set(gca,'CameraPosition',[2 2 2], 'FontSize', 35);
    campos(cameraPosition)
    
    % Final
    figure
    target3_local = pinv(H.H_ab)*[target3;1];
    needle_shape_tip_local = [Y_C12s_BFrame(:,end) [target3_local(1:3); Y_C12s_BFrame(4:end,end)]];
    BigMag.draw_field(0:10e-3:30e-3, 0:30:360, 0:30:360,3);
    plot3(target1(1), target1(2), target1(3), 'k.', 'MarkerSize', targetSize);
    plot3(target2(1), target2(2), target2(3), 'k.', 'MarkerSize', targetSize);
    plot3(target3(1), target3(2), target3(3), 'k.', 'MarkerSize', targetSize);
    plot_rod(H.H_ab, Y_C12s_BFrame, 2e-3,0.5, [0 0 1])
    plot_rod(H.H_ab, needle_shape_tip_local, 0.5e-3, 1, [0 0 0])
    plot3(Y_C1f_AFrame(1,1), Y_C1f_AFrame(2,1), Y_C1f_AFrame(3,1), 'k.', 'MarkerSize', 35);
    xlabel('$x$ (m)', 'Interpreter', 'Latex');  ylabel('$y$ (m)', 'Interpreter', 'Latex');  zlabel('$z$ (m)', 'Interpreter', 'Latex');
    set(gca,'CameraPosition',[2 2 2], 'FontSize', 35);
    campos(cameraPosition)
    
end










%% Functions
% Load in some material parameters
function material = load_materials()
material.pdms.E = 2e6;
material.pdms.nu = 0.4;
material.pdms.G = material.pdms.E/(2*(1+material.pdms.nu));
material.pdms.rho = 1000;

material.nitinol.E = 50e9;
material.nitinol.nu = 0.33;
material.nitinol.G = material.nitinol.E/(2*(1+material.nitinol.nu));
material.nitinol.rho = 6450;

material.stainlessSteel.E = 193e9;
material.stainlessSteel.nu = 0.3;
material.stainlessSteel.G = material.stainlessSteel.E/(2*(1+material.stainlessSteel.nu));
material.stainlessSteel.rho = 7500;

material.smpRubber.E = 7.02e6;
material.smpRubber.nu = 0.35;
material.smpRubber.G = material.smpRubber.E/(2*(1+material.smpRubber.nu));
material.smpRubber.rho = 2000;

material.smpGlass.E = 0.5e9;
material.smpGlass.nu = 0.35;
material.smpGlass.G = material.smpGlass.E/(2*(1+material.smpGlass.nu));
material.smpGlass.rho = 2000;
end




% Input: Shape solution (Y), length of rod, order of polynomial.
function P_s_sym = getPolynomialFromShape(Y, length, order)

% Length of discretized segments
ds = length / (size(Y,2)-1);

% Number of centerline points
numpoints = length/ds;

% Base position
p0 = Y(1:3,1);

% Initialize matrices
mat = NaN(numpoints*3, order*3);
vec = NaN(numpoints*3, 1);

% Fill matrices
for i = 1:numpoints
    % Arc-length
    s = i * ds;
    % Position
    p_i = Y(1:3,i+1);
    % Rows to fill
    rows = (i-1)*3+1:(i-1)*3+3;
    % Fill position vector
    vec(rows,1) = p_i - p0;
    % Fill matrix
    for j = 1:order
        cols = (j-1)*3+1:(j-1)*3+3;
        mat(rows, cols) = s^j * eye(3);
    end
end

% Compute coefficients
coeffs = pinv(mat)*vec;

% Form symbolic polynomial
syms s;
% Initialize
mat = NaN(3,size(coeffs,1));
vec = sym(NaN(size(mat,2),1));
for i = 1:order
    % Columns
    cols = (i-1)*3+1:(i-1)*3+3;
    % Fill matrix of coefficients
    mat(:,cols) = diag(coeffs(cols,1));
    % Fill vector of powers of s
    vec(cols,1) = sym([1;1;1]) * s^i;
end

% Polynomial
P_s_sym = sym(p0) + sym(mat) * vec;

%P_s = matlabFunction(P_s_sym);
end




% Input: rotation matrix (R) and position (p)
function T = getTransform(R,p)
T = [R p; zeros(1,3) 1];
end



% Input: 3D vector (a)
function S = getSkew(a)
S = [0 -a(3) a(2);
    a(3) 0 -a(1);
    -a(2) a(1) 0];
end

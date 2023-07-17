%% Import voxels
% Read voxels in Actuation Frame (AFrame)
voxels_AFrame = readmatrix('voxels_workspace.csv');

% Plot allvoxels
figure; 
plot3(voxels_AFrame(1,:), voxels_AFrame(2,:), voxels_AFrame(3,:),'r.'); 
grid on; hold on; axis equal; 
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)')

%% Compute base position of the AGC in the workspace
% Search range of voxel positions to find base position of the AGC
rangeX0 = [-0.06 -0.055];
rangeY0 = [-0.005 0.005];
rangeZ0 = [-0.005 0.005];
basePos = getTargetPosition(voxels_AFrame, rangeX0, rangeY0, rangeZ0);
plot3(basePos(1), basePos(2), basePos(3), 'k.', 'MarkerSize', 25)

%% Find Target Positions in the workspace

% Target 1
rangeX1 = [0.03 0.06];
rangeY1 = [-0.02 0.02];
rangeZ1 = [0.02 0.04];
target1 = getTargetPosition(voxels_AFrame, rangeX1, rangeY1, rangeZ1);
plot3(target1(1), target1(2), target1(3), 'k.', 'MarkerSize', 25)

% Target 2
rangeX2 = [0.03 0.06];
rangeY2 = [-0.04 -0.02];
rangeZ2 = [0 0.015];
target2 = getTargetPosition(voxels_AFrame, rangeX2, rangeY2, rangeZ2);
plot3(target2(1), target2(2), target2(3), 'k.', 'MarkerSize', 25)

% Target 3
rangeX3 = [0.02 0.04];
rangeY3 = [0.03 0.05];
rangeZ3 = [-0.04 0];
target3 = getTargetPosition(voxels_AFrame, rangeX3, rangeY3, rangeZ3);
plot3(target3(1), target3(2), target3(3), 'k.', 'MarkerSize', 25)


%% Write target coordinates to file
Folder = cd;
Folder = fullfile(Folder, '..');
save(fullfile([Folder '\Program2_Control\MatData'], 'basePos.mat'), 'basePos')
save(fullfile([Folder '\Program2_Control\MatData'], 'target1.mat'), 'target1')
save(fullfile([Folder '\Program2_Control\MatData'], 'target2.mat'), 'target2')
save(fullfile([Folder '\Program2_Control\MatData'], 'target3.mat'), 'target3')


%% Functions
% Get central position of voxels within specified range of X- Y-
% Z-positions
function P = getTargetPosition(voxels, rangeX, rangeY, rangeZ)
% Extract minimum and maximum value of xyz-component
xmin = min(rangeX); xmax = max(rangeX);
ymin = min(rangeY); ymax = max(rangeY);
zmin = min(rangeZ); zmax = max(rangeZ);

% Initialize matrix of voxel positions within range
centerVoxels = [];
% Loop over all voxels
for i = 1:size(voxels,2)
    % Extract voxel coordinates
    x = voxels(1,i);
    y = voxels(2,i);
    z = voxels(3,i);
    % If voxel coordinates are within range, store the voxel position
    if ((x >= xmin && x <= xmax) && (y >= ymin && y <= ymax) && (z >= zmin && z <= zmax))
        centerVoxels = [centerVoxels [x;y;z]];
    end
end

% Compute central position coordinate of the voxels
P = mean(centerVoxels, 2);
end
function plot_rod(T, shape, radius_outer, faceAlpha, color)

% Get ds
p0 = shape(1:3,1);  p1 = shape(1:3,2);
ds = norm(p1-p0);

% Get basic cylinder
[X0, Y0, Z0] = cylinder(radius_outer);

% Scale length of cylinder with ds
Z0 = 1.1*ds * Z0;

% Loop over shape
for i = 2 : size(shape,2)
    
    % Pose of desired cylinder
    pi = shape(1:3,i-1);
    Ri = quat2rotm(shape(4:7,i-1)');
    
    % Copy basic cylinder as new
    X = X0; Y = Y0; Z = Z0;
    
    % Rotate basic cylinder
    for j = 1:size(X0,1)
        for k = 1:size(X0,2)
            pos = [X(j,k); Y(j,k); Z(j,k)];
            pos = Ri * pos;
            X(j,k) = pos(1);
            Y(j,k) = pos(2);
            Z(j,k) = pos(3);
        end
    end
    
    % Translate cylinder
    X = X + pi(1);
    Y = Y + pi(2);
    Z = Z + pi(3);
    
    % To global frame
    row1 = T * [X(1,:); Y(1,:); Z(1,:); ones(1,size(X,2))];
    row2 = T * [X(2,:); Y(2,:); Z(2,:); ones(1,size(X,2))];
    X = [row1(1,:); row2(1,:)];
    Y = [row1(2,:); row2(2,:)];
    Z = [row1(3,:); row2(3,:)];
    
    % Plot
    mesh(X,Y,Z,'FaceAlpha',faceAlpha,'EdgeColor','none', 'FaceColor',color);
    hold on;
end





end


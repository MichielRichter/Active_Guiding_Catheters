classdef BigMag_Object
    %BIGMAG Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        coils;
        angles;
        currents;
    end
    
    methods
        function obj = BigMag_Object()
            
            % Add necessary folders to path
            addpath('./BigMag/Classes', './BigMag/Classes/Fields');
            
            %Values obtained from BigMagArrayForward.cpp
            R1 = rotz(180) * roty(-34.21);
            R2 = roty(-33.22);
            R3 = rotz(180)*roty(33);
            R4 = roty(33.14);
            R5 = roty(-90);
            R6 = roty(90);
            
            %Values obtained from BigMagArrayForward.cpp
            p1 = R1 * [-0.1021; 0; 0];
            p2 = R2 * [-0.1033; 0; 0];
            p3 = R3 * [-0.101;  0; 0];
            p4 = R4 * [-0.103;  0; 0];
            p5 = R5 * [-0.0855; 0; 0];
            p6 = R6 * [-0.0855; 0; 0];
            
            % Make coils
            c0 = IroncoreCoil(R1, p1, '0'); %bottom left
            c1 = IroncoreCoil(R2, p2, '1'); %bottom right
            c2 = IroncoreCoil(R3, p3, '2'); %top left
            c3 = IroncoreCoil(R4, p4, '3'); %top right
            c4 = AircoreCoil(R5, p5, '4');  %bottom
            c5 = AircoreCoil(R6, p6, '5');  %top
            
            % Assign coils to BigMag
            obj.coils = {c0 c1 c2 c3 c4 c5};
            
            % Assign initial motor angles [degrees]
            obj.angles = {0 0};
            obj.currents = {0 0 0 0 0 0};
            
        end
        
        
        % I = [I1;I2;...], gamma = [top; bottom] in DEGREES
        function obj = update(obj, I, gamma)
            % Update currents
            obj.currents{1} = I(1);
            obj.currents{2} = I(2);
            obj.currents{3} = I(3);
            obj.currents{4} = I(4);
            obj.currents{5} = I(5);
            obj.currents{6} = I(6);
            
            % Update positions
            T_bottom = [rotz(gamma(1)) zeros(3,1); zeros(1,3) 1];
            T_top = [rotz(gamma(2)) zeros(3,1); zeros(1,3) 1];
            
            % Coil IDs
            top_ids = ["2","3","5"];
            bottom_ids = ["0","1","4"];
            
            % Update coil poses
            for i = 1:size(obj.coils,2)
                % Get coil ID
                coil_id = obj.coils{i}.id;
                
                % Check if ID is part of bottom hemisphere of coils
                if ismember(coil_id, bottom_ids)
                    obj.coils{i}.transform = T_bottom * obj.coils{i}.transform_init;
                    % Else check if ID is part of top hemisphere of coils
                elseif ismember(coil_id, top_ids)
                    obj.coils{i}.transform = T_top * obj.coils{i}.transform_init;
                end
            end
        end
        
        
        % Get magnetic field at specified locations
        % P = [ [p1x; p1y; p1z] [p2x; p2y; p2z] ... ]
        function fields = getFieldAtPositions(obj, P)
            
            % Initialize
            fields = zeros(size(P));
            
            % Loop over positions
            for i = 1 : size(P,2)
                % Position
                p = P(:,i);
                
                % Loop over coils
                for j = 1 : size(obj.coils, 2)
                    fields(:,i) = fields(:,i) + obj.coils{j}.getUnitFieldatLocation(p) * ...
                        obj.currents{j};
                end
            end
        end
        
        % Get magnetic field gradient at specified location
        function fieldGradients = getFieldGradientAtPosition(obj, p)
            
            % Initialize
            fieldGradients = zeros(3,3);
            
            % Loop over coils
            for i = 1 : size(obj.coils, 2)
                fieldGradient = obj.coils{i}.getUnitGradientatLocation(p) * ...
                    obj.currents{i};
                fieldGradients = fieldGradients + fieldGradient;
            end
        end
        
        
        % Get unit-field map (3x6 matrix)
        function unitFieldMap = getUnitFieldMap(obj, p)  
            % Initialize
            unitFieldMap = NaN(3,6);
            % Fill
            for i = 1 : size(obj.coils, 2)
                field = obj.coils{i}.getUnitFieldatLocation(p);
                unitFieldMap(:,i) = field;
            end  
        end
        
        % Get unit-gradient map (5x6 matrix)
        function unitGradientMap = getUnitGradientMap(obj, p)
           % Initialize
           unitGradientMap = NaN(5,6); 
           % Fill
           for i = 1 : size(obj.coils, 2)
                gradient = obj.coils{i}.getUnitGradientatLocation(p);
                unitGradientMap(:,i) = [gradient(1,1); gradient(1,2); gradient(1,3); gradient(2,2); gradient(2,3)];
           end  
        end
               
        
        % Draw BigMag
        function h = draw(obj, scale)
            
            h = figure('units','normalized','outerposition',[0 0 1 1], 'Name','BigMag');
            
            % Loop over coils
            for i = 1 : size(obj.coils,2)
                coil = obj.coils{i};
                T = coil.transform;
                id = coil.id;
                opaque = 1; % Boolean to make coil frame opaque
                
                obj.plotFrame(T, scale, id, opaque, h);
            end
            
            xlabel('$x$ (m)', 'Interpreter', 'Latex');
            ylabel('$y$ (m)', 'Interpreter', 'Latex');
            zlabel('$z$ (m)', 'Interpreter', 'Latex');
            set(gca,'CameraPosition',[2 2 2], 'FontSize', 35);
            
            lim = 0.2;
            axis square
            xlim([-lim lim]);
            ylim([-lim lim]);
            zlim([-lim lim]);
            campos([1.5655    2.4780    0.2806])
            camlight;
            lighting flat;
            hold off;
        end
        
        % Draw BigMag field
        function h = draw_field(obj, radius, angles_horizontal, angles_vertical, factor_scale)
            
            % Draw field arrows at radial and angular positions
            % radius - 1xN vector of radial positions
            % angles_horizontal - 1xM vector of horizontal angular positions [deg]
            % angles_vertical - 1xK vector of vertical angular positions [deg]
            
            % Make mesh of XYZ coordinates
            [Rad, A, B] = meshgrid(radius, angles_horizontal, angles_vertical);
            
            % Initialize matrices
            pos_X = NaN(size(Rad));
            pos_Y = NaN(size(Rad));
            pos_Z = NaN(size(Rad));
            field_X = NaN(size(Rad));
            field_Y = NaN(size(Rad));
            field_Z = NaN(size(Rad));
            
            % Maximum field norm
            field_norm_max = 0;
            % Loop
            for i1 = 1 : size(Rad,1)
                for i2 = 1 : size(Rad,2)
                    for i3 = 1 : size(Rad,3)
                        % Get spherical position
                        radius = Rad(i1, i2, i3);
                        angle_horizontal = A(i1, i2, i3);
                        angle_vertical = B(i1, i2, i3);
                        % Compute cartesian positions
                        pos_x = radius * sind(angle_vertical) * cosd(angle_horizontal);
                        pos_y = radius * sind(angle_vertical) * sind(angle_horizontal);
                        pos_z = radius * cosd(angle_vertical);
                        % Compute field
                        field = obj.getFieldAtPositions([pos_x; pos_y; pos_z]) * 1e3; % Tesla to milliTesla
                        % Compute field norm
                        field_norm = norm(field);
                        % Fill
                        pos_X(i1,i2,i3) = pos_x;
                        pos_Y(i1,i2,i3) = pos_y;
                        pos_Z(i1,i2,i3) = pos_z;
                        field_X(i1,i2,i3) = field(1);
                        field_Y(i1,i2,i3) = field(2);
                        field_Z(i1,i2,i3) = field(3);
                        if (field_norm > field_norm_max)
                            field_norm_max = field_norm;
                        end
                    end
                end
            end
            
            hold on;
            q = quiver3(pos_X, pos_Y, pos_Z, field_X, field_Y, field_Z, ...
                'AutoScale', 'on', 'AutoScaleFactor', factor_scale, ...
                'LineWidth', 1.25);
            obj.color_quiver_plot(q,'jet','mT',1);
            axis equal;
        end
        
    end
    
    
    
    
    
    
    methods (Access = private)
        
        function plotFrame(obj, T, scale, id, opaque, fig)
            
            figure(fig)
            
            % Get orientation and position
            R = T(1:3,1:3);
            p = T(1:3,4);
            
            % Make vectors for plotting coil frame
            x = p + R*[0.02 0 0]';
            y = p + R*[0 0.02 0]';
            z = p + R*[0 0 0.02]';
            
            % Plot
            if norm(p) ~= 0
                
                % Vertical coils are air-cored coils and are drawn differently
                vertical_coil_ids = ["4","5"];
                if ismember(id,vertical_coil_ids)
                    [X_core,Y_core,Z_core] = cylinder(scale*0.01,40); Z_core = Z_core/scale;
                    [X_wire1,Y_wire1,Z_wire1] = cylinder(scale*0.02,40); Z_wire1 = Z_wire1/scale;
                    [X_wire2,Y_wire2,Z_wire2] = cylinder(scale*0.016,40); Z_wire2 = Z_wire2/scale;
                    [X_wire3,Y_wire3,Z_wire3] = cylinder(scale*0.013,40); Z_wire2 = Z_wire3/scale;
                    core_length = scale*0.06; Z_core = core_length*Z_core;
                    wire_length_1 = scale*0.020;  Z_wire1 = wire_length_1*Z_wire1;
                    wire_length_2 = scale*0.020;  Z_wire2 = wire_length_2*Z_wire2;
                    wire_length_3 = scale*0.020;  Z_wire3 = wire_length_3*Z_wire3;
                    core_disp = R*scale*[-0.002;0;0];
                    
                    % Plot iron-cored diagonal coils
                else
                    % Mesh for iron core
                    [X_core,Y_core,Z_core] = cylinder(scale*0.01,40);
                    % Mesh for first layer of copper wiring
                    [X_wire1,Y_wire1,Z_wire1] = cylinder(scale*0.02,40);
                    % Mesh for second layer of copper wiring
                    [X_wire2,Y_wire2,Z_wire2] = cylinder(scale*0.016,40);
                    % Mesh for third layer of copper wiring
                    [X_wire3,Y_wire3,Z_wire3] = cylinder(scale*0.013,40);
                    % Assign length of the core
                    core_length = scale*0.06; Z_core = core_length*Z_core;
                    % Assign lengths of the wires
                    wire_length_1 = scale*0.045;  Z_wire1 = wire_length_1*Z_wire1;
                    wire_length_2 = scale*0.045;  Z_wire2 = wire_length_2*Z_wire2;
                    wire_length_3 = scale*0.045;  Z_wire3 = wire_length_3*Z_wire3;
                    % Displace core slightly for plotting
                    core_disp = R*scale*[-0.002;0;0];
                end
                
                % If wires
                if ~ismember(id, vertical_coil_ids)
                    wire_disp_1 = R*[-0.002-0.020;0;0];
                    wire_disp_2 = R*[-0.002-0.015;0;0];
                    wire_disp_3 = R*[-0.002-0.010;0;0];
                else
                    wire_disp_1 = R*[-0.002;0;0];
                end
                
                % Rotate the mesh for the core according to the coil
                % pose
                for j = 1:size(X_core,1)
                    for k = 1:size(X_core,2)
                        pos = [X_core(j,k); Y_core(j,k); Z_core(j,k)];
                        pos = R*roty(-90)*pos + core_disp;
                        X_core(j,k) = pos(1);
                        Y_core(j,k) = pos(2);
                        Z_core(j,k) = pos(3);
                    end
                end
                % Rotate the mesh for the wires according to the coil pose
                for j = 1:size(X_wire1,1)
                    for k = 1:size(X_wire1,2)
                        pos = [X_wire1(j,k); Y_wire1(j,k); Z_wire1(j,k)];
                        pos = R*roty(-90)*pos + wire_disp_1;
                        X_wire1(j,k) = pos(1);
                        Y_wire1(j,k) = pos(2);
                        Z_wire1(j,k) = pos(3);
                    end
                end
                
                % If coil is iron-cored, rotate mesh for the wires
                % according to the coil pose
                if ~ismember(id, vertical_coil_ids)
                    % Wires 2
                    for j = 1:size(X_wire2,1)
                        for k = 1:size(X_wire2,2)
                            pos = [X_wire2(j,k); Y_wire2(j,k); Z_wire2(j,k)];
                            pos = R*roty(-90)*pos + wire_disp_2;
                            X_wire2(j,k) = pos(1);
                            Y_wire2(j,k) = pos(2);
                            Z_wire2(j,k) = pos(3);
                        end
                    end
                    
                    % Wires 3
                    for j = 1:size(X_wire3,1)
                        for k = 1:size(X_wire3,2)
                            pos = [X_wire3(j,k); Y_wire3(j,k); Z_wire3(j,k)];
                            pos = R*roty(-90)*pos + wire_disp_3;
                            X_wire3(j,k) = pos(1);
                            Y_wire3(j,k) = pos(2);
                            Z_wire3(j,k) = pos(3);
                        end
                    end
                    X_wire2 = X_wire2+p(1);   Y_wire2 = Y_wire2+p(2);   Z_wire2 = Z_wire2+p(3);
                    X_wire3 = X_wire3+p(1);   Y_wire3 = Y_wire3+p(2);   Z_wire3 = Z_wire3+p(3);
                end
                
                % Displace meshes according to coil pose
                X_core = X_core+p(1);   Y_core = Y_core+p(2);   Z_core = Z_core+p(3);
                X_wire1 = X_wire1+p(1);   Y_wire1 = Y_wire1+p(2);   Z_wire1 = Z_wire1+p(3);
                
                
                
                
                % Draw mesh of the first set of wires
                mesh(X_wire1,Y_wire1,Z_wire1,'FaceAlpha',1,'EdgeColor','none', 'FaceColor',[0.72, 0.45, 0.20]); hold on;
                % Fill edges of the mesh
                fill3(X_wire1(1,:),Y_wire1(1,:),Z_wire1(1,:),[0.72, 0.45, 0.20]);
                fill3(X_wire1(2,:),Y_wire1(2,:),Z_wire1(2,:),[0.72, 0.45, 0.20]);
                
                % If coil is iron-cored, repeat for the core and
                % second/third set of wires
                if ~ismember(id, vertical_coil_ids)
                    mesh(X_core,Y_core,Z_core,'FaceAlpha',1,'EdgeColor','none', 'FaceColor',[0.75, 0.8, 0.8]);
                    fill3(X_core(1,:),Y_core(1,:),Z_core(1,:),[0.75, 0.8, 0.8]);
                    fill3(X_core(2,:),Y_core(2,:),Z_core(2,:),[0.75, 0.8, 0.8]);
                    
                    mesh(X_wire2,Y_wire2,Z_wire2,'FaceAlpha',1,'EdgeColor','none', 'FaceColor',[0.72, 0.45, 0.20]);
                    fill3(X_wire2(1,:),Y_wire2(1,:),Z_wire2(1,:),[0.72, 0.45, 0.20]);
                    fill3(X_wire2(2,:),Y_wire2(2,:),Z_wire2(2,:),[0.72, 0.45, 0.20]);
                    
                    mesh(X_wire3,Y_wire3,Z_wire3,'FaceAlpha',1,'EdgeColor','none', 'FaceColor',[0.72, 0.45, 0.20]);
                    fill3(X_wire3(1,:),Y_wire3(1,:),Z_wire3(1,:),[0.72, 0.45, 0.20]);
                    fill3(X_wire3(2,:),Y_wire3(2,:),Z_wire3(2,:),[0.72, 0.45, 0.20]);
                end
            end
            
            if opaque
                p1 = plot3([p(1) x(1)],[p(2) x(2)],[p(3) x(3)],'r','LineWidth',2.5);
                p2 = plot3([p(1) y(1)],[p(2) y(2)],[p(3) y(3)],'g','LineWidth',2.5);
                p3 = plot3([p(1) z(1)],[p(2) z(2)],[p(3) z(3)],'b','LineWidth',2.5);
                
                p1.Color(4) = 0.35;
                p2.Color(4) = 0.35;
                p3.Color(4) = 0.35;
            else
                plot3([p(1) x(1)],[p(2) x(2)],[p(3) x(3)],'r','LineWidth',2.5);
                plot3([p(1) y(1)],[p(2) y(2)],[p(3) y(3)],'g','LineWidth',2.5);
                plot3([p(1) z(1)],[p(2) z(2)],[p(3) z(3)],'b','LineWidth',2.5);
            end
            
            % Annotate coil number
            if ~isempty(id)
                text(x(1),x(2),x(3),['\{' num2str(id) '\}'], 'FontSize', 20)
            end
            
        end
        
        
        % Give a quiver plot a colormap
        function [mags,cbar] = color_quiver_plot(obj, q, varargin)
            % q = handle to quiver
            % Varargin = desired colormap
            %// Compute the magnitude of the vectors
            mags = sqrt(sum(cat(2, q.UData(:), q.VData(:), ...
                reshape(q.WData, numel(q.UData), [])).^2, 2));
            
            if nargin > 1
                colormap(varargin{1})
            end
            
            %// Get the current colormap
            currentColormap = colormap(gca);
            
            %// Now determine the color to make each arrow using a colormap.
            % Divide the magnitudes over a set of bins equal to the gca colormap size
            [~, ~, ind] = histcounts(mags, size(currentColormap, 1));
            
            %// Now map this to a colormap to get RGB
            cmap = uint8(ind2rgb(ind(:), currentColormap) * 255);
            cmap(:,:,4) = 255;
            cmap = permute(repmat(cmap, [1 3 1]), [2 1 3]);
            
            %// We repeat each color 3 times (using 1:3 below) because each arrow has 3 vertices
            set(q.Head, ...
                'ColorBinding', 'interpolated', ...
                'ColorData', reshape(cmap(1:3,:,:), [], 4).');   %'
            
            %// We repeat each color 2 times (using 1:2 below) because each tail has 2 vertices
            set(q.Tail, ...
                'ColorBinding', 'interpolated', ...
                'ColorData', reshape(cmap(1:2,:,:), [], 4).');
            
            cbar = colorbar;
            cbar.Location = 'westoutside';
            caxis([floor(min(mags)), ceil(max(mags))])
            
            if nargin > 2
                ylabel(cbar,varargin{2});
            end
            
            % Remove z-axis
            if nargin > 3 && varargin{3}
                fig = gca;
                fig.ZAxis.Visible = 'off';
                fig.ZGrid = 'off';
                fig.Color = 'none';
            end
        end
        
    end
end


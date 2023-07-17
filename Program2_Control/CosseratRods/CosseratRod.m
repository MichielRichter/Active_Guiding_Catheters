classdef CosseratRod
    % This class contains all dimensions and material properties of a
    % cylindrical rod
    
    properties
        % Inherent properties
        r_i; % Inner radius
        r_o; % Outer radius
        a; % Cross-sectional area
        l; % Length
        E; % Elastic modulus
        G; % Shear modulus
        rho; % Density
        
        % Derived properties
        I; % Second moment of area
        Iz; % Polar moment of area;
        Kb; % Bending stiffness matrix
        Ks; % Shear stiffness matrix
        
        % Other properties
        p0 = [0;0;0]; % Base position
        q0 = [1;0;0;0]; % Base orientation quaternion
        theta0 = 0; % Base rotation angle
        
        % Variable properties
        v_star = @(s) [0;0;1];
        u_star = @(s) [0;0;0];
    end
    
    
    %% Public methods
    methods (Access = public)
        % Constructor
        function obj = CosseratRod(rod_radius_inner, rod_radius_outer, rod_length, ...
                material_modulus_elastic, material_modulus_shear, material_density)
            % Assign arguments as object properties
            obj.r_i = rod_radius_inner;
            obj.r_o = rod_radius_outer;
            obj.l = rod_length;
            obj.E = material_modulus_elastic;
            obj.G = material_modulus_shear;
            obj.rho = material_density;
            
            % Assign derived properties
            rod_cross_sectional_area = pi * (rod_radius_outer^2 - rod_radius_inner^2);
            obj.a = rod_cross_sectional_area;
            rod_second_moment_area_xx = pi * (rod_radius_outer^4 - rod_radius_inner^4) / 4;
            rod_second_moment_area_yy = pi * (rod_radius_outer^4 - rod_radius_inner^4) / 4;
            rod_polar_second_moment_area = rod_second_moment_area_xx + rod_second_moment_area_yy;
            obj.I = rod_second_moment_area_xx;
            obj.Iz = rod_polar_second_moment_area;
            obj.Kb = diag( [material_modulus_elastic*rod_second_moment_area_xx, material_modulus_elastic*rod_second_moment_area_yy, material_modulus_shear*rod_polar_second_moment_area] );
            obj.Ks = diag( [material_modulus_shear*rod_cross_sectional_area, material_modulus_shear*rod_cross_sectional_area, material_modulus_elastic*rod_cross_sectional_area] );
            
        end
        
        % Set predefined precurvature vector handle
        function obj = setPrecurvatureHandle(obj, handle)
            obj.u_star = handle;
        end
        
    end
    
    %% Private methods
    methods (Access = private)
        
        
    end
end

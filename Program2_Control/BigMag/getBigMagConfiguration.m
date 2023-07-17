function [BigMag_positions, BigMag_currents] = getBigMagConfiguration(BigMag_object, magnet_position, magnet_moment, torque_desired)

%{
    Find the positions and currents of BigMag that exerts the desired
    torque on the magnet and and force in the magnet direction. That is:
    - Loop over forces in direction of magnet axis, such that the BigMag
    currents are minimized to achieve desired torque
%}
%% Compute
angles0 = [0;0];        forces0 = [0;0;0];

anglesLB = [-180; -180];    anglesUB = [180; 180];
forcesLB = -0.005*ones(3,1); forcesUB = abs(forcesLB); % Allow force of some magnitude (0.005 = 5 mN)

x0 = [angles0; forces0];
LB = [anglesLB; forcesLB];
UB = [anglesUB; forcesUB];

x = fmincon(@fun, x0, [], [], [], [], LB, UB, @nonlcon);



%% Output
% Extract BigMag positions and assign BigMag state
BigMag_positions = x(1:2);
BigMag_object = BigMag_object.update(zeros(6,1), BigMag_positions);

% Obtain unit-current field and gradient maps of BigMag
BigMag_unit_field_map = BigMag_object.getUnitFieldMap(magnet_position); % Unit field map
BigMag_unit_gradient_map = BigMag_object.getUnitGradientMap(magnet_position); % Unit gradient map

% Construct matrix of unit-current torques and forces on dipole
matrix_linear = [skew(magnet_moment)*BigMag_unit_field_map;
    getM(magnet_moment)*BigMag_unit_gradient_map];

% Compute required currents through the coils of BigMag
force_desired = x(3:5);
BigMag_currents = pinv(matrix_linear) * [torque_desired; force_desired];


%% Functions

    function I = fun(x)
        % Extract values
        BigMag_positions = x(1:2);
        force_desired = x(3:end);
        
        % Update BigMag position (currents do not matter here)
        BigMag_object = BigMag_object.update(zeros(6,1), BigMag_positions);
        
        % Compute BigMag unit-current field and gradient maps
        BigMag_unit_field_map = BigMag_object.getUnitFieldMap(magnet_position); % Unit field map
        BigMag_unit_gradient_map = BigMag_object.getUnitGradientMap(magnet_position); % Unit gradient map
        
        % Construct matrix of unit-current field and gradient maps
        matrix_linear = [skew(magnet_moment)*BigMag_unit_field_map;
            getM(magnet_moment)*BigMag_unit_gradient_map];
        
        % Compute minimum required currents to achieve desired torque
        % (target) and force
        currents_required = pinv(matrix_linear) * [torque_desired; force_desired];
        
        % Output of objective function: norm of current
        I = norm(currents_required);
    end

    function [c, ceq] = nonlcon(x)
        % Extract values
        BigMag_positions = x(1:2);
        force_desired = x(3:end);
        
        % Update BigMag position (currents do not matter here)
        BigMag_object = BigMag_object.update(zeros(6,1), BigMag_positions);
        
        % Compute BigMag unit-current field and gradient maps
        BigMag_unit_field_map = BigMag_object.getUnitFieldMap(magnet_position); % Unit field map
        BigMag_unit_gradient_map = BigMag_object.getUnitGradientMap(magnet_position); % Unit gradient map
        
        % Construct matrix of unit-current field and gradient maps
        matrix_linear = [skew(magnet_moment)*BigMag_unit_field_map;
            getM(magnet_moment)*BigMag_unit_gradient_map];
        
        % Compute minimum required currents to achieve desired torque
        % (target) and force
        currents_required = pinv(matrix_linear) * [torque_desired; force_desired];
        
        % Output
        c = abs(currents_required) - 5; % minimize the maximum required current as much as possible towards 5A (8A is the current limit for each coil)
        ceq = [];
    end

    function S = skew(a)
        S = [0 -a(3) a(2);
            a(3) 0 -a(1);
            -a(2) a(1) 0];
    end

    function M = getM(u)
        M = [u(1) u(2) u(3) 0 0;
            0 u(1) 0 u(2) u(3);
            -u(3) 0 u(1) -u(3) u(2)];
    end

end


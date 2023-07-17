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



%% Test
% Extract and assign BigMag state
BigMag_positions = x(1:2);
force_desired = x(3:end);
BigMag_object = BigMag_object.update(zeros(6,1), BigMag_positions);

% Compute unit-current field and gradient maps
BigMag_unit_field_map = BigMag_object.getUnitFieldMap(magnet_position); % Unit field map
BigMag_unit_gradient_map = BigMag_object.getUnitGradientMap(magnet_position); % Unit gradient map

% Compute currents
matrix_linear = [skew(magnet_moment)*BigMag_unit_field_map;
    getM(magnet_moment)*BigMag_unit_gradient_map];
currents_required = pinv(matrix_linear) * [torque_desired; force_desired];
BigMag_object = BigMag_object.update(currents_required, BigMag_positions);
torqueforce = matrix_linear * currents_required;
BigMag_currents = currents_required;


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

%% Obsolete
%{
% Initialize forces magnitude [N]
forces_magnitude = 0;% : 0.0005 : 0.01;
forces_desired = forces_magnitude .* (force_direction_desired/norm(force_direction_desired));


% Reset BigMag configuration
BigMag_positions = [0; 0];  dtheta = 0.1; % [deg]
BigMag_currents = zeros(6,1);

% Loop over desired forces
iter = 1;   solutions = NaN(8, size(forces_desired,2));
for force_desired = forces_desired
    
    % Update BigMag
    BigMag_object = BigMag_object.update(BigMag_currents, BigMag_positions);
    
    % Get unit-current field and gradient maps
    BigMag_unit_field_map = BigMag_object.getUnitFieldMap(magnet_position);
    BigMag_unit_gradient_map = BigMag_object.getUnitGradientMap(magnet_position);
    
    % Compute currents
    matrix_linear = [skew(magnet_moment)*BigMag_unit_field_map;
        getM(magnet_moment)*BigMag_unit_gradient_map];
    currents_required = pinv(matrix_linear) * [torque_desired; force_desired];
    
    % Initialize error, equal to norm of currents
    E_minus = norm(currents_required);  E_plus = 1e5;   E_change = abs(E_minus - E_plus);
    %E_minus = max(abs(currents_required));  E_plus = 1e5;   E_change = abs(E_minus - E_plus);
    
    % Optimize BigMag positions with Levenberg-Marquardt convex algorithm
    lambda = 1e-2;
    while(E_change > 1e-6)
        
        % Initialize Jacobian
        J = NaN(1,2);
        
        % Fill Jacobian
        for i = 1 : size(BigMag_positions,1)
            BigMag_positions(i) = BigMag_positions(i) + dtheta; % Increment BigMag position
            BigMag_object = BigMag_object.update(BigMag_currents, BigMag_positions); % Update BigMag
            BigMag_unit_field_map = BigMag_object.getUnitFieldMap(magnet_position); % Unit field map
            BigMag_unit_gradient_map = BigMag_object.getUnitGradientMap(magnet_position); % Unit gradient map
            % Compute new required currents
            matrix_linear = [skew(magnet_moment)*BigMag_unit_field_map;
                getM(magnet_moment)*BigMag_unit_gradient_map];
            currents_required = pinv(matrix_linear) * [torque_desired; force_desired];
            J(i) = (norm(currents_required) - E_minus) / dtheta; % Fill Jacobian entry
            %J(i) = (max(abs(currents_required)) - E_minus) / dtheta; % Fill Jacobian entry
            BigMag_positions(i) = BigMag_positions(i) - dtheta; % Decrement BigMag position
        end
        
        % Compute position increment
        lhs = J' * J;
        rhs = -J' * E_minus;
        lambda = lambda / 1.1;
        lhs = lhs + lambda * eye(size(J,2));
        step_position = lhs \ rhs;
        
        % Increment position
        BigMag_positions = BigMag_positions + step_position;
        BigMag_object = BigMag_object.update(BigMag_currents, BigMag_positions);
        
        % New unit-current field and gradient maps
        BigMag_unit_field_map = BigMag_object.getUnitFieldMap(magnet_position); % Unit field map
        BigMag_unit_gradient_map = BigMag_object.getUnitGradientMap(magnet_position); % Unit gradient map
        
        % Compute currents
        matrix_linear = [skew(magnet_moment)*BigMag_unit_field_map;
            getM(magnet_moment)*BigMag_unit_gradient_map];
        currents_required = pinv(matrix_linear) * [torque_desired; force_desired];
        
        % New residual
        E_plus = norm(currents_required);
        %E_plus = max(abs(currents_required));
        
        % Reset BigMag position
        BigMag_positions = BigMag_positions - step_position;
        BigMag_object = BigMag_object.update(BigMag_currents, BigMag_positions);
        
        % If solution is an improvement
        if E_plus < E_minus
            % Update position and error values
            BigMag_positions = BigMag_positions + step_position;
            E_change = abs(E_minus - E_plus);
            E_minus = E_plus;
        else
            while E_plus > E_minus
                % Increment lambda
                lambda = 1.05*lambda;
                % Recompute shift vector
                lhs = lhs + lambda * eye(size(J,2));
                step_position = lhs \ rhs;
                % Update position
                BigMag_positions = BigMag_positions + step_position;
                BigMag_object = BigMag_object.update(BigMag_currents, BigMag_positions);
                % Compute new residual
                BigMag_unit_field_map = BigMag_object.getUnitFieldMap(magnet_position); % Unit field map
                BigMag_unit_gradient_map = BigMag_object.getUnitGradientMap(magnet_position); % Unit gradient map
                matrix_linear = [skew(magnet_moment)*BigMag_unit_field_map;
                    getM(magnet_moment)*BigMag_unit_gradient_map];
                currents_required = pinv(matrix_linear) * [torque_desired; force_desired];
                E_plus = norm(currents_required);
                %E_plus = max(abs(currents_required));
                % Reset BigMag position
                BigMag_positions = BigMag_positions - step_position;
                BigMag_object = BigMag_object.update(BigMag_currents, BigMag_positions);
            end
            % Update position and error values
            BigMag_positions = BigMag_positions + step_position;
            E_change = abs(E_minus - E_plus);
            E_minus = E_plus;
        end
    end
    % Complete
    BigMag_currents = currents_required;
    
    % Add to solution
    solutions(:, iter) = [BigMag_positions; BigMag_currents];
    
    % Increment
    iter = iter + 1;
end


%{
BigMag_positions = zeros(2,1);
BigMag_currents = zeros(6,1);

solutions_angles = [];
solutions_current = [];
solutions_forces = [];
for i = 1 : 3 : 360
    for j = 1 : 3 : 360
        BigMag_positions = [i; j];
        % Update BigMag
        BigMag_object = BigMag_object.update(BigMag_currents, BigMag_positions);
        
        % Get unit-current field and gradient maps
        BigMag_unit_field_map = BigMag_object.getUnitFieldMap(magnet_position);
        BigMag_unit_gradient_map = BigMag_object.getUnitGradientMap(magnet_position);
        
        % Compute currents
        matrix_linear = skew(magnet_moment)*BigMag_unit_field_map;
        currents_required = pinv(matrix_linear) * torque_desired;
        force = getM(magnet_moment)*BigMag_unit_gradient_map*currents_required;
        
        if (max(abs(currents_required)) <= 8)
            solutions_angles = [solutions_angles [i;j]];
            solutions_current = [solutions_current currents_required];
            solutions_forces = [solutions_forces norm(force)];
        end
    end
end

%% Find solution that minimizes forces
minimum_force_norm = min(solutions_forces);
minimum_force_idx = find(solutions_forces == minimum_force_norm);
BigMag_positions = solutions_angles(:, minimum_force_idx);
BigMag_currents = solutions_current(:, minimum_force_idx);
%}

%% Find the BigMag configuration with minimum currents
solutions_current = solutions(3:end, :);
solutions_current_norm = vecnorm(solutions_current, 2);
minimum_current_norm = min(solutions_current_norm);
minimum_current_idx = find(solutions_current_norm == minimum_current_norm);

% Assign output
BigMag_positions = solutions(1:2, minimum_current_idx);
BigMag_currents = solutions(3:end, minimum_current_idx);
%}
end


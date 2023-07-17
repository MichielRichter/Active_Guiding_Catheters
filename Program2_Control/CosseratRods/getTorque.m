function [tau_mag, Y_C1f_star, Y_C1f] = getTorque(C1f, C1s, C2, p0, q0, nL, mL, p_target)

    % Number of rod discretizations (arbitrary)
    number_of_discretizations = 60;
    
    % Optimization vectors for BVPs of the AGC and CTR
    xi_C1f = zeros(6,1);
    xi_C12s = zeros(7,1);
    
    % Make ctr for rubber-phase AGC
    C1f_model = ConcentricTubeRobot([C1f]);
    
    % Solve shape of rubber-pase AGC
    Y_C1f_star = C1f_model.solve([p0; q0; xi_C1f], [nL; mL], 1e-12, number_of_discretizations); 
    
    % Compute angular error 
    alpha = computeError(C1f, C1s, C2, p0, q0, nL, mL, p_target);
    
    % Change with previous error
    E_minus = norm(alpha);
    E_plus = 1e5;
    E_change = abs(E_minus-E_plus);
    
    % Loop
    tau_mag = mL;
    lambda = 1e-4;
    while (E_change > deg2rad(0.5))
        
        % Compute Jacobian matrix
        J = getJacobian(C1f, C1s, C2, p0, q0, nL, tau_mag, p_target, alpha);
        
        % Compute shift vector for torque
        lhs = J' * J;
        rhs = -J' * alpha;
        lambda = 0.1*lambda;
        lhs = lhs + lambda * eye(size(J,2));
        delta_tau_mag = lhs \ rhs;
        delta_tau_mag = setLimit(delta_tau_mag, 1e-3);
        
        % Increment torque
        tau_mag = tau_mag + [delta_tau_mag; 0];
        
        % Compute new error
        alpha = computeError(C1f, C1s, ...
            C2, p0, q0, nL, tau_mag, p_target);
        E_plus = norm(alpha);
        
        % Reduce torque
        tau_mag = tau_mag - [delta_tau_mag; 0];
        
        % If error has reduced, keep going, else, increase lambda
        if E_plus < E_minus
            % update torque
            tau_mag = tau_mag + [delta_tau_mag; 0];
            E_change = abs(E_minus - E_plus)
            E_minus = E_plus;
            
        else 
            loop = 1; % Loop counter
            while E_plus > E_minus
                % Increase lambda
                lambda = 1.5 * lambda;
                % Compute new shift vector for torque
                lhs = lhs + lambda * eye(size(J,2));
                delta_tau_mag = lhs \ rhs;
                delta_tau_mag = setLimit(delta_tau_mag, 5e-3); % Limit permissible shift vector magnitude
                % Update torque
                tau_mag = tau_mag + [delta_tau_mag; 0];
                % Compute new error
                alpha = computeError(C1f, C1s, ...
                    C2, p0, q0, nL, tau_mag, p_target);
                E_plus = norm(alpha);
                % Decrement torque
                tau_mag = tau_mag - [delta_tau_mag; 0];
                
                % Increase loop counter and break if it exceeds an arbitrary
                % threshold
                loop = loop + 1;
                if (loop > 25)
                   break; 
                end
            end
            % If an improvement is found, update torque and errors
            tau_mag = tau_mag + [delta_tau_mag; 0];
            E_change = abs(E_minus - E_plus)
            E_minus = E_plus;
        end
    end
    
    % Solve final shape of rubber sheath
    Y_C1f = C1f_model.solve([p0; q0; xi_C1f], [nL; tau_mag], 1e-12, number_of_discretizations);
    
    
    
    %% Functions
    
    % Compute a Jacobian matrix that relates increments in tip-torque on
    % the AGC to changes in the angular error between the CTR and target 
    function J = getJacobian(C1f, C1s, C2, p0, q0, nL, mL, p_target, error_old)
        
        % Initialize Jacobian matrix
        J = NaN( size(error_old,1), size(mL,1)-1 );
        
        % Define step in torque to try
        increment = 1e-4;
        tau_new = mL + increment * eye(size(mL,1));
        
        % Initialize cell arrays to hold instances of AGC Cosserat rod / CTR objects.
        array_C1f = cell(1,size(mL,1));
        array_C1s = cell(1,size(mL,1));
        array_C2 = cell(1,size(mL,1));
        array_C12s = cell(1,size(mL,1));
        array_handle_polynomials = cell(1,size(mL,1));
        array_handle_shapetopose = cell(1,size(mL,1));
        array_handle_error = cell(1,size(mL,1));
        % Fill cell arrays to hold instances of AGC Cosserat rod / CTR
        % objects, necessary for parallel computation of the Jacobian
        % matrix.
        for i = 1 : size(mL,1)
            array_C1f{i} = C1f; % CosseratRod
            array_C1s{i} = C1s; % CosseratRod
            array_C2{i} = C2; % CosseratRod
            array_C12s{i} = ConcentricTubeRobot([C1f]); % ConcentricTubeRobot
            array_handle_polynomials{i} = @getPolynomialFromShape; % Defined below
            array_handle_shapetopose{i} = @getPoseFromShape; % Defined below
            array_handle_error{i} = @getError; % Defined below
        end
        
        % Increment torque values in a parallel fashion and compute its
        % effect on the change in angular error
        length = C1f.l;
        parfor i = 1 : size(tau_new,2)-1 % parallel for-loop
            
            % Initialize current torque
            tau = tau_new(:,i);
            
            % Extract C1f ConcentricTubeRobot
            ctr_C1f = array_C12s{i};
            
            % Solve shape of rubber-phase AGC
            ctr_C1f_Y = ctr_C1f.solve([p0; q0; xi_C1f], [nL; tau], 1e-12, number_of_discretizations);
            
            % Extract function handle to fit polynomial to shape solution
            handle_polynomial = array_handle_polynomials{i};
            
            % Fit symbolic polynomial
            ctr_C1f_p_sym = handle_polynomial(ctr_C1f_Y, length, 4);
            
            % Fit Frenet-Serret frames to shape polynomial
            [frenetSerret_R, frenetSerret_u] = getCurvature(ctr_C1f_p_sym, 0);
            
            % Assign precurvature vector to glass sheath
            array_C1s{i} = array_C1s{i}.setPrecurvatureHandle(frenetSerret_u);
            
            % Create CTR object
            ctr_C12s = ConcentricTubeRobot([array_C1s{i}, array_C2{i}]);
            
            % CTR boundary conditions (assume the field is still on)
            y0_new = [p0; rotm2quat(frenetSerret_R(0))'; xi_C12s; 0];
            yL_new = [nL; tau(1:2); tau(3); 0];
            
            % Solve shape of CTR
            ctr_C12s_Y = ctr_C12s.solve(y0_new, yL_new, 1e-12, number_of_discretizations);
            
            % Extract function handle to obtain tip position and direction
            % vector
            handle_shapetopose = array_handle_shapetopose{i};
            
            % Compute tip position and direction vector of CTR
            [ctr_C12s_pL, ctr_C12s_pdotL] = handle_shapetopose(ctr_C12s_Y);
            
            % Extract function handle to compute angular error
            handle_error = array_handle_error{i};
            
            % Compute new angular error
            error_new = handle_error(ctr_C12s_pL, ctr_C12s_pdotL, p_target);
            
            % Fill jacobian
            J(:,i) = (error_new - error_old) / increment;
        end
    end


    % Function to compute the angular error between the CTR tip direction
    % vector and the tip-target vector
    function alpha = getError(p_12, pdot_12, p_tar)
        p_tar_12 = p_tar - p_12;
        alpha = acos(dot( p_tar_12/norm(p_tar_12), pdot_12/norm(pdot_12)));
    end

    % Function to compute the CTR tip direction vector and tip-target
    % vector, and to compute the angular error
    function alpha = computeError(C1f, C1s, C2, p0, q0, nL, mL, p_target)
        
        % Initialize model for rubber-phase AGC
        ctr_C1f = ConcentricTubeRobot([C1f]);
        
        % Solve the BVP and obtain the shape solution for the rubber-phase
        % AGC
        ctr_C1f_Y = ctr_C1f.solve([p0; q0; xi_C1f], [nL; mL], 1e-12, number_of_discretizations);
        
        % Fit centerline curve to the shape solution
        ctr_C1f_p_sym = getPolynomialFromShape(ctr_C1f_Y, C1f.l, 4);
        
        % Compute curvature vector of the centerline curve
        [frenetSerret_R, frenetSerret_u] = getCurvature(ctr_C1f_p_sym, 0);
        
        % Assign curvature vector as a new precurvature vector
        C1s = C1s.setPrecurvatureHandle(frenetSerret_u);
        
        % Create CTR object
        ctr_C12s = ConcentricTubeRobot([C1s, C2]);
        
        % CTR boundary conditions
        y0 = [p0; rotm2quat(frenetSerret_R(0))'; xi_C12s; 0];
        yL = [nL; mL(1:2); mL(3); 0];
        
        % Solve shape of CTR
        ctr_C12s_Y = ctr_C12s.solve(y0, yL, 1e-12, number_of_discretizations);
        
        % Compute tip position and direction vector of CTR
        [ctr_C12s_pL, ctr_C12s_pdotL] = getPoseFromShape(ctr_C12s_Y);

        % Compute the angular error
        alpha = getError(ctr_C12s_pL, ctr_C12s_pdotL, p_target);
        
        % Update guesses
        xi_C1f = ctr_C1f_Y(8:end,1);
        xi_C12s = ctr_C12s_Y(8:end-1,1);
    end

    % Function to compute the position and direction vector based on a
    % shape solution
    function [p, p_s] = getPoseFromShape(Y)
        p = Y(1:3,end);
        p_s = Y(1:3,end) - Y(1:3,end-1);
        p_s = p_s / norm(p_s);
    end

    % Function to compute a symbolic centerline curve from a shape solution
    function p_s_sym = getPolynomialFromShape(Y, length, order)
        
        % Length of discretized segments
        ds = length / (size(Y,2)-1);
        
        % Number of centerline points
        numpoints = length/ds;
        
        % Base position
        p_base = Y(1:3,1);
        
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
            vec(rows,1) = p_i - p_base;
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
        p_s_sym = sym(p_base) + sym(mat) * vec;
    end

    % Function to limit the norm of the input vector
    function a = setLimit(a, maxval)
        if max(abs(a)) > maxval
            a = maxval * (a/max(abs(a)));
        else
            a = a;
        end
    end
end


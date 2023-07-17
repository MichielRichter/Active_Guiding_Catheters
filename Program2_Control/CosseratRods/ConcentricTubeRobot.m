classdef ConcentricTubeRobot
    % This class contains a collection of CosseratRod objects. It assumes
    % that all rods have a length equal or smaller than the first tube.
    
    %%
    properties
        % Inherent properties
        ctr_rods;
        N_rods;
        D; % discretizations per rod
        bc_dist_nm;
        
        % Gravity direction in the local frame
        gravity = 9.81 * [1;0;0];
        
    end
    
    %%
    methods (Access = public)
        function obj = ConcentricTubeRobot(ctr_rods)
            %{
                INPUTS
                ctr_rods - array of CosseratRod object. Object 1 is
                considered the outer tube (longest length). All
                subsequent objects are ordered with decreasing length.
            %}
            obj.ctr_rods = ctr_rods;
        end
        
        % Solve the concentric tube robot shape
        % Input: y0 - Base state vector
        % bc_dist_nm - Distal boundary condition for internal force (n) and moment (m)
        % discr_per_rod - Number of discretizations from the tip of one rod to the next
        function [Y, guess_solved] = solve(obj, y0, bc_dist_nm, thresh_sol, discr_per_rod)
            %{
                Boundary condition requires a desired distal value for
                unknown state variables: (n, m_b_xy, m_b_1_z, ..., m_b_n_z)
            %}
            obj.bc_dist_nm = bc_dist_nm;
            
            % Number of rods in the CTR
            obj.N_rods = size(obj.ctr_rods, 2);
            
            % Number of discretizations per unique rod length
            obj.D = discr_per_rod;
            
            % Shoot a solution using provided initial state vector y0
            [Y, rods_final_cols] = obj.shoot(y0);
            
            % Get boundary values (n, mbxy, mb1z, ..., mbnz) where rods
            % terminate
            boundary_values = obj.get_boundary_values(Y, rods_final_cols);
            
            % Get boundary residuals
            boundary_residuals = obj.get_boundary_residuals(boundary_values, bc_dist_nm);
            
            % Find solution
            E_minus = norm(boundary_residuals);
            E_plus = 1e5;
            E_change = abs(E_minus - E_plus);
            lambda = 1e-4;
            while (E_change > thresh_sol)
                % Get convex optimization increment vector
                J = obj.get_jacobian(y0, boundary_residuals);
                lhs = J' * J;
                rhs = -J' * boundary_residuals;
                lambda = lambda / 2;
                lhs = lhs + lambda * eye(size(J,2));
                delta_beta = lhs \ rhs;
                
                % Update base state vector
                y0(8:8+4 + obj.N_rods*1, 1) = y0(8:8+4 + obj.N_rods*1, 1) + delta_beta;
                
                % Shoot new solution
                [Y, rods_final_cols] = obj.shoot(y0);
                
                % Get residuals
                boundary_values = obj.get_boundary_values(Y, rods_final_cols);
                boundary_residuals = obj.get_boundary_residuals(boundary_values, bc_dist_nm);
                E_plus = norm(boundary_residuals);
                
                % Reset base state vector
                y0(8:8+4 + obj.N_rods*1, 1) = y0(8:8+4 + obj.N_rods*1, 1) - delta_beta;
                
                
                % If solution is an improvement
                if E_plus < E_minus
                    % Update base state vector
                    y0(8:8+4 + obj.N_rods*1, 1) = y0(8:8+4 + obj.N_rods*1, 1) + delta_beta;
                    E_change = abs(E_minus - E_plus);
                    E_minus = E_plus;
                else % If solution is not an improvement
                    % Increase lambda until an improvement is found
                    while E_plus > E_minus
                        % Increment lambda
                        lambda = 4*lambda;
                        % Recompute shift vector
                        lhs = lhs + lambda * eye(size(J,2));
                        delta_beta = lhs \ rhs;
                        % Update base state vector
                        y0(8:8+4 + obj.N_rods*1, 1) = y0(8:8+4 + obj.N_rods*1, 1) + delta_beta;
                        % Shoot new solution
                        [Y, rods_final_cols] = obj.shoot(y0);
                        % Get new residual norm
                        boundary_values = obj.get_boundary_values(Y, rods_final_cols);
                        boundary_residuals = obj.get_boundary_residuals(boundary_values, bc_dist_nm);
                        E_plus = norm(boundary_residuals);
                        % Reset base state vector
                        y0(8:8+4 + obj.N_rods*1, 1) = y0(8:8+4 + obj.N_rods*1, 1) - delta_beta;
                    end
                    y0(8:8+4 + obj.N_rods*1, 1) = y0(8:8+4 + obj.N_rods*1, 1) + delta_beta;
                    E_change = abs(E_minus - E_plus);
                    E_minus = E_plus;
                end
            end
            guess_solved = y0(8:8+4 + obj.N_rods*1, 1);
        end
    end
    
    %%
    methods (Access = private)

        % Forward integration with Runge-Kutta fourth order method
        function [Y, rods_final_cols] = shoot(obj, y0)

            % Get the Cosserat rod objects in the Concentric tube
            rods = obj.ctr_rods;
            
            % Initialize solution matrix & find unique rod lengths
            rod_numbers_unique = [1]; % First rod (outer tube) always has unique length
            Y_size_cols = 1+obj.D; 
            rod_prev_length = rods(1).l;
            for rod_number = 2:obj.N_rods
                rod_length = rods(rod_number).l;
                if (rod_length < rod_prev_length)
                    Y_size_cols = Y_size_cols + obj.D;
                    rod_numbers_unique = [rod_numbers_unique, rod_number];
                end
                rod_prev_length = rod_length;
            end
            Y = NaN(size(y0,1), Y_size_cols);
            Y(:,1) = y0;
            
            % Set terminating colums for each rod in the concentric tube
            rods_final_cols = NaN(1, obj.N_rods);
            last_final_col = Y_size_cols;
            for i = 1:obj.N_rods % Loop over the rods
                % Check if the rod is unique
                rod_unique = ismember(i, rod_numbers_unique);
                
                % If rod is unique add new columns
                if rod_unique
                    rods_final_cols(i) = last_final_col;
                    last_final_col = last_final_col - obj.D;
                else % Otherwise the column does not change
                    rods_final_cols(i) = last_final_col + obj.D;
                end
            end
            
            % Shoot a solution
            s = 0; 
            subsegment = 1;
            for rod_number = size(rod_numbers_unique,2):-1:1 % Discretize shortest rod length first 
                rod = rods(rod_number);
                ds = (rod.l - s) / obj.D; % Adjust step size
                
                % Increment along rod length
                for si = 0:ds:rod.l-s-ds
                    % Solve next state
                    y1 = Y(:, subsegment);
                    k1 = ds * obj.ODE(rods, y1, s);
                    y2 = y1 + 0.5*k1;
                    k2 = ds * obj.ODE(rods, y2, s);
                    y3 = y1 + 0.5*k2;
                    k3 = ds * obj.ODE(rods, y3, s);
                    y4 = y1 + k3;
                    k4 = ds * obj.ODE(rods, y4, s);
                    Y(:,subsegment+1) = y1 + (k1+2*k2+2*k3+k4)/6;
                    
                    % Increment along subsegment and arc length
                    subsegment = subsegment + 1;
                    s = s + ds;
                end
                %rods = rods(1:end-1); % Remove accountable rods
            end
        end
        
        
        % Arc-length derivative of state vector (y)
        function dot_y = ODE(obj, rods, yi, s)
            % Extract values from base state vector
            p = yi(1:3, 1);
            q = yi(4:7, 1);
            n = yi(8:10, 1);
            mbxy = yi(11:12, 1);
            mbiz = yi(13:13+obj.N_rods-1, 1);
            
            
            thetaiz = yi(13+obj.N_rods:end, 1);
            thetaiz = [0; thetaiz]; % Add theta for first rod
            
            % Initialize vectors for arc length derivatives
            dot_theta_i = zeros(1,size(thetaiz,1));
            dot_mbi_z = zeros(1,size(rods,2));
            
            % Derived values
            Rq = quat2rotm(q');
            Q = [-q(2:end)'; q(1)*eye(3)-obj.skew(q(2:end))];
            
            % Strain vectors
            v = [0;0;1];
            u1xy_numerator = mbxy;
            u1xy_denominator = 0;
            uiz = zeros(1, size(rods,2));
            for i = 1:size(rods,2)
                % Get the cosserat rod
                rod = rods(i);
                
                % Orientation w.r.t rod 1
                Rz = rotz( rad2deg(thetaiz(i)) );
                
                % Precurvature vector of rod
                u_star = rod.u_star(s);
                uxy_star = u_star(1:2);
                uz_star = u_star(3);
                
                
                % Rod material properties
                E = rod.E;
                I = rod.I;
                G = rod.G;
                J = rod.Iz;
                
                % Account for rods where s exceed their length
                if s >= rod.l
                    uxy_star = zeros(2,1);
                    uz_star = 0;
                    E = 0;
                    I = 0;
                    G = 0;
                    J = 0;
                end
                
                % Curvature vectors of rod 1 about local xy-axes
                u1xy_numerator = u1xy_numerator + E*I*Rz(1:2,1:2)*uxy_star;
                u1xy_denominator = u1xy_denominator + E*I;
                
                
                % Curvature vectors of all rods about local z-axis
                uiz_numerator = mbiz(i) + G*J*uz_star;
                uiz_denominator = G*J;
                uiz(i) = uiz_numerator / uiz_denominator;
                
                % Account for rods where s exceed their length
                if s >= rod.l
                    uiz(i) = 0;
                end
            end
            
            u1xy = u1xy_numerator ./ u1xy_denominator;
            u1 = [u1xy; uiz(1)];
            
            % Arc length derivatives
            ui = zeros(3, size(rods,2));
            mbi = zeros(3, size(rods,2));
            mb = 0;
            fi = zeros(3, size(rods,2));
            
            for i = 1:size(rods,2)
                % Get the cosserat rod
                rod = rods(i);
                
                % Orientation w.r.t rod 1
                Rz = rotz( rad2deg(thetaiz(i)) );
                
                if i == 1
                    v = rod.v_star(s);
                    dot_p1 = Rq * v;
                    dot_q1 = 0.5 * Q * Rq * u1;
                end
                
                % Precurvature vector of rod
                u_star = rod.u_star(s);
                
                % Arc length derivatives for each rod
                if s <= rod.l
                    dot_theta_i(i) = uiz(i) - uiz(1);
                    ui(:,i) = Rz*u1 + dot_theta_i(i)*[0;0;1];
                    mbi(:,i) = rod.Kb * (ui(:,i) - u_star);
                    mb = mb + Rz * mbi(:,i);
                    dot_mbi_z(i) = -[0,0,1] * obj.skew(ui(:,i)) * mbi(:,i);
                    fi(:,i) = rod.rho * rod.a * obj.gravity;
                    
                    zeros(3,1); % For now all external forces are zero
                else
                    dot_theta_i(i) = 0;
                    ui(:,i) = zeros(3,1);
                    mbi(:,i) = zeros(3,1);
                    dot_mbi_z(i) = 0;
                    fi(:,i) = zeros(3,1);
                end
                
            end
            
            f = sum(fi,2);
            dot_n = -f;
            dot_mbxy = -obj.skew(u1)*mb - obj.skew([0;0;1])*Rq'*n;
            dot_mbxy = dot_mbxy(1:2);
            
            % At this point we have: dot_p1, dot_q1, dot_n, dot_mbxy,
            % dot_mbi_z (i=1,2,...,n), dot_theta_i (i=1,2,...,n)
            % Make arc length derivative state vector
            dot_y = [dot_p1; dot_q1; dot_n; dot_mbxy; dot_mbi_z'; dot_theta_i(2:end)'];
            
        end
        
        
        % Jacobian matrix
        function J = get_jacobian(obj, y0, residuals_old)
            % Initialize Jacobian
            guess = y0(8:8+4 + obj.N_rods*1, 1);
            J = NaN( length(residuals_old), length(guess) );
            
            % Fill Jacobian
            beta = guess;
            increment = 1e-8;
            for i = 1:length(guess)
                % Increment base state vector
                beta(i) = beta(i) + increment;
                y0(8:8+4 + obj.N_rods*1, 1) = beta;
                
                % Shoot new solution
                [Y, rods_final_cols] = obj.shoot(y0);
                
                % Get boundary residuals
                boundary_values = obj.get_boundary_values(Y, rods_final_cols);
                residuals_new = obj.get_boundary_residuals(boundary_values, obj.bc_dist_nm);
                
                % Fill Jacobian
                J(:,i) = (residuals_new - residuals_old) / increment;
                
                % Decrement
                beta(i) = beta(i) - increment;
            end
        end
        
        
        % Boundary value residuals
        function residuals = get_boundary_residuals(obj, boundary_values, boundary_conditions)
            % A residual for each tube
            residuals = boundary_conditions - boundary_values;
        end
        
        % Boundary values
        function boundary_values = get_boundary_values(obj, Y, rods_final_cols)
            boundary_values = NaN(5 + obj.N_rods*1,1);
            
            % First assign n mbxy
            boundary_values(1:5) = Y(8:12, end);
            
            % Loop over remaining rods for mbiz i=1,2,...
            for i = 1:obj.N_rods
                boundary_values(5+i) = Y(12+i, rods_final_cols(i));
            end
        end
        
        
        % Skew-symmetric matrix
        function skew_symmetric_matrix = skew(obj, a)
            skew_symmetric_matrix = [  0   -a(3)  a(2) ;
                a(3)   0   -a(1) ;
                -a(2)  a(1)   0  ];
        end
        
    end
end


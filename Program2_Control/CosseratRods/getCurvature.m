function [frenetSerret_R, frenetSerret_u] = getCurvature(p_s_sym, writeToFile)

% Get 0-3rd derivative of shape polynomial
r = p_s_sym;
r_prime = diff(r, sym('s'));
r_prime_prime = diff(r_prime, sym('s'));
r_prime_prime_prime = diff(r_prime_prime, sym('s'));

% Compute symbolic axes for frenet-serret frames
t = r_prime/norm(r_prime);
r = cross(cross(r_prime,r_prime_prime)/norm(cross(r_prime,r_prime_prime)),t);
b = cross(r_prime,r_prime_prime)/norm(cross(r_prime,r_prime_prime));

% Compute symbolic bending and twist
kappa = norm(cross(r_prime,r_prime_prime))/(norm(r_prime)^3);
tau = dot(cross(r_prime,r_prime_prime),r_prime_prime_prime)/(norm(cross(r_prime,r_prime_prime))^2);

R_rbt = [r b t]; % Local coordinate frames (NT is the bending plane, B the axis of rotation)

% Make arc-length derivatives of Frenet-Serret frames
dts = kappa*r;
drs = -kappa*t + tau*b;
dbs = -tau*r;
R_rbt_prime = [drs dbs dts];

% Compute pre-curvature based on Frenet-Serret frames
u_skew = R_rbt'*R_rbt_prime;
u = [u_skew(3,2); u_skew(1,3); u_skew(2,1)];

% Output and write
frenetSerret_R = matlabFunction(R_rbt);
frenetSerret_u = matlabFunction([u_skew(3,2); u_skew(1,3); u_skew(2,1)]);

if (writeToFile)    matlabFunction(u,'File','precurvature_function');  end


end
clear;
close all;
%% 1) Define numerical variables

% link lengths (m)
l1 = 0.5;
l2 = 0.5;
l3 = 0.3;

%masses (Kg)
m1 = 0.05;
m2 = 0.5;
m3 = 0.3;
m4 = 0.5;
m5 = 0.05;
m6 = 0.5;

% gravity
g = 9.81; %m/s^2

% degrees of freedom
n = 5;

%% 2) Define real symbolic variables
syms t q1 q2 q3 q4 q5 q1dot q2dot q3dot q4dot q5dot real
syms x1 x2 x1dot x2dot real
syms tau1 tau2 tau3 tau4 real

q = [q1; q2; q3; q4; q5];
qdot = [q1dot; q2dot; q3dot; q4dot; q5dot];
tau = [tau1; tau2; tau3; tau4];
x = [x1; x2];
xdot = [x1dot; x2dot];
qbar = [q; x];
qbardot = [qdot; xdot];

%% 3) Find total kinetic energy for unpinned robot
%%% a) positions of six masses [ri = ....]
% Note - there should also be a handwritten submission for this derivation
% Note - assumes stance foot is at x 

% m1
r1 = x;
% r2
r2 = r1 + l1 * [cos(q1); sin(q1)];
% r3
r3 = r2 + l2 * [cos(q1+q2); sin(q1+q2)];
% r6
r6 = r3 + l3 * [cos(q1+q2+q5); sin(q1+q2+q5)];

% r4
r4 = r3 + l2 * [cos(q1+q2+q3); sin(q1+q2+q3)];
% r5
r5 = r4 + l1 * [cos(q1+q2+q3+q4); sin(q1+q2+q3+q4)];

%%% b)  Symbolic expressions for r and rdot

r1dot = jacobian(r1, qbar)*qbardot;
r2dot = jacobian(r2, qbar)*qbardot;
r3dot = jacobian(r3, qbar)*qbardot;
r4dot = jacobian(r4, qbar)*qbardot;
r5dot = jacobian(r5, qbar)*qbardot;
r6dot = jacobian(r6, qbar)*qbardot;
%% 

%%% c) Kinetic energy
K = 0.5 * (m1 * (r1dot'*r1dot) + m2 * (r2dot'*r2dot) ...
    + m3 * (r3dot'*r3dot) + m4 * (r4dot'*r4dot) ...
    + m5 * (r5dot'*r5dot) + m6 * (r6dot'*r6dot));
K = simplify(K);
%%
%%% d) Dbar
% Dbar=simplify(hessian(K,qbardot));
Dbar = hessian(K,qbardot);


%% 
%%% e) D matrix 
% This matrix D(q) is the kinetic energy of the pinned robot
% Submatrix of Dbar; first 5 rows and columns

D = Dbar(1:5, 1:5); 

%% 4
% Write the total potential energy of the pinned robot in terms of q and create a variable P containing
% it. Create a column vector G containing the symbolic gradient (i.e., the transposed Jacobian) of P,
h1 = 0;
h2 = l1*sin(q1);
h3 = h2 + l2*sin(q1+q2);
h4 = h3 + l2*sin(q1+q2+q3);
h5 = h4 + l1*sin(q1+q2+q3+q4);
h6 = h3 + l3*sin(q1+q2+q5);

h_mat = [h1,h2, h3, h4, h5, h6];
M_mat = [m1, m2, m3, m4, m5, m6];

P = g*dot(M_mat, h_mat);
G = jacobian(P,q)';

%% 5
% Define the input matrix B of the robot, of dimension 5 × 4. Note that the first configuration variable,
% q1 in unactuated, while q2, . . . , q5 are directly actuated.

B = [zeros(1, 4); eye(4,4);];
%% 6
% Using symbolic differentiation and the formula given in class for the Christoffel coefficients of D
% (note: D, not Dbar), find the Coriolis matrix C(q, ˙q) and name it C.

C = sym(zeros(n,n));
for i=1:n
    for j=1:n
        for k=1:n
            C(k,j)=C(k,j)+1/2*(diff(D(k,j),q(i))+diff(D(k,i),q(j))-diff(D(i,j),q(k)))*qdot(i);
        end
    end
end

%% 7 
% The impact map ∆ requires Dbar(q), which you’ve already computed, and the matrix-valued function
% E(q) = [(dpx)q I_2]
% where (dpx)q is the Jacobian of px(q), the vector with head at the swing foot and tail at x (see
% Figure 1). Define symbolic variables px and E containing px(q) and the expression above for E(q).

px = r5 - r1;
E=[jacobian(px,q), eye(2)];


%% 8
% Turn the symbolic expressions in D,Dbar,px,E,C,G into Matlab functions with appropriate inputs.

Dfun=matlabFunction(D,'Vars',{q});
Dbarfun = matlabFunction(Dbar, 'Vars', {qbar});
pxfun = matlabFunction(px, 'Vars', {q});
Efun = matlabFunction(E, 'Vars', {q});
Cfun = matlabFunction(C, 'Vars',{q, qdot});
Gfun = matlabFunction(G, 'Vars',{q});

Dinv = inv(D);
Dinvfun = matlabFunction(Dinv,'Vars',{q});
%% 9)
% Create a structure array named data containing these objects: Dfun,Dbarfun,Efun,Cfun,Gfun,B.
% For instance, the command data.D=Dfun will create a field D in the structure containing the function
% Dfun. Later on, you will pass this structure array to the ode45 function for numerical integration,
% and to the impact_map function for computation of the impact map.

data = struct;
data.D = Dfun;
data.Dbar = Dbarfun;
data.E = Efun;
data.C = Cfun;
data.G = Gfun;
data.B = B;
data.Dinv = Dinvfun;

% My changes
data.l = [l1 l2 l2 l1 l3];

%% 1) Define numerical value

% Define numerical qref
% q2ref = pi/6;
% q3ref = pi + pi/6;
% q4ref = -pi/6; 
% q5ref = -pi/6;
% qref = [q2ref;q3ref;q4ref;q5ref];

% Define two control gains Kp,Kd for the controller (2)
% placing the roots of the polynomial λ2 + Kdλ + Kp at {−20, −20}
Kd = 40;
Kp = 400;
H = [zeros(4,1) eye(4)];

data.H = H;
data.Kp = Kp;
data.Kd = Kd;
% data.qref = qref;

%% VHC Design
% beta = 0.316637;
% v1 = -0.894373;
% v2 = 1.9;
% a = [   2.824955653589793   0.424637665380008  -7.354342992280046  30.184724649940087 -35.944227984560072  13.322482661520024]';


% a =  [2.333387716653704   0.462881344738286  -6.085601718089357  31.898187171469729 -35.794497497124738  11.135134254596011]';
% beta = 0.808159502255349;
% v1 = -3.515306672443490;
% v2 = 3.883565855889314;


% Improved gait
% a =   [2.762805677485372   0.450224259171958  -7.160116751212509  30.175234835431542 -35.794643882019443  13.086875017408621]';
% beta = 0.378786915911426;
% v1 = -1.088503639682148;
% v2 = 2.113860411893566;
% 
% beta = 0.394039670996995;
% v1 = -1.131459980189369;
% v2 = 2.170720312833271;

% beta = 0.618395074187240;
% v1 = -2.049452736343263;
% v2 = 2.579343536371052;

% Pass VHC and V and M stability
beta = 0.394039670996995;
v1 = -1.131459980189369;
v2 = 2.170720312833271;

% a =   [2.974396028703572   0.293445937914470  -7.596127085702170  30.032522282572153 -36.311084631206221  13.915401249644194]';
% beta = 0.167166054446172;
% v1 = -0.467693955986393;
% v2 = 1.414876334896876;

% beta = 0.754803159389710;
% v1 = -3.031022546502935;
% v2 = 3.188188537607566;

%% Important values for sigma definition
qv_minus = [0.5*(pi - beta); pi + beta];
qv_plus = [0.5*(pi + beta); pi - beta];
% I = jacobian(Delta(3:4), qdot);
% I = double(subs(I, [q1, q2], q_minus'));
qv_tilde = qv_plus - qv_minus;
qv_tilde_1 = qv_tilde(1);


%% Define Sigma
syms theta

% knee angle
% alpha_0 = pi - pi/18;
alpha_max = pi;
alpha_delta = pi/9;
alpha = alpha_max - alpha_delta;

knee_contract = pi/8;

gamma = (pi-alpha )/2;

% Save to data structure
data.alpha = alpha;
data.gamma = gamma;
% Polynomial
% % k = 6;
% a_flipped = flip(a);
% phi_sym = poly2sym(a_flipped, theta);

%%
% Solve for polynomials
%% Pre and post impact

% target positions for biped
q5_desired = -pi/12;

% qv_minus = [0.5*(pi - beta); pi + beta];
% qv_plus = [0.5*(pi + beta); pi - beta];

% q1_plus = qv_plus(1) - gamma;
% q2_plus = pi - alpha;
% q3_plus = qv_plus(2);
% q4_plus = pi + alpha;
% 
% q1_minus = qv_minus(1) - gamma;
% q2_minus = pi - alpha;
% q3_minus = qv_minus(2);
% q4_minus = pi + alpha;
% q5_minus = q5_desired;
% 
% q_minus = [q1_minus; q2_minus; q3_minus; q4_minus; q5_minus];
% q_plus = impact_map([q_minus; zeros(6,1)], data);
% q_plus = q_plus(1:5); 
% 
% q_tilde = q_plus - q_minus;

%% P3
% 
% 
% b_vector = [q_plus(2); q_minus(2); pi; f_v1; v1; v2];
% a = phi_matrix \ b_vector;
% 

%% P4
% A4 = [1 0 0 0 0 0;
%      1 1 1 1 1 1;
%      1 0.5 0.5^2 0.5^3 0.5^4 0.5^5];
% %
% 
% q4_plus = q_plus(4);
% q4_minus = q_minus(4);
% q4_midpoint = pi + alpha - knee_contract;
% 
% b4 = [q4_plus; q4_minus; q4_midpoint];
% q4_a_coeffs = linsolve(A4, b4);
% phi4_sym = poly2sym(flip(q4_a_coeffs), theta);

%% P5
% A5 = [1 0 0 0 0 0;
%      1 1 1 1 1 1;
%      1 0.5 0.5^2 0.5^3 0.5^4 0.5^5];
% q5_plus = q_plus(5);
% q5_midpoint = q5_minus + beta/2;
% b5 = [q5_plus; q5_minus; q5_midpoint];
% q5_a_coeffs = linsolve(A5,b5);
% phi5_sym = poly2sym(flip(q5_a_coeffs), theta);
%%

% q1_plus = qv_plus(1) - gamma;
% q2_plus = pi - alpha;
% q3_plus = qv_plus(2);
% q4_plus = pi + alpha;
% 
% q1_minus = qv_minus(1) - gamma;
% q2_minus = pi - alpha;
% q3_minus = qv_minus(2);
% q4_minus = pi + alpha;
% q5_minus = q5_desired;

q2ref= pi - alpha;
q4ref = pi + alpha;

q_minus = [(pi-beta)/2 - (pi-alpha)/2; 
            pi - alpha;
            pi + beta;
            pi + alpha;
            q5_desired+beta];
q_plus = [(beta+pi)/2 - (pi-alpha)/2; 
            pi - alpha;
            pi - beta;
            pi + alpha;
            q5_desired;];

% q_plus = impact_map([q_minus; zeros(6,1)], data);
% q_plus = q_plus(1:5); 
q_tilde = q_plus - q_minus;
q_tilde_1 = q_tilde(1);

%% Compute I
% q_minus = sigma_fun(1);
% q_plus = sigma_fun(0);

delta1_qdot = [eye(5) zeros(5,4)]*( [ Dbarfun(q_minus) -Efun(q_minus).' ; Efun(q_minus) zeros(2,2) ]\[Dbarfun(q_minus)*[eye(5); zeros(2,5)]; zeros(2,5)] ) * qdot;

delta2_R = [1  1  1  1 0;...
        0  0  0 -1 0; ...
        0  0 -1  0 0; ...
        0 -1  0  0 0; ...
        0  0 -1  0 1];
delta2_d = [-3; 2; 2; 2; 1].*pi;

q_final = delta2_R*q + delta2_d;
qdot_final = delta2_R*delta1_qdot;
Delta = [q_final; qdot_final];

% I = subs(qdot_final, qdot, q_minus);
I = jacobian(qdot_final, qdot);

% f_v1 = 0.424637665380007;

f_v1 = -(q_tilde_1 * (I(3,3)*v1 - I(3,1)*q_tilde_1))/(I(1,3)*v1 - ...
    I(1,1)*q_tilde_1);

% q1_minus = qv_minus(1) - gamma;
% q2_minus = pi - alpha;
% q3_minus = qv_minus(2);
% q4_minus = pi + alpha;
% q5_minus = q5_desired;

q2_mid= pi - alpha;
q4_mid = pi + alpha -knee_contract;
q5_mid = q5_desired + beta/2;

A_VHC2 = [1 0 0 0 0 0;
         1 1 1 1 1 1;
         1 0.5 0.5^2 0.5^3 0.5^4 0.5^5];
b_VHC2 = [q_plus(2); q_minus(2);q2_mid];

A_VHC3 = [1 0 0 0 0 0;
         1 1 1 1 1 1;
         1 0.5 0.5^2 0.5^3 0.5^4 0.5^5;
         0 1 0 0 0 0;
         0 1 2 3 4 5;
         0 1 2*0.5 3*0.5^2 4*0.5^3 5*0.5^4;];
b_VHC3 = [q_plus(3); q_minus(3); pi; f_v1; v1; v2];

A_VHC4 = [1 0 0 0 0 0;
         1 1 1 1 1 1;
         1 0.5 0.5^2 0.5^3 0.5^4 0.5^5];
b_VHC4 = [q_plus(4); q_minus(4); q4_mid];

A_VHC5 = [1 0 0 0 0 0;
         1 1 1 1 1 1;
         1 0.5 0.5^2 0.5^3 0.5^4 0.5^5];
b_VHC5 = [q_plus(5); q_minus(5); q5_mid];



% solve for a
a2 = linsolve(A_VHC2,b_VHC2);
a3 = double(linsolve(A_VHC3,b_VHC3));
a4 = linsolve(A_VHC4,b_VHC4);
a5 = linsolve(A_VHC5,b_VHC5);


%%


% % virtual q1 and q2 (from acrobot)
% qv1_sym = qv_plus(1)- theta*qv_tilde_1; 
% qv2_sym = phi_sym;
% 
% 
% sigma_q1 = poly2sym(a1_fl, theta);
% sigma_q2 = pi - alpha;
phi1 = @(theta) q_plus(1) - theta*q_tilde_1;

a2_flipped = flip(a2);
phi2=@(theta) polyval(a2_flipped,theta);
phiprime2=@(theta) polyval(polyder(a2_flipped),theta);
phipprime2=@(theta) polyval(polyder(polyder(a2_flipped)),theta);

a3_flipped = flip(a3);
phi3=@(theta) polyval(a3_flipped,theta);
phiprime3=@(theta) polyval(polyder(a3_flipped),theta);
phipprime3=@(theta) polyval(polyder(polyder(a3_flipped)),theta);

a4_flipped = flip(a4);
phi4=@(theta) polyval(a4_flipped,theta);
phiprime4=@(theta) polyval(polyder(a4_flipped),theta);
phipprime4=@(theta) polyval(polyder(polyder(a4_flipped)),theta);

a5_flipped = flip(a5);
phi5=@(theta) polyval(a5_flipped,theta);
phiprime5=@(theta) polyval(polyder(a5_flipped),theta);
phipprime5=@(theta) polyval(polyder(polyder(a5_flipped)),theta);



sigma_fun = @(theta) [phi1(theta); phi2(theta); phi3(theta); phi4(theta); phi5(theta)];

sigma_prime_fun = @(theta) [-q_tilde_1; phiprime2(theta); phiprime3(theta); phiprime4(theta); phiprime5(theta) ];
delta_numerator = sigma_prime_fun(0)' * I * sigma_prime_fun(1);
delta_denominator = sigma_prime_fun(0)' * sigma_prime_fun(0);
delta = delta_numerator / delta_denominator;

data.sigma_fun = sigma_fun;

data.phiprime2 = phiprime2;
data.phiprime3 = phiprime3;
data.phiprime4 = phiprime4;
data.phiprime5 = phiprime5;
%%
% a_flipped=flip(a);
% phi=@(theta) polyval(a_flipped,theta);
% phiprime=@(theta) polyval(polyder(a_flipped),theta);
% phipprime=@(theta) polyval(polyder(polyder(a_flipped)),theta);
% 
% data.phi = phi;
% data.phiprime= phiprime;
% data.phipprime = phipprime;
% 
% % virtual q1 and q2 (from acrobot)
% qv1_sym = qv_plus(1)- theta*qv_tilde_1; 
% qv2_sym = phi_sym;
% 
% 
% sigma_q1 = poly2sym(a1_fl, theta);
% sigma_q2 = pi - alpha;
% sigma_q3 = qv2_sym;
% sigma_q4 = phi4_sym;
% sigma_q5 = phi5_sym;

%%


phi1_sym =  q_plus(1) - theta*q_tilde_1;
phi2_sym = poly2sym(a2_flipped, theta);
phi3_sym = poly2sym(a3_flipped, theta);
phi4_sym = poly2sym(a4_flipped, theta);
phi5_sym = poly2sym(a5_flipped, theta);
% 
sigma_sym = [phi1_sym; phi2_sym; phi3_sym; phi4_sym; phi5_sym];
sigma_prime_sym = diff(sigma_sym, theta);
sigma_pprime_sym = diff(sigma_prime_sym, theta);
% 
% sigma_fun = matlabFunction(sigma_sym,'Vars',theta);
% sigma_fun_1 = matlabFunction(sigma_q1,'Vars',theta);
% sigma_fun_3 = matlabFunction(sigma_q3,'Vars',theta);
% sigma_fun_4 = matlabFunction(sigma_q4,'Vars',theta);
% sigma_fun_5 = matlabFunction(sigma_q5,'Vars',theta);

%%
% sigma_prime_sym = diff(sigma_sym, theta);
% sigma_pprime_sym = diff(sigma_prime_sym, theta);
% 
% sigma_prime_fun = matlabFunction(sigma_prime_sym,'Vars',theta);
% TODO - re-examine

% q_plus = sigma_fun(0);
% q_minus = sigma_fun(1);
% q_tilde = q_plus - q_minus;

% data.sigma_fun = sigma_fun;
% data.sigma_prime_fun = sigma_prime_fun;
data.q_tilde = q_tilde;
data.q_plus = q_plus;


data.qv_plus = qv_plus; 
data.qv_tilde = qv_tilde; 


%% Check conditions
% condition_11a = le(v1,2*q_tilde_1);
% if(not(condition_11a))
%     disp("Error: inequality 11a) is not satisfied")
%     % return;
% end
% 
% condition_11b = ge(v2, 2*q_tilde_1);
% if(not(condition_11b))
%     disp("Error: inequality 11b) is not satisfied")
%     % return;
% end

%% Symbolic error calculations
% e_sym = q - sigma_sym;

%% Compute the quantity B⊥ D(σ(θ))σ′(θ) (you’ll find that it is actually a constant, independent of θ), and verify that it is not zero.

B_perp = [1 0 0 0 0];
% sigma = [theta;qref];
% % verify_num = B_perp*Dfun(theta,q2ref,q3ref,q4ref,q5ref)*[diff(theta);0;0;0;0];
% verify_num = B_perp*Dfun([theta; qref])*[diff(theta);0;0;0;0];

B_perp_D_sigma_prime_sym = B_perp*Dfun(sigma_sym)*sigma_prime_sym;
B_perp_D_sigma_prime_fun=matlabFunction(B_perp_D_sigma_prime_sym,'Vars',{theta});

theta_vals =linspace(0,1,1001);
B_perp_D_sigma_prime_vals = B_perp_D_sigma_prime_fun(theta_vals);

figure()
plot(theta_vals', B_perp_D_sigma_prime_vals');
xlabel("\theta");
ylabel("B^\perp D \sigma'");

% check for zero crossings:
has_zero_crossings = min(B_perp_D_sigma_prime_vals(:)) < 0 && max(B_perp_D_sigma_prime_vals(:)) > 0;
if(has_zero_crossings)
    disp("Error: Curve is not a regular VHC, B^perp D sigma'(theta) has zero crossings")
    % return;
else
    disp("Passes check: Curve is a regular VHC B^perp D sigma'(theta) has no zero crossings")
end
%% Compute Psi1 Psi2
denominator=B_perp_D_sigma_prime_sym;
% Psi1=simplify(-subs((B_perp*G)/denominator,q,sigma_sym));
Psi1 = B_perp_D_sigma_prime_fun(sigma_sym);
%%
Psi2=-subs(B_perp*(D*sigma_pprime_sym+subs(C,[q;qdot],[sigma_sym;sigma_prime_sym])*sigma_prime_sym)/denominator,q,sigma_sym);


Psi1fun = matlabFunction(Psi1 ,'Vars',{theta});
Psi2fun = matlabFunction(Psi2 ,'Vars',{theta});

%% Find virtual mass M(theta) and virtual potential V(theta)

data2.Psi1=Psi1fun;
data2.Psi2=Psi2fun;

ops2=odeset('RelTol',1e-4,'AbsTol',1e-4);
[Theta,X]=ode45(@mass_potential,linspace(0,1,1000),[1;0;0;0;0; 0],ops2,data2);
M=spline(Theta,X(:,1));
V=spline(Theta,X(:,2));

figure()
plot(Theta,X(:,1));
hold on
plot(Theta,X(:,2));
legend('M','V')
title("Virtual Mass and Potential");
xlabel("\theta")
%% Check existence and stability conditions

% Computing parameters for inequalities
M_minus = ppval(M,1);
V_minus = ppval(V,1);
V_max = max(X(:,2));

delta_numerator = sigma_prime_fun(0)' * I * sigma_prime_fun(1);
delta_denominator = sigma_prime_fun(0)' * sigma_prime_fun(0);
delta = delta_numerator / delta_denominator;
% Checking inequalities:
% 12a) 
d_squared_over_M_minus = delta^2/M_minus;

existence = d_squared_over_M_minus > 0 && d_squared_over_M_minus < 1;
stability_value = (V_minus * delta^2)/(M_minus - delta^2)+V_max;
stability = (stability_value) < 0;

if existence && stability
    disp("Passes check: Virtual Mass / Potential stability and existence conditions are met");
else
    disp("Error: Virtual Mass / Potential stability and existence conditions fail");
end

% MY CODE END

%% Verify curve is contained in W for all theta (0,0.5) U (0.5,1)
theta_vals =linspace(0,1,1001);

q1_vals = phi1 (theta_vals);
q3_vals = phi3(theta_vals);

limit_line = -2*q1_vals+2*pi - 2*gamma; %+alpha/2;
figure()
plot(q1_vals, q3_vals);
hold on;
plot(q1_vals, limit_line);
hold off;

% yline(3*pi);
legend(["VHC", "-2q1+2pi"])
xlabel("q1");
ylabel("q3");
%%
% q_vals = [q1_vals; q4_vals];
% figure()
% title("q vs theta")
% plot(theta_vals, q_vals);
% xlabel("\theta");
% ylabel("q");
% legend(["q1", "q4"]);
%%
% 
% index_of_05=find(theta_vals == 0.5);
% 
% q1_vals = q_vals(1,:);
% q2_vals = q_vals(2,:);
% 
% 
% % check conditions of first half of curve
% interval1_indices = q1_vals < pi/2 & q1_vals > 0;
% q2_vals_interval1 = q2_vals(interval1_indices);
% 
% in_W1 = all(limit_line(interval1_indices) <= q2_vals_interval1) && all(q2_vals_interval1 < 3*pi);
% 
% if(in_W1)
%     disp("Passes check: VHC is contained in W1")
% else
%     disp("Error: VHC is not contained in W1")
%     % return;
% end
% 
% % check conditions of first half of curve
% interval2_indices = q1_vals > pi/2 & q1_vals < pi;
% q2_vals_interval2 = q2_vals(interval2_indices);
% 
% q2_vals_2 = q_vals(2,(index_of_05+1):end);
% limit2 = limit_line(:,(index_of_05+1):end);
% 
% in_W2 = all(-pi < q2_vals_interval2) && all(q2_vals_interval2 <= limit_line(interval2_indices));
% 
% if(in_W2)
%     disp("Passes check: VHC is contained in W2")
% else
%     disp("Error: VHC is not contained in W2")
%     % return;
% end

%% Check if constraints induce stable hybrid limit cycle
% denominator = B_perp*D


%% Section 5 New: Numerical Simulation

ops = odeset(RelTol=1e-10,AbsTol=1e-10,Events=@ground_impact);

% initial condition
q0=[pi/3;pi/4;pi-pi/6;-pi/3;-pi/6]; 
qdot0=zeros(5,1);

theta_0 = 0;
% q0 =  sigma_fun(theta_0);
% % % theta_dot_a = 3.928935686646192;
% % qdot0 = zeros(5,1);
% % sigma_prime_fun = @(theta) [-q_tilde_1; phiprime2(theta); phiprime3(theta); phiprime4(theta); phiprime5(theta) ];
% % qdot0(1) = -q_tilde_1;
% % qdot0(3) = phiprime3(theta_0);
% theta_dot_a = delta * sqrt(-2*V_minus/(M_minus-delta^2));
% qdot0 = sigma_prime_fun(theta_0) * 4;

dt = 1/60;
number_steps = 6;

T=[];
X=[];
Te=[];
Ie=[];
Xe=[];
post_impact_state=[q0;qdot0];
% Simulate number_steps steps
for step=1:number_steps
    fprintf('\n...step %d...\n',step);
    % [t,x,te,xe,ie]=ode45(@(t,x) acrobot(t,x,data),0:dt:10,post_impact_state,ops);
    [t,x,te,xe,ie]=ode45(@(t,x) biped(t,x,data),0:dt:10,post_impact_state,ops);    
    % Application of the impact map
    impact_state=xe';
    post_impact_state=impact_map(impact_state,data);
    T{step}=t;
    X{step}=x;
    Ie{step}=ie;
    Te{step}=te;
    Xe{step}=xe;
end

%% Animation of the simulation results


fprintf('\n Animation is ready...\n')
figure();

Axis=[-1 4 0 2];
axis equal;
step =0;
Time=text(-1+2,1.8,['time= ','0',' secs,',' step= ',num2str(step)]);
axis(Axis);
q=q0;

q1=q(1);q2=q(2);q3=q(3);q4=q(4);q5=q(5);
xdata=0;
ydata=0;
l=[l1 l2 l2 l1 l3];
Q=[q1 q1+q2 q1+q2+q3 q1+q2+q3+q4 q1+q2+q5];
for j=1:4
    xdata=[xdata xdata(end)+l(j)*cos(Q(j))];
    ydata=[ydata ydata(end)+l(j)*sin(Q(j))];
end
xdata=[xdata xdata(3)+l(5)*cos(Q(5))];
ydata=[ydata ydata(3)+l(5)*sin(Q(5))];

link1=line([xdata(1) xdata(2)],[ydata(1) ydata(2)],'color','red','linewidth',2);
link2=line([xdata(2) xdata(3)],[ydata(2) ydata(3)],'color','red','linewidth',2);
link3=line([xdata(3) xdata(4)],[ydata(3) ydata(4)],'linewidth',2);
link4=line([xdata(4) xdata(5)],[ydata(4) ydata(5)],'linewidth',2);
link5=line([xdata(3) xdata(6)],[ydata(3) ydata(6)],'linewidth',2);

% Other variables
ref=0; % This variable keeps track of the position of the stance foot accross multiple steps
time_passed=0;step=1;
%%
fprintf('\n Creating video file...\n')
v = VideoWriter("walking.avi");
open(v)

animation_slowdown_factor=10; % >1 means slow down
for step=1:length(Ie)
    t=T{step};
    x=X{step};
    xe=Xe{step};
    xe=xe(end,:);
    x_data_5 = 0;
    for k=2:length(t)
        t0=clock;
        drawnow;
        q=x(k,1:5)';
        q1=q(1);q2=q(2);q3=q(3);q4=q(4);q5=q(5);
        Q=[q1 q1+q2 q1+q2+q3 q1+q2+q3+q4 q1+q2+q5];
        xdata=0;
        ydata=0;
        for j=1:4
            xdata=[xdata xdata(end)+l(j)*cos(Q(j))];
            ydata=[ydata ydata(end)+l(j)*sin(Q(j))];
        end
        xdata=ref+[xdata xdata(3)+l(5)*cos(Q(5))];
        ydata=[ydata ydata(3)+l(5)*sin(Q(5))];
        
        x_data_5 = xdata(5);

        set(link1,'xdata',[xdata(1) xdata(2)],'ydata',[ydata(1) ydata(2)]);
        set(link2,'xdata',[xdata(2) xdata(3)],'ydata',[ydata(2) ydata(3)]);
        set(link3,'xdata',[xdata(3) xdata(4)],'ydata',[ydata(3) ydata(4)]);
        set(link4,'xdata',[xdata(4) xdata(5)],'ydata',[ydata(4) ydata(5)]);
        set(link5,'xdata',[xdata(3) xdata(6)],'ydata',[ydata(3) ydata(6)]);
        set(Time,'String',['time= ',num2str(round(time_passed+t(k),1)),' secs,',' step= ',num2str(step)]);
        current_axis=gca;
        if ref>.95*current_axis.XLim(end)
            current_axis.XLim=[.95*ref .95*ref+5];
            Time.Position=[.95*ref+2 1.8 0];
            Axis=axis;
        else
            axis(Axis)
        end
        while etime(clock,t0)<animation_slowdown_factor*(t(k)-t(k-1))
        end
        frame = getframe(gcf);
        writeVideo(v,frame)
    end
    time_passed=time_passed+t(end);
    % ref=ref+l*(cos(xe(1))+cos(xe(1)+xe(2)));
    ref=x_data_5;
end

close(v);
%%

function xdot=mass_potential(theta,x,data2)
M=x(1);
V=x(2);

Mdot = -2*M*data2.Psi2(theta);
Vdot = -data2.Psi1(theta)*M;
xdot=[Mdot;Vdot];
end
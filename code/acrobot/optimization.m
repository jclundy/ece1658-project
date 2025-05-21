%% Setup

%initialization
k_dim = 6;

beta_0 = 0.316637;
v1_0 = -0.894373;
v2_0 = 1.9;
a_0 = [   2.824955653589793   0.424637665380008  -7.354342992280046  30.184724649940087 -35.944227984560072  13.322482661520024]';

% v1_0 = -2;
% v2_0 = 2;
% a_0 = ones(k_dim,1);
% beta_0 = 0.5;
% X0 = X_optimized;
% 
% beta_test = 0.495242040516821;
% v1_test = -0.717699683565239;
% v2_test = 2.981043469410864;
% a_test = [2.740209691408301   3.310831198181539 -29.545021722053210  90.440054573900213 -99.999999999996092  36.743711443416068 ];
% X0 = [a_test, beta_test, v1_test, v2_test]';


% a_0 = [2.747552982544943   0.460240412763022  -7.174554868120053  30.384908371960194 -36.015257100612303  13.132742525684138]';
% % a_unflipped = a;
% beta_0 = 0.394039670996995;
% v1_0 = -1.131459980189369;
% v2_0 = 2.170720312833271;

% a_0 = [2.523197579294350   0.496968421420248  -3.015253857554625  16.912811540665107 -18.534331717251145   5.376595760336740]';
% beta_0 = 0.618395074187240;
% v1_0 = -2.049452736343263;
% v2_0 = 2.579343536371052;

% 
% a_0 =   [2.126418375652961   0.161492962044396  18.827649003614532 -66.783543202174783  96.589901601987734 -46.765153591876597 ]';  
% beta_0 = 1.015173838411510;
% v1_0 = -9.999999999998076;
% v2_0 = 2.582325256489500;

% a_0 = [   1.201316970125183 -12.853020374666150  68.622164548477954 -99.999999892254550  70.630399558391687 -22.520156184639660]';
% beta_0 = 1.937450684508994;
% v1_0 = -5.687864725175802;
% v2_0 = 9.046803716923238;

X0 = [a_0; beta_0; v1_0; v2_0];

%% Psi1 and Psi2

syms th;
k = 6;
a_var = sym("a", [k 1]);
Phi_var = sym("Phi", [1 k]);
for i=1:k
    Phi_var(i) = th^(i-1);
end

Phi_var = Phi_var * a_var;
Phi_prime_var = diff(Phi_var, th);
Phi_pprime_var = diff(Phi_prime_var, th);

Phi_func_a_theta = matlabFunction(Phi_var, 'Vars',{a_var, th});
Phiprime_func_a_theta = matlabFunction(Phi_prime_var, 'Vars',{a_var, th});
Phipprime_func_a_theta = matlabFunction(Phi_pprime_var, 'Vars',{a_var, th});

beta_fun = @(X,k)(X(k+1));
q_minus_fun = @(X,k)[0.5*(pi - beta_fun(X,k)); pi + beta_fun(X,k)];
q_plus_fun = @(X,k)[0.5*(pi + beta_fun(X,k)); pi - beta_fun(X,k)];

Sigma_func_a_theta_qp1_qt1 = @(a, th, qp1, qt1) [qp1- th*qt1; Phi_func_a_theta(a, th)];
Sigmaprime_func_a_theta_qp1_qt1 = @(a, th, qt1) [-qt1; Phiprime_func_a_theta(a, th)];
Sigmapprime_func_a_theta_qp1_qt1 = @(a, th) [0; Phipprime_func_a_theta(a, th)];

testData0.Sigma_fun = Sigma_func_a_theta_qp1_qt1;
testData0.Sigmaprime_fun = Sigmaprime_func_a_theta_qp1_qt1;
testData0.Sigmapprime_fun = Sigmapprime_func_a_theta_qp1_qt1;

% 
% Sigma_func_X_theta = @(X, th) [q_plus(1)- theta.*q_tilde_1; phi(theta)];
% sigma=@(theta) [q_plus(1)- theta.*q_tilde_1; phi(theta)];
% sigmaprime = @(theta) [-q_tilde_1; phiprime(theta)];
% sigmapprime= @(theta) [0; phipprime(theta)];

%% 
syms qp1 qt1

sigma_var = [qp1- th*qt1; Phi_var];
sigma_prime_var = diff(sigma_var, th);
sigma_pprime_var = diff(sigma_prime_var, th);

% Psi1
denominator=simplify(B_perp*D*sigma_prime_var);
Psi1_var=simplify(-B_perp*G/denominator);
Psi1_var= subs(Psi1_var, [q1; q2], sigma_var);
Psi1_fun = matlabFunction(Psi1_var, 'Vars', {a_var, th, qp1, qt1});
testData0.Psi1fun=Psi1_fun;

% Psi2
Psi2_var=-B_perp*(D*sigma_pprime_var+subs(C,qdot,sigma_prime_var)*sigma_prime_var)/denominator;    
Psi2_var=simplify(subs(Psi2_var,[q1; q2],sigma_var));
Psi2_fun = matlabFunction(Psi2_var, 'Vars',  {a_var, th, qp1, qt1});
testData0.Psi2fun  = Psi2_fun;

%%
% B_perp_D_sigma_prime 
B_perp_D_sigma_prime = subs(denominator, [q1; q2], sigma_var);
B_perp_D_sigma_prime_fun =matlabFunction(B_perp_D_sigma_prime, 'Vars', {a_var, th, qp1, qt1});

testData0.B_perp_D_sigma_prime_fun = B_perp_D_sigma_prime_fun;
%%
% I matrix
% syms q1 q2 q1dot q2dot;
Imat = jacobian(Delta(3:4), [q1dot; q2dot]);
Ifun = @(q_min)subs(Imat, [q1; q2],q_min);
testData0.Ifun = Ifun;

%%
% Other parameters
% try 
% 90 deg: pi/2
% 120 deg 2*pi/3
% 
% targetBeta = pi/12; %15 deg
% targetBeta = pi/3;
targetBeta = beta_0;

testData0.Delta = Delta;
testData0.C = C;
testData0.G = G;
testData0.B_perp = [1, 0];
testData0.D = D;
testData0.epsilon = 1e-2;
testData0.N = 1000;
testData0.beta = targetBeta;
testData0.V_max = 0;
testData0.stability_value = 0;
% f = @(X_opt)parameterfun(X_opt,testData0);

%%

options = optimoptions('fmincon', 'Display','iter');
options.MaxIterations = 1000;
% options.StepTolerance = 1e-7;
problem.options = options;
problem.solver='fmincon';
problem.objective = @(X_opt)objective_function(X_opt, testData0);
problem.x0 = X0;
problem.nonlcon = @(X_opt)generate_nonlinear_constraint(X_opt, testData0);
problem.ub = [100*ones(1,k_dim), 3*pi/4, 10, 10];
problem.lb = [-100*ones(1,k_dim), 0, -10, -10];

%%
% Test functions

objective_function(X0, testData0)
%%
transversality_condition(X0, testData0)

%%
% Test M,V
[a_coeffs, k] = get_a_k(X0);
[qt1, qp1] = get_qt1_qp1(X0,k);
solver_data.Psi1 = @(theta) testData0.Psi1fun(a_coeffs, theta, qp1, qt1);
solver_data.Psi2 = @(theta) testData0.Psi2fun(a_coeffs, theta, qp1, qt1);

[M_test,V_test, MV_vals] = compute_M_V(solver_data);

%%

% Transversality condition
% testData1 = testData0;
% testData1.epsilon = 1e-3;
c3 = transversality_condition(X0, testData0)

% c3 = transversality_condition(X_optimized, testData0)
%%

% Compute virtual Mass and Potential

% Extract parameters
% [qt1, qp1] = get_qt1_qp1(X,k);
% solver_data.Psi1 = @(theta) data.Psi1fun(a, theta, qp1, qt1);
% solver_data.Psi2 = @(theta) data.Psi2fun(a, theta, qp1, qt1);
% 
% [M,V, MV_vals] = compute_M_V(solver_data);

% Stability conditions from (12)
[c4, c5, c6] = existence_and_stability_conditions(X0, M_test, V_test, MV_vals, testData0)
%%

M_minus_test = ppval(M_test,1)
V_minus_test = ppval(V_test,1)
V_max_test = max(MV_vals(:,2))
% data.V_max = V_max_test;
%%
sigma_prime_0 = testData0.Sigmaprime_fun(a_coeffs, 0, q_tilde_1)
sigma_prime_1 = testData0.Sigmaprime_fun(a_coeffs, 1, q_tilde_1)
%Sigmaprime_func_a_theta_qp1_qt1 = @(a, th, qt1)
delta_numerator_test = sigma_prime_0' * Ifun(q_minus) * sigma_prime_1
delta_denominator_test = sigma_prime_0' * sigma_prime_0
delta_test = delta_numerator_test / delta_denominator_test
%%
% Checking inequalities:
% 12a) S
% existence = d_squared_over_M_minus > 0 && d_squared_over_M_minus < 1;
stability_value = (V_minus_test * delta^2)/(M_minus_test - delta^2)+V_max_test
data.stability_value = stability_value
% stability = (stability_value) < 0;

d_squared_over_M_minus = delta^2/M_minus_test
c1 = -d_squared_over_M_minus % <= 0;
c2 = d_squared_over_M_minus - 1 %<= 0;
c3 = (V_minus_test * delta^2)/(M_minus_test - delta^2)+V_max_test %<= 0

%%
c7 = in_W_inequalities(X0,testData0)
% c = [c1; c2; c3; c4; c5; c6; c7];

%%

c = generate_inequality_constraint(X0,testData0)
ceq = generate_equality_constraint(X0, testData0)
%%
[c, ceq] = generate_nonlinear_constraint(X0, testData0)
%% Run optimzer

problem.x0 = X0;

X_optimized = fmincon(problem);

%%
X_optimized'
X0'
generate_inequality_constraint(X_optimized,testData0)'
generate_equality_constraint(X_optimized, testData0)'
%% Functions


function [a, k] = get_a_k(X)
k=max(size(X)) - 3;
a = X(1:k);
end

function [qt1, qp1] = get_qt1_qp1(X,k)
beta = X(k+1);
q_minus = [0.5*(pi - beta); pi + beta];
q_plus = [0.5*(pi + beta); pi - beta];
q_tilde = q_plus - q_minus;
qt1 = q_tilde(1);
qp1 = q_plus(1);
end

function cost=objective_function(X, data)

[a,k] = get_a_k(X);
beta = X(k+1);
v1 = X(k+2);
v2 = X(k+3);
% Inequalities for (11)
q_minus = [0.5*(pi - beta); pi + beta];
q_plus = [0.5*(pi + beta); pi - beta];
q_tilde = q_plus - q_minus;
q_tilde_1 = q_tilde(1);


a_flipped = flip(a);
phi=@(theta) polyval(a_flipped,theta);
sigma=@(theta) [q_plus(1)- theta*q_tilde_1; phi(theta)];

% sigma_norm= @(theta)norm(sigma(theta));
sigma_norm = @(theta)sqrt((q_plus(1)- theta*q_tilde_1).^2 + phi(theta).^2);


% sum of values not in W1 or W2
theta_vals = linspace(0,1,data.N);
q_vals = data.Sigma_fun(a, theta_vals, q_plus(1), q_tilde_1);
limit_line = -2*q_vals(1,:)+2*pi;
q1_vals = q_vals(1,:);
q2_vals = q_vals(2,:);

% Safe-set Cost from q1-q2 graph
interval1_indices = q1_vals < pi/2;
q2_vals_interval1 = q2_vals(interval1_indices);

a1 = sum(q2_vals_interval1 - limit_line(interval1_indices));

interval2_indices = q1_vals > pi/2;
q2_vals_interval2 = q2_vals(interval2_indices);

a2 = sum(limit_line(interval2_indices) - q2_vals_interval2);

safe_set_cost = abs(a1 - a2);

% safe_set_error = cost1 + cost2 + cost3 + cost4;
qp1 = q_plus(1);
qt1 = q_tilde_1;
B_perp_D_sigma_prime_vals = data.B_perp_D_sigma_prime_fun(a, theta_vals, qp1, qt1);

epsilon_delta = (data.epsilon - min(B_perp_D_sigma_prime_vals.^2));
transversality_cost = (epsilon_delta > 0) * epsilon_delta;

v_cost = v1'*v1 + v2'*v2;

w1 = 1;
w2 = 10;
w3 = 5;
% w3 = 0;
w4 = 0.1;
w5 = 2;
w6 = 50;
beta_error = (data.beta - beta)^2;

[qt1, qp1] = get_qt1_qp1(X,k);
solver_data.Psi1 = @(theta) data.Psi1fun(a, theta, qp1, qt1);
solver_data.Psi2 = @(theta) data.Psi2fun(a, theta, qp1, qt1);

[M,V, MV_vals] = compute_M_V(solver_data);

M_minus = ppval(M,1);
V_minus = ppval(V,1);
V_max = max(MV_vals(:,2));
data.V_max = V_max;

sigma_prime_0 = data.Sigmaprime_fun(a, 0, q_tilde_1);
sigma_prime_1 = data.Sigmaprime_fun(a, 1, q_tilde_1);
%
I = data.Ifun(q_minus);
delta_numerator = sigma_prime_0' * I * sigma_prime_1;
delta_denominator = sigma_prime_0' * sigma_prime_0;
delta = delta_numerator / delta_denominator;

% Checking inequalities:
% 12a) 
% existence = d_squared_over_M_minus > 0 && d_squared_over_M_minus < 1;
stability_value = (V_minus * delta^2)/(M_minus - delta^2)+V_max;
% stability_value
% transversality_cost
% cost = integral(sigma_norm,0,1) + data.stability_value + cost_B_perp_D_sigma_prime;
% cost = w1 * integral(sigma_norm,0,1) + w2 * data.stability_value + w3 * beta_error + w4 * safe_sebetat_error;

% cost = w1 * integral(sigma_norm,0,1) + w2 * data.stability_value + w3 * beta_error + w4 * safe_set_cost + w5 * v_cost + w6 * transversality_cost;

% cost = w1 * integral(sigma_norm,0,1) + 
% data.stability_value
% cost = w2 * data.stability_value + w3 * beta_error +  w5 * v_cost  + w6 * transversality_cost;
cost = w1 * integral(sigma_norm,0,1) + w2 * stability_value + w3 * beta_error +  w6 * transversality_cost; %+ w4 * safe_set_cost;
% 
% stabilitycost = w2 * data.stability_value
% beta_cost = w3 * beta_error
% trans_cost = w6 * transversality_cost
cost = double(cost);

end

function [f_v1] = f_of_v1(v1, q_tilde_1, I)
f_v1 = -(q_tilde_1 * (I(2,2)*v1 - I(2,1)*q_tilde_1))/(I(1,2)*v1 - ...
    I(1,1)*q_tilde_1);
end

%% 

function [ceq] = generate_equality_constraint(X, data)
[a,k] = get_a_k(X);
beta = X(k+1);
v1 = X(k+2);
v2 = X(k+3);

q_minus = [0.5*(pi - beta); pi + beta];
q_plus = [0.5*(pi + beta); pi - beta];
q_tilde = q_plus - q_minus;

I = data.Ifun(q_minus);

f_v1 = f_of_v1(v1, q_tilde(1), I);

% beq = [q_plus(2); q_minus(2); pi; f_v1; v1; v2; data.beta];
beq = [q_plus(2); q_minus(2); pi; f_v1; v1; v2];

a_flipped = flip(a);
phi=@(theta) polyval(a_flipped,theta);
phiprime=@(theta) polyval(polyder(a_flipped),theta);

c1 = phi(0);
c2 = phi(1);
c3 = phi(0.5);
c4 = phiprime(0);
c5 = phiprime(1);
c6 = phiprime(0.5);

% Ax = [c1; c2; c3; c4; c5; c6; beta];
Ax = [c1; c2; c3; c4; c5; c6];

ceq = Ax - beq;
end
%%
function [c] = in_W_inequalities(X, data)

    [a,k] = get_a_k(X);
    [qt1, qp1] = get_qt1_qp1(X,k);
    
    theta_vals = linspace(0,1,data.N);
    q_vals = data.Sigma_fun(a, theta_vals, qp1, qt1);
    limit_line = -2*q_vals(1,:)+2*pi;
        
    q1_vals = q_vals(1,:);
    q2_vals = q_vals(2,:);
       
    % check conditions of first half of curve
    interval1_indices = q1_vals < pi/2;
    q2_vals_interval1 = q2_vals(interval1_indices);
    
    c1 = max(limit_line(interval1_indices)- q2_vals_interval1); %<= 0;
    c2 = max(q2_vals_interval1 - 3*pi); %<=0
    % in_W1 = all(limit_line(interval1_indices) <= q2_vals_interval1) && all(q2_vals_interval1 < 3*pi);
    
    % check conditions of first half of curve
    interval2_indices = q1_vals > pi/2;
    q2_vals_interval2 = q2_vals(interval2_indices);
    
    c3 = max(-pi - q2_vals_interval2); %<= 0
    c4 = max(q2_vals_interval2 - limit_line(interval2_indices));
    % in_W2 = all(-pi < q2_vals_interval2) && all(q2_vals_interval2 <= limit_line(interval2_indices));
    % q2_vals_interval2 < = limit_line(interval2_indices)
    c = [c1;c2;c3;c4];
    c=reshape(c,numel(c),1);

end

function [c] = generate_inequality_constraint(X,data)


    [a,k] = get_a_k(X);
    beta = X(k+1);
    v1 = X(k+2);
    v2 = X(k+3);
    % Inequalities for (11)
    q_minus = [0.5*(pi - beta); pi + beta];
    q_plus = [0.5*(pi + beta); pi - beta];
    q_tilde = q_plus - q_minus;
    q_tilde_1 = q_tilde(1);

    c1 = v1-2*q_tilde_1; % < 0;
    c2 = 2*q_tilde_1 - v2; % < 0;

    % Transversality condition
    c3 = transversality_condition(X, data);

    % Compute virtual Mass and Potential

    % Extract parameters
    [qt1, qp1] = get_qt1_qp1(X,k);
    solver_data.Psi1 = @(theta) data.Psi1fun(a, theta, qp1, qt1);
    solver_data.Psi2 = @(theta) data.Psi2fun(a, theta, qp1, qt1);

    [M,V, MV_vals] = compute_M_V(solver_data);

    % Stability conditions from (12)
    [c4, c5, c6] = existence_and_stability_conditions(X, M, V, MV_vals, data);
    c7 = in_W_inequalities(X,data);
    c = [c1; c2; c3; c4; c5; c6; c7];
 end


 function [c] = transversality_condition(X, data)
    epsilon = data.epsilon;
    N = data.N;
    
    [a,k] = get_a_k(X);
    [qt1, qp1] = get_qt1_qp1(X,k);
      
    theta_vals = linspace(0,1,N+1);
    B_perp_D_sigma_prime_vals = data.B_perp_D_sigma_prime_fun(a, theta_vals, qp1, qt1);

    % % For no zero crossings, 
    % % if max(vals) > 0, min needs to be >0 too
    % % otherwise max <0 => min needs to be < 0 
    % c = 0;
    % max = max(B_perp_D_sigma_prime_vals);
    % min = min(B_perp_D_sigma_prime_vals);
    % if(max > 0)
    %     c = -min; % min > 0 => -min < 0
    % else
    %     % min < 0
    %     c = min;
    % end

    c = epsilon - min(B_perp_D_sigma_prime_vals.^2);
 end

function [M,V,X] = compute_M_V(data)
    % ode_solver_ops = data.ops;
    ode_solver_ops=odeset('RelTol',1e-4,'AbsTol',1e-4);
    % ode_solver_time = data.ode_solver_time; 
    solver_time = linspace(0,1,1000);

    % Note, data parameter needs to contain Psi2 and Psi1, used in mass_potential()
    [Theta,X]=ode45(@mass_potential,solver_time,[1;0],ode_solver_ops,data);
    M=spline(Theta,X(:,1));
    V=spline(Theta,X(:,2));

end

function xdot=mass_potential(theta,x,data2)
    M=x(1);
    V=x(2);

    % Compute Mdot, Vdot
    Mdot = -2*M*data2.Psi2(theta);
    Vdot = -data2.Psi1(theta)*M;
    xdot=[Mdot;Vdot];
end

function [c1, c2, c3] = existence_and_stability_conditions(X, M,V, MV_vals, data)
    Delta = data.Delta;
    [a,k] = get_a_k(X);
    beta = X(k+1);

    q_minus = [0.5*(pi - beta); pi + beta];
    q_plus = [0.5*(pi + beta); pi - beta];
    q_tilde = q_plus - q_minus;
    q_tilde_1 = q_tilde(1);

    I = data.Ifun(q_minus);

    % Computing parameters for inequalities
    M_minus = ppval(M,1);
    V_minus = ppval(V,1);
    V_max = max(MV_vals(:,2));
    data.V_max = V_max;

    sigma_prime_0 = data.Sigmaprime_fun(a, 0, q_tilde_1);
    sigma_prime_1 = data.Sigmaprime_fun(a, 1, q_tilde_1);
    %
    delta_numerator = sigma_prime_0' * I * sigma_prime_1;
    delta_denominator = sigma_prime_0' * sigma_prime_0;
    delta = delta_numerator / delta_denominator;
    
    % Checking inequalities:
    % 12a) 
    % existence = d_squared_over_M_minus > 0 && d_squared_over_M_minus < 1;
    stability_value = (V_minus * delta^2)/(M_minus - delta^2)+V_max;
    data.stability_value = stability_value;
    % stability = (stability_value) < 0;

    d_squared_over_M_minus = delta^2/M_minus;
    c1 = -d_squared_over_M_minus; % <= 0;
    c2 = d_squared_over_M_minus - 1; %<= 0;
    c3 = (V_minus * delta^2)/(M_minus - delta^2)+V_max; %<= 0

end

function [c, ceq] = generate_nonlinear_constraint(X,data)
    c = double(generate_inequality_constraint(X, data));
    ceq = double(generate_equality_constraint(X, data));
    % c = -0.1;
    % ceq = 0;
end

% Test running solver
% TODO - pass in extra function parameters 

function [c, ceq] = parameterfun(X_opt, testData)
    [c, ceq] = generate_nonlinear_constraint(X_opt, testData);
end
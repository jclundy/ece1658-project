%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ECE1658
%% FINAL PROJECT
%% VHC for walking acrobot
%% PARTIAL CODE TO GET STUDENTS STARTED ON THE PROJECT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
clc

number_steps=25;
symbolic_computations=1;
%% Physical parameters
l=1;
lc=0.5;
m=1;
Iz=1/12*m*l^2;
g=9.81;
psi=deg2rad(2);
incline_degrees =0;
%% Control parameters
Kp=20^2;
Kd=20*2;

%% Symbolic computations
if symbolic_computations
    % Define symbolic variables
    fprintf('\n Initializing symbolic variables...\n')
    % syms l lc Iz m g real
    syms t q1 q2 x1 x2 q1dot q2dot x1dot x2dot tau real
    
    q=[q1;q2];x=[x1;x2]; qbar=[q;x];
    qdot=[q1dot;q2dot];  xdot=[x1dot;x2dot]; qbardot=[qdot;xdot];
    
    fprintf('\n Symbolic model computation...\n')
    
    % Define centres of mass of two links
    rc1=x+lc*[cos(q1);sin(q1)];
    rc2=x+l*[cos(q1);sin(q1)]+lc*[cos(q1+q2);sin(q1+q2)];
    
    % Compute time derivatives of centres of mass
    rc1dot=jacobian(rc1,qbar)*qbardot;
    rc2dot=jacobian(rc2,qbar)*qbardot;
    
    % Define the total kinetic energy of the robot
    K=1/2*m*(rc1dot'*rc1dot+rc2dot'*rc2dot)+1/2*Iz*(q1dot^2+(q1dot+q2dot)^2);
    K=simplify(K);
    
    % Extract the square symmetric matrix of the kinetic energy
    Dbar=simplify(hessian(K,qbardot));
    
    % Extract the matrix of the kinetic energy of the pinned robot
    D = Dbar(1:2,1:2);
    
    % Define the potential energy of the pinnedrobot
    P = m*g*(lc*sin(q1)+l*sin(q1)+lc*sin(q1+q2));
    psi=deg2rad(incline_degrees);
    gvector=[cos(psi) -sin(psi);sin(psi) cos(psi)]*[0; -1];
    P=-(m*g*lc*[cos(q1);sin(q1)]+m*g*[l*cos(q1)+lc*cos(q1+q2);l*sin(q1)+lc*sin(q1+q2)])'*gvector;
    P=simplify(P);
    % Input matrix of pinned robot
    B=[0;1];
    
    % Computation of matrix C(q,qdot) of pinned robot
    C = sym(zeros(2,2));
    for i=1:2
        for j=1:2
            for k=1:2
                C(k,j)=C(k,j)+1/2*(diff(D(k,j),q(i))+diff(D(k,i),q(j))-diff(D(i,j),q(k)))*qdot(i);
            end
        end
    end
    
    % Computation of gradient of the potential for the pinned robot
    
    G = jacobian(P,q)';
    
    % Computation of qddot (pinned robot)
    
    %     qddot = simplify(inv(D)*(-C*qdot-G+B*tau));
    Dfun=matlabFunction(D,'Vars',{q});
    Gfun=matlabFunction(G,'Vars',{q});
    Cfun=matlabFunction(C,'Vars',{[q;qdot]});
    fprintf('\n Impact map computation...\n')
    
    %% Impact map
    
    % Relabelling map
    
    T=-eye(2)*q + [pi;2*pi];
    Delta2=[T;-eye(2)*qdot];
    
    % First impact map
    
    px=l*[cos(q1)+cos(q1+q2);sin(q1)+sin(q1+q2)];
    xo=sym(zeros(2,1));
    
    E=[jacobian(px,q), eye(2)];
    Deltaqdot=[eye(2),zeros(2,4)]*inv([Dbar -E';E zeros(2,2)])*...
        [Dbar*[eye(2);jacobian(xo,q)];zeros(2,2)];
    Delta1=eval([q;Deltaqdot*qdot]);
    Delta1=simplify(Delta1);
    
    % Composition of two maps
    Delta=simplify(subs(Delta2,[q;qdot],Delta1));
    Deltafun=matlabFunction(Delta,'Vars',{[q;qdot]});
    save('walking_acrobot_model','D','Dfun','Cfun','Gfun','Deltafun','B');
else
    fprintf('\nLoading robot model...\n')
    load('walking_acrobot_model');
end
%% HERE WRITE YOUR CODE FOR THE VHC DESIGN
% The outcome of this part should be a parameter vector a, whose components
% a_1,ldots,a_k define the polynomial phi_a(theta)=a_1 + a_2 theta + ... +
% a_k theta^(k-1)

fprintf('\nDetermining vhc...\n')
% MY CODE START
% Inputs:
beta = 0.316637;
v1 = -0.894373;
v2 = 1.9;
solve_for_a = true;

% 
% beta = 0.5;
% v1= -1.4888;
% v2 = 1.0;

%
% beta = 0.162296883755300;
% v1 = -0.452800816302641;
% v2 = 1.396485516441244;

% Test 5
% beta = 0.246609826879610;
% v1 = -0.693177151860256;
% v2 = 1.666579596179071;
% a=[2.894982826712655	0.373752829732151	-7.483681788032130	30.062465411634317	-36.009621357182866	13.550304557504635]';
% a_unflipped = a;

%Test 6
% beta = 0.499999999986875;
% v1 = -0.995448687388377;
% v2= 5.813130793156454;

% Test 7
% beta = 0.544733708648113;
% v1 = -5.716594282131112;
% v2 = 4.570889081408835;

% Test 8
% beta=   0.414532310409151;
% v1 =  -9.997856753489156;
% v2 =   4.246310770641121;

% a =[2.845983328271486	0.410715120185887	-7.390164735802243	30.121369810639177	-35.925403324803803	13.374701780400764]';
% a_unflipped = a;
% beta = 0.295609325340617;
% v1 = -0.833609316648596;
% v2 = 1.828470386352731;

% Test 9
% a = [2.705532844725604   3.249590490455112 -29.368744694171887  90.173966152421912 -99.999999999544428  36.857593300299314]';
% a_unflipped = a;

% beta = 0.506466748923118;
% v1 = -0.714215160545731;
% v2 = 3.040744692074665;

% Test 10 (Walks)
% 
% a = [2.747552982544943   0.460240412763022  -7.174554868120053  30.384908371960194 -36.015257100612303  13.132742525684138]';
% a_unflipped = a;
% beta = 0.394039670996995;
% v1 = -1.131459980189369;
% v2 = 2.170720312833271;
% solve_for_a = false;

% 
% Test 11
% a = [2.740209691408301   3.310831198181539 -29.545021722053210  90.440054573900213 -99.999999999996092  36.743711443416068]'
% a_unflipped = a;
% beta = 0.495242040516821;
% v1 = -0.717699683565239;
% v2 = 2.981043469410864;
% solve_for_a = false;


% Test 12 (Walks)

% a = [2.653019389824864   0.481405414272296  -6.916477458928292  30.672729656408208 -36.166730111091027  12.906219015599449]';
% a_unflipped = a;
% beta = 0.488573262421356;
% v1 = -1.469185866864844;
% v2 = 2.519303590793809;
% solve_for_a = false;

% Test 13
% a =   [2.333387716653704   0.462881344738286  -6.085601718089357  31.898187171469729 -35.794497497124738  11.135134254596011]';
% beta = 0.808159502255349;
% v1 = -3.515306672443490;
% v2 = 3.883565855889314;
% a_unflipped = a;
% solve_for_a = false;

% Test 15 (WALKS)
% a = [2.863384231886026   0.395031602410334  -7.420079016807185  30.113209321434038 -35.978296233016877  13.446551169002513 ]';
% beta = 0.278208421662461;
% v1 = -0.785927552601494;
% v2 = 1.772758700720685;
% a_unflipped = a;
% solve_for_a = false;

% Test 16
% a = [2.523197579294350   0.496968421420248  -3.015253857554625  16.912811540665107 -18.534331717251145   5.376595760336740]';
% beta = 0.618395074187240;
% v1 = -2.049452736343263;
% v2 = 2.579343536371052;
% a_unflipped = a;
% solve_for_a = false;

% % Test 17
% a = [2.386763852426788   0.475902935045416  -2.685767702108755  17.572946394641502 -18.413314940874709   4.559726241924818]';
% beta = 0.754803159389710;
% v1 = -3.031022546502935;
% v2 = 3.188188537607566;
% a_unflipped = a;
% solve_for_a = false;

% Test 18
% a = [ 2.097604173255261   0.196932024260636  -4.327382992392389  27.216031626857344 -22.990739887712394   1.993136189669809]';
% beta = 1.043988480329504;
% v1 = -8.807017682450226;
% v2 = 5.409057867426868;
% a_unflipped = a;
% solve_for_a = false;


% Test 19
% a =   [2.974396028703572   0.293445937914470  -7.596127085702170  30.032522282572153 -36.311084631206221  13.915401249644194]';
% beta = 0.167166054446172;
% v1 = -0.467693955986393;
% v2 = 1.414876334896876;
% a_unflipped = a;
% solve_for_a = false;

% Test 20
% a = [1.460164462161743   0.503094067164003  -0.915030931273630   8.183158311432738   0.978593084513387  -4.942968435246349]';
% beta = 0.660406925169085;
% v1 = -8.136294720705347;
% v2 = 3.699650621018530;
% a_unflipped = a;
% solve_for_a = false;

% Test 21
% a =   [2.126418375652961   0.161492962044396  18.827649003614532 -66.783543202174783  96.589901601987734 -46.765153591876597 ]';  
% beta = 1.015173838411510;
% v1 = -9.999999999998076;
% v2 = 2.582325256489500;
% a_unflipped = a;
% solve_for_a = false;


% Test 24
% a =   [2.984704762477705   0.301043926451871  -7.682681849956985  30.124779617042300 -36.404903088549368  13.975302957238812]';
% beta = 0.156820383583068;
% v1 = -0.432117377252332;
% v2 = 1.376936758100289;
% a_unflipped = a;
% solve_for_a = false;

% % Test 25
% a = [2.976594328528107   0.335714390339802  -7.652632710718229  29.809529383833862 -35.879285360909698  13.716579646056296]';
% beta = 0.166356570988570;
% v1 = -0.455197378240314;
% v2 = 1.389073882022778;
% a_unflipped = a;
% solve_for_a = false;

% Test 26 - 'improved gait'

a =   [2.762805677485372   0.450224259171958  -7.160116751212509  30.175234835431542 -35.794643882019443  13.086875017408621]';
beta = 0.378786915911426;
v1 = -1.088503639682148;
v2 = 2.113860411893566;
solve_for_a = false;
a_unflipped = a;

% Test 27
% a =     [2.0945    0.7650    5.0991  -19.3702   25.6086  -11.3095]';
% beta = 0.6320;
% v1 = -1.1722;
% v2 = 1.0073;
% a_unflipped = a;
% solve_for_a = false;

% Test 28

% a = [2.7747    4.9882  -35.3886   95.9009 -100.0000   35.4681]';
% beta = 0.4501;
% v1 = -0.6981;
% v2 = 2.6457;
% a_unflipped = a;
% solve_for_a = false;


% Important values for sigma definition
q_minus = [0.5*(pi - beta); pi + beta];
q_plus = [0.5*(pi + beta); pi - beta];
I = jacobian(Delta(3:4), [q1dot; q2dot]);
I = double(subs(I, [q1, q2], q_minus'));
q_tilde = q_plus - q_minus;
q_tilde_1 = q_tilde(1);

% Check v1, v2 satisfy inequalities from equation 11
condition_11a = le(v1,2*q_tilde_1);

if(not(condition_11a))
    disp("Error: inequality 11a) is not satisfied")
    % return;
end

condition_11b = ge(v2, 2*q_tilde_1);
if(not(condition_11b))
    disp("Error: inequality 11b) is not satisfied")
    % return;
end

val_1 = I(1,1)/I(1,2) * q_tilde_1;
val_2 = (I(2,1) + 2*I(1,1)) / (2*I(1,2)+I(2,2))*q_tilde_1;

condition_11c_1 = v1 > val_1 && v1 > val_2;
condition_11c_2 = v1 < val_1 && v1 < val_2;
condition_11c = condition_11c_1 || condition_11c_2;
if(not(condition_11c))
    disp("Error: inequality 11c) is not satisfied")
    % return;
end

% Check if solution to f(v1) exists
f_v1_exists = not(I(1,:) * [-q_tilde_1; v1] == 0);
if(not(f_v1_exists))
    disp("Error: no solution to f(v1)")
    % return;
end

f_v1 = -(q_tilde_1 * (I(2,2)*v1 - I(2,1)*q_tilde_1))/(I(1,2)*v1 - ...
    I(1,1)*q_tilde_1);

syms th;
phi_vector = [1, th, th^2, th^3, th^4, th^5];
if solve_for_a
    % Solving for vector of polynomial coefficients: a
    k = 6;
   
    phi_vector_prime = diff(phi_vector);
    phi_matrix = zeros(6,6);
    
    phi_matrix(1,:) = subs(phi_vector, th, 0);
    phi_matrix(2,:) = subs(phi_vector, th, 1);
    phi_matrix(3,:) = subs(phi_vector, th, 0.5);
    
    phi_matrix(4,:) = subs(phi_vector_prime, th, 0);
    phi_matrix(5,:) = subs(phi_vector_prime, th, 1);
    phi_matrix(6,:) = subs(phi_vector_prime, th, 0.5);
    
    % My own check - check that Phi matrix is invertible
    invertibility_tolerance = 1e-3;
    Phi_matrix_invertible = gt(abs(det(phi_matrix)),invertibility_tolerance);
    if(not(Phi_matrix_invertible))
        disp("Error: Phi matrix is not invertible")
        return;
    end
    % Solve for a
    b_vector = [q_plus(2); q_minus(2); pi; f_v1; v1; v2];
    a = phi_matrix \ b_vector;
    a_unflipped = a;
end

% MY CODE END

%% HERE WE DEFINE THE FUNCTION phi_a AND ITS DERIVATIVES
a=flip(a);
phi=@(theta) polyval(a,theta);
phiprime=@(theta) polyval(polyder(a),theta);
phipprime=@(theta) polyval(polyder(polyder(a)),theta);

% Using phi and its derivatives, below you should define functions sigma,
% sigmaprime, sigmapprime.

% MY CODE START
% sigma_theta = [q_plus(1)- th.*q_tilde_1; phi];
% sigma = matlabFunction(sigma_theta,'Vars',{th});
sigma=@(theta) [q_plus(1)- theta.*q_tilde_1; phi(theta)];

% sigma_theta_prime = diff(sigma_theta);
% sigmaprime = matlabFunction(sigma_theta_prime,'Vars',{th});
sigmaprime = @(theta) [-q_tilde_1; phiprime(theta)];

% sigma_theta_pprime = diff(sigma_theta);
% sigmapprime = matlabFunction(sigma_theta_pprime,'Vars',{th});
sigmapprime= @(theta) [0; phipprime(theta)];
% MY CODE END

% This is the data structure to be passed to various functions
% You will need to add extra information to it.
data.Kp=Kp;
data.Kd=Kd;
data.D=Dfun;

% MY CODE START
Dinv = inv(D);
Dinvfun = matlabFunction(Dinv,'Vars',{[q1;q2]});
data.Dinv = Dinvfun;
% MY CODE END

data.C=Cfun;
data.G=Gfun;
data.B=B;
data.phi=phi;
data.phiprime=phiprime;
data.phipprime=phipprime;
data.sigma=sigma;
data.sigmaprime=sigmaprime;
data.sigmapprime=sigmapprime;

% MY CODE START
data.q_plus = q_plus;
data.q_minus = q_minus;
data.q_tilde = q_tilde;
% MY CODE END

%% HERE WRITE CODE TO TEST WHETHER YOUR VHC WORKS
% MY CODE START

% Plot resulting VHC
theta_vals = linspace(0,1,101);
q_vals = sigma(theta_vals);
figure();
pt = plot(q_vals(1,:), q_vals(2,:));
xlabel("q1");
ylabel("q2");
title("VHC")
%% Verify it passes through q+, q-, q_bar

% check it passes through q+
if(isapprox(sigma(0),q_plus, 'verytight'))
    disp("Passes check: VHC passes through q+")
else
    disp("Error: VHC does not pass through q+")
    % return;
end

% check it passes through q-
if(isapprox(sigma(1),q_minus, 'verytight'))
    disp("Passes check: VHC passes through q-")
else
    disp("Error: VHC does not pass through q-")
    % return;
end

% check it passes through qbar
qbar_vhc = [pi/2; pi];
if(isapprox(sigma(0.5),qbar_vhc, 'verytight'))
    disp("Passes check: VHC passes through qbar")
else
    disp("Error: VHC does not pass through qbar")
    % return;
end
%% Verify it has correct slopes
slopes = [sigmaprime(0), sigmaprime(0.5), sigmaprime(1);];
slopes = slopes(2,:);

expected_slopes = [f_v1, v2, v1];

if(isapprox(slopes,expected_slopes, 'tight'))
    disp("Passes check: tangent vectors at theta {0, 0.5, 1} have correct slope")
else
    disp("Error: tangent vectors at theta {0, 0.5, 1} do not have correct slope")
    % return;
end


%% Verify curve is contained in W for all theta (0,0.5) U (0.5,1)

limit_line = -2*q_vals(1,:)+2*pi;
figure()
plot(q_vals(1,:), q_vals(2,:));
hold on;
plot(q_vals(1,:), limit_line);
% yline(3*pi);
legend(["VHC", "-2q1+2pi"])
xlabel("q1");
ylabel("q2");

figure()
title("q vs theta")
plot(theta_vals, q_vals);
xlabel("\theta");
ylabel("q");
legend(["q1", "q2"]);

index_of_05=find(theta_vals == 0.5);

q1_vals = q_vals(1,:);
q2_vals = q_vals(2,:);


% check conditions of first half of curve
interval1_indices = q1_vals < pi/2 & q1_vals > 0;
q2_vals_interval1 = q2_vals(interval1_indices);

in_W1 = all(limit_line(interval1_indices) <= q2_vals_interval1) && all(q2_vals_interval1 < 3*pi);

if(in_W1)
    disp("Passes check: VHC is contained in W1")
else
    disp("Error: VHC is not contained in W1")
    % return;
end

% check conditions of first half of curve
interval2_indices = q1_vals > pi/2 & q1_vals < pi;
q2_vals_interval2 = q2_vals(interval2_indices);

q2_vals_2 = q_vals(2,(index_of_05+1):end);
limit2 = limit_line(:,(index_of_05+1):end);

in_W2 = all(-pi < q2_vals_interval2) && all(q2_vals_interval2 <= limit_line(interval2_indices));

if(in_W2)
    disp("Passes check: VHC is contained in W2")
else
    disp("Error: VHC is not contained in W2")
    % return;
end
%% Check whether or not curve is a regular VHC
B_perp = [1, 0];
B_perp_D = B_perp * D;

% a_unflipped = phi_matrix \ b_vector;
phi_a = phi_vector * a_unflipped;
sigma_th = [q_plus(1)- th*q_tilde_1; phi_a];
sigma_th_prime = diff(sigma_th);

B_perp_D_sigma_prime=simplify(subs(B_perp_D,[q1;q2],eval(sigma_th))*eval(sigma_th_prime));
B_perp_D_sigma_prime_fun=matlabFunction(B_perp_D_sigma_prime,'Vars',{th});

B_perp_D_sigma_prime_vals = B_perp_D_sigma_prime_fun(theta_vals);

figure()
plot(theta_vals', B_perp_D_sigma_prime_vals');
xlabel("\theta");
ylabel("B^\perp D \sigma'");

% check for zero crossings:
has_zero_crossings = min(B_perp_D_sigma_prime_vals(:)) < 0 && max(B_perp_D_sigma_prime_vals(:)) > 0;
if(has_zero_crossings)
    disp("Error: Curve is not a regular VHC, B^perp D sigma'(theta) has zero crossings")
    return;
else
    disp("Passes check: Curve is a regular VHC B^perp D sigma'(theta) has no zero crossings")
end
%% Check if constraints induce stable hybrid limit cycle

%% Find Psi1 and Psi2
% denominator=simplify(B_perp*D*[1;phiprime])

% phi=@(theta) polyval(a,theta);
% phiprime=@(theta) polyval(polyder(a),theta);
% phipprime=@(theta) polyval(polyder(polyder(a)),theta);


% denominator = @(theta) B_perp*D*[1;phiprime(theta)];
% G_fun = matlabFunction(subs(G,,'Vars',)
% % B_perp_D_sigma_prime_fun=matlabFunction(B_perp_D_sigma_prime,'Vars',{theta});
% Psi1 = @(theta) B_perp*G/denominator(theta)

denominator=subs(B_perp*D*sigma_th_prime,q,sigma_th);
Psi1=simplify(-subs((B_perp*G)/denominator,q,sigma_th));

sigma_th_pprime = diff(sigma_th_prime);

Psi2=-subs(B_perp*(D*sigma_th_pprime+subs(C,[q;qdot],[sigma_th;sigma_th_prime])*sigma_th_prime)/denominator,q,sigma_th);
Psi2=simplify(Psi2);

Psi1fun = matlabFunction(Psi1 ,'Vars',{th});
Psi2fun = matlabFunction(Psi2 ,'Vars',{th});

data2.Psi1=Psi1fun;
data2.Psi2=Psi2fun;

%% Find virtual mass M(theta) and virtual potential V(theta)
ops2=odeset('RelTol',1e-4,'AbsTol',1e-4);
[Theta,X]=ode45(@mass_potential,linspace(0,1,1000),[1;0],ops2,data2);
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

delta_numerator = sigmaprime(0)' * I * sigmaprime(1);
delta_denominator = sigmaprime(0)' * sigmaprime(0);
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

%% NOW YOU CAN SIMULATE THE ROBOT. PLACE YOUR CONTROLLER INSIDE THE FUNCTION acrobot AT THE END OF THIS SCRIPT
ops= odeset('reltol',1e-7,'abstol',1e-7,'Events',@ground_impact);
dt=1/60; % 60 fps; time increment in simulations and animations

fprintf('\n Simulating...\n')
%% DEFINE THE INITIAL CONDITION [q0;qdot0];

% MY CODE START
% q0 = [pi/3; 5* pi/4];
% q0 = q_minus;
theta_a = 0; theta_b = 1;
theta_0 = 0;
% theta_0 = 0.001;
% theta_0 = 0.7;

theta_dot_a = delta * sqrt(-2*V_minus/(M_minus-delta^2));

max_peturb_angle = 0;
% max_peturb_angle = 10*pi/180; % 1 degree
disturbance = max_peturb_angle*[1; -1];

q0 = sigma(theta_0) + disturbance;

perturb_factor_v = 1;
qdot0 = sigmaprime(theta_0) * theta_dot_a * perturb_factor_v;
% qdot0 = [0;0];
% 
% q0 = sigma(theta_0);
% qdot0 = sigmaprime(theta_0) * theta_dot_a;

% MY CODE END

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
    [t,x,te,xe,ie]=ode45(@(t,x) acrobot(t,x,data),0:dt:10,post_impact_state,ops);    
    % Application of the impact map
    impact_state=xe(end,:)';
    post_impact_state=Deltafun(impact_state);
    T{step}=t;
    X{step}=x;
    Ie{step}=ie;
    Te{step}=te;
    Xe{step}=xe;
end



%% Animation of the simulation results

fprintf('\n Creating video file...\n')
v = VideoWriter("acrobot.avi");
open(v)

fprintf('\n Setting up animation...\n')
figure();

ref=0;time_passed=0;step=1;
Axis=[-1 4 0 2];
Time=text(-1+2,1.8,['time= ','0',' secs,',' step= ',num2str(step)]);
axis(Axis);
stance_leg=line([ref l*cos(q0(1))],[0 l*sin(q0(1))],'color','red','linewidth',2);
swing_leg=line([ref+l*cos(q0(1)) ref+l*cos(q0(1))+l*cos(q0(1)+q0(2))],...
    [l*sin(q0(1)) l*sin(q0(1))+l*sin(q0(1)+q0(2))],'linewidth',2);
fprintf('\n Animation is ready...\n')

animation_slowdown_factor=1; % >1 means slow down
for step=1:length(Ie)
    t=T{step};
    x=X{step};
    xe=Xe{step};
    xe=xe(end,:);
    for k=2:length(t)
        t0=clock;
        drawnow;
        q=x(k,1:2)';
        xdata1=[ref ref+l*cos(q(1))];
        xdata2=[ref+l*cos(q(1)) ref+l*cos(q(1))+l*cos(q(1)+q(2))];
        ydata1=[0 l*sin(q(1))];
        ydata2=[l*sin(q(1)) l*sin(q(1))+l*sin(q(1)+q(2))];
        set(stance_leg,'xdata',xdata1,'ydata',ydata1);
        set(swing_leg,'xdata',xdata2,'ydata',ydata2);
        set(Time,'String',['time= ',num2str(round(time_passed+t(k),1)),' secs,',' step= ',num2str(step)]);
        current_axis=gca;
        if ref>.95*current_axis.XLim(end)
            current_axis.XLim=[.95*ref .95*ref+5];
            Time.Position=[.95*ref+2 1.8 0];
            Axis=axis;
        else
            axis(Axis)
        end
        frame = getframe(gcf);
        writeVideo(v,frame)
        while etime(clock,t0)<animation_slowdown_factor*(t(k)-t(k-1))
        end
    end
    time_passed=time_passed+t(end);
    ref=ref+l*(cos(xe(1))+cos(xe(1)+xe(2)));
end

close(v);


%% FUNCTIONS
function xdot=acrobot(t,x,data)
q1=x(1);
q2=x(2);
q1dot=x(3);
q2dot=x(4);
q=x(1:2);
qdot=x(3:4);
Kp=data.Kp;
Kd=data.Kd;
D=data.D;
Dinv = data.Dinv;
C=data.C;
G=data.G;
B=data.B;
phi=data.phi;
phiprime=data.phiprime;
phipprime=data.phipprime;
% DEFINE YOUR CONTROLLER HERE
% tau=....

% MY CODE START
q_tilde = data.q_tilde;
q_tilde_1 = q_tilde(1);
q_plus = data.q_plus;

% feedback gains
Kp = 400;
Kd = 40;

theta = q_plus(1)/q_tilde_1 - q1/q_tilde_1;
e = q2 - phi(theta);
edot = q2dot + phiprime(theta)*q1dot/q_tilde_1;

K = [Kp, Kd];
H = [sin(e); edot];
v = K*H;
 

LgLf_h = [phiprime(theta)/q_tilde_1, 1]*Dinv(q)*B;
Lf2_h =  -[phiprime(theta)/q_tilde_1, 1]*Dinv(q)*(C(x)*qdot+G(q)) - phipprime(theta)*q1dot^2/q_tilde_1^2;
tau = (-Lf2_h - v)/LgLf_h;

% MY CODE END

qddot=inv(D(q))*(-C(x)*qdot-G(q)+B*tau);
xdot=[qdot;qddot];
end

function [value,isterminal,direction]=ground_impact(t,x)
q1=x(1);
q2=x(2);
% impact occurs when q2 = -2*q1+2*pi
value=q2+2*q1-2*pi;

% We exclude the scuffing point from the impact conditions
if abs(q1-pi/2)<0.01
    isterminal=0;
else
    isterminal=1;
end

% We distinguish between impact on S^+ or S^- by changing the way in which
% ode45 monitors the zero crossing
if abs(q1)<pi/2
    direction=-1;
else
    direction=1;
end
end

function xdot=mass_potential(theta,x,data2)
M=x(1);
V=x(2);

Mdot = -2*M*data2.Psi2(theta);
Vdot = -data2.Psi1(theta)*M;
xdot=[Mdot;Vdot];
end
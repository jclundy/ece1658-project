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
q2ref = pi/6;
q3ref = pi + pi/6;
q4ref = -pi/6; 
q5ref = -pi/6;
qref = [q2ref;q3ref;q4ref;q5ref];

% Define two control gains Kp,Kd for the controller (2)
% placing the roots of the polynomial λ2 + Kdλ + Kp at {−20, −20}
Kd = 40;
Kp = 400;
H = [zeros(4,1) eye(4)];

data.H = H;
data.Kp = Kp;
data.Kd = Kd;
data.qref = qref;

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
a =   [2.762805677485372   0.450224259171958  -7.160116751212509  30.175234835431542 -35.794643882019443  13.086875017408621]';
beta = 0.378786915911426;
v1 = -1.088503639682148;
v2 = 2.113860411893566;


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
alpha0 = alpha_max - alpha_delta;

alpha_2_sym = alpha0 + alpha_delta*sin(theta*pi);
alpha_4_sym = alpha0 - alpha_delta*sin(theta*pi);

gamma_2_sym = (pi-alpha_2_sym )/2;
gamma_4_sym = (pi-alpha_4_sym )/2;

% Save to data structure
data.alpha = alpha0;

% Polynomial
% k = 6;
a_flipped = flip(a);
phi_sym = poly2sym(a_flipped, theta);


%%
% a_flipped=flip(a);
phi=@(theta) polyval(a_flipped,theta);
phiprime=@(theta) polyval(polyder(a_flipped),theta);
phipprime=@(theta) polyval(polyder(polyder(a_flipped)),theta);

data.phi = phi;
data.phiprime= phiprime;
data.phipprime = phipprime;

% virtual q1 and q2 (from acrobot)
qv1_sym = qv_plus(1)- theta*qv_tilde_1; 
qv2_sym = phi_sym;

% target positions for biped
q5_desired = -pi/12;

sigma_q1 = qv1_sym - gamma_2_sym ;
sigma_q2 = pi - alpha_2_sym ;
sigma_q3 = gamma_4_sym  + qv2_sym;
sigma_q4 = pi + alpha_4_sym ;
sigma_q5 = q5_desired;

sigma_sym = [sigma_q1; sigma_q2; sigma_q3; sigma_q4; sigma_q5];

sigma_fun = matlabFunction(sigma_sym,'Vars',theta);
sigma_fun_1 = matlabFunction(sigma_q1,'Vars',theta);
sigma_fun_3 = matlabFunction(sigma_q3,'Vars',theta);
sigma_fun_4 = matlabFunction(sigma_q4,'Vars',theta);

% TODO - re-examine
data.sigma_fun = sigma_fun;
data.q_tilde = qv_tilde;
data.q_plus = qv_plus;

%%
sigma_prime_sym = diff(sigma_sym, theta);
sigma_pprime_sym = diff(sigma_prime_sym, theta);

sigma_prime_fun = matlabFunction(sigma_prime_sym,'Vars',theta);
%% Compute I
q_minus = sigma_fun(1);
q_plus = sigma_fun(0);

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

%% Compute the quantity B⊥ D(σ(θ))σ′(θ) (you’ll find that it is actually a constant, independent of θ),
% and verify that it is not zero.

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
%% Verify curve is contained in W for all theta (0,0.5) U (0.5,1)

q1_vals = sigma_fun_1 (theta_vals);
q3_vals = sigma_fun_3 (theta_vals);

limit_line = -2*q1_vals+2*pi;
figure()
plot(q1_vals, q3_vals);
hold on;
plot(q1_vals, limit_line);
hold off;

% yline(3*pi);
legend(["VHC", "-2q1+2pi"])
xlabel("q1");
ylabel("q4");
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

ops = odeset(RelTol=1e-8,AbsTol=1e-8,Events=@ground_impact);

% initial condition
q0=[pi/3;pi/4;pi-pi/6;-pi/3;-pi/6]; 
qdot0=zeros(5,1);

% theta_0 = 0;
% q0 = sigma_fun(theta_0);
% theta_dot_a = 3.928935686646192;
% qdot0 = sigma_prime_fun(theta_0) * theta_dot_a;

dt = 1/60;
number_steps = 13;

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

animation_slowdown_factor=1; % >1 means slow down
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
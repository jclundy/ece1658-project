% Create an ODE function biped.m accepting the structure array data. The function declaration
% should be this: function xdot=biped(t,x,data). Here, x is the robot state, [q;  ̇q] and xdot is its
% time derivative. The function therefore implements the closed-loop vector field.
% Within this function, extract from x the subvectors q and qdot, and extract from the structure
% data the variables you need for the ODE. You need to compute τ as in (2) and compute  ̈q =
% D−1(q)(−C(q,  ̇q)  ̇q − ∇q P(q) + Bτ). Once you’ve done that, you’ll set xdot = [qdot;qddot].

function xdot=biped(t,x,data)
    % extract variables and functions needed for the ODE
    q = x(1:end/2);
    qdot = x(end/2+1:end);
    Dfun = data.D;
    Cfun = data.C;
    Gfun = data.G;
    B = data.B;
    H = data.H;
    Kp = data.Kp;
    Kd = data.Kd;
    % qref = data.qref;



    %% My changes

    % Compute positions of masses
    l = data.l;
    l1 = l(1); l2 = l(2); l3 = l(3);
    q1 = q(1); q2 = q(2); q3 = q(3); q4 = q(4); q5 = q(5);
    r1 = [0;0];
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

    % Trig to compute qv1, qv2
    v_3_1 = r3 - r1;
    x1 = v_3_1(1); 
    y1 = v_3_1(2); 
    qv1 = atan2(y1, x1);
    % 
    % v_3_5 = r3 - r5;
    % x2 = v_3_5(1);
    % y2 = v_3_5(2);
    % qv2 = atan2(y2, x2) - qv1;

    % Compute theta
    q_tilde = data.q_tilde;
    q_tilde_1 = q_tilde(1);
    q_plus = data.q_plus;
    q_plus_1 = q_plus(1);

    % theta = q_plus_1/q_tilde_1 - qv1 / q_tilde_1;
    % theta = q_plus_1/q_tilde_1 - q(1) / q_tilde_1;
    % theta = double(theta);

    qv_plus_1 = data.qv_plus(1);
    qv_tilde_1 = data.qv_tilde(1);
    theta = qv_plus_1/qv_tilde_1 - qv1 / qv_tilde_1;
    theta = double(theta);
    %% 'Alternate' theta estimation
    % alpha = data.alpha;
    % gamma = (pi-alpha)/2;
    % qv1 = q1 + gamma
    % qv1 = q1 + q2/2;
    % theta = (q_plus_1 - qv1)/q_tilde_1;
    %%  Compute setpoints  
    % Compute setpoints
    sigma_fun = data.sigma_fun;
    qref_bar = sigma_fun(theta);
    qref = qref_bar(2:end);
    % qref = sigma_fun(theta);
    % sigma_prime_fun = data.sigma_prime_fun;
    
    % sigma_prime = sigma_prime_fun(theta);
    % phiprime2 = sigma_prime(2);
    % phiprime3 = sigma_prime(3);
    % phiprime4 = sigma_prime(4);
    % phiprime5 = sigma_prime(5);

    phiprime2 = data.phiprime2;
    phiprime3 = data.phiprime3;
    phiprime4 = data.phiprime4;
    phiprime5 = data.phiprime5;
    % 
    % H = [phiprime2(theta)/q_tilde(1) 1 0 0 0;
    %         phiprime3(theta)/q_tilde(1) 0 1 0 0;
    %         phiprime4(theta)/q_tilde(1) 0 0 1 0;
    %         phiprime5(theta)/q_tilde(1) 0 0 0 1];

    %%
    % matlab warning: INV(A)*b can be slower and less accurate than A\b.
    % Consider using A\b for INV(A)*b

    tau = ( H*( Dfun(q)\B ) )\( H*( Dfun(q) \ ( Cfun(q,qdot)*qdot + Gfun(q) ) ) - Kp*sin(H*q-qref) - Kd*H*qdot  );

    %%
    % % My control
    % q1dot = qdot(1); q2dot = qdot(2); q3dot = qdot(3); q4dot = qdot(4); q5dot = qdot(5);
    % 
    % % Get functions
    % phiprime = data.phiprime;
    % phipprime = data.phipprime;
    % Dinvfun = data.Dinv;
    % 
    % % Run control
    % % K = [Kp, Kd];
    % e = H*q-qref;
    % 
    % edot = [-q2dot; phiprime(theta)*q1dot/q_tilde_1-q3dot; -q4dot; -q5dot];
    % 
    % v = Kp*e + Kd * edot;
    % 
    % Phi_1 = [[0; phiprime(theta)/q_tilde_1; 0; 0], eye(4)];
    % 
    % % LgLf_h = [phiprime(theta)/q_tilde_1, 1]*Dinv(q)*B;
    % Dinv= Dinvfun(q);
    % LgLf_h = Phi_1 * Dinv*B;
    % 
    % % Lf2_h = -[phiprime(theta)/q_tilde_1, 1]*Dinv(q)*(C(x)*qdot+G(q)) - phipprime(theta)*q1dot^2/q_tilde_1^2;
    % Phi_2 = [0; phipprime(theta)*q1dot^2/q_tilde_1^2; 0; 0];
    % Lf2_h = -Phi_1 * Dinv*(Cfun(q,qdot)*qdot + Gfun(q)) + Phi_2;
    % 
    % tau = LgLf_h\(-Lf2_h - v);
%%
    qddot = Dfun(q) \ ( -Cfun(q,qdot)*qdot - Gfun(q) + B*tau );
    xdot= [qdot; qddot];

end

%% 3-DOF arm dynamics (PARA parameters)
syms q1 q2 q3 dq1 dq2 dq3 mp real
q  = [q1; q2; q3];
dq = [dq1; dq2; dq3];

% Parameters
l1 = 0.414; l2 = 0.476; l3 = 0.464;
m2 = 2.77;  m3 = 1.21;
lc2 = l2/2; lc3 = l3/2;
I1 = 0.5;   I2 = 0.0535; I3 = 0.0222;   % I1 is a placeholder
g0 = 9.81;
q23 = q2 + q3;  dq23 = dq2 + dq3;

% Kinetic energy
K1 = 0.5*I1*dq1^2;
K2 = 0.5*(m2*lc2^2 + I2)*(dq2^2 + cos(q2)^2*dq1^2);
K3 = 0.5*m3*( l2^2*dq2^2 + lc3^2*dq23^2 + 2*l2*lc3*cos(q3)*dq2*dq23 ...
            + (l2*cos(q2) + lc3*cos(q23))^2*dq1^2 ) ...
   + 0.5*I3*(dq23^2 + cos(q23)^2*dq1^2);
Kp = 0.5*mp*( l2^2*dq2^2 + l3^2*dq23^2 + 2*l2*l3*cos(q3)*dq2*dq23 ...
            + (l2*cos(q2) + l3*cos(q23))^2*dq1^2 );
K = K1 + K2 + K3 + Kp;

% Potential energy
V = g0*( m2*(l1 + lc2*sin(q2)) ...
       + m3*(l1 + l2*sin(q2) + lc3*sin(q23)) ...
       + mp*(l1 + l2*sin(q2) + l3*sin(q23)) );

% Extract M, C, g
M = simplify(hessian(K, dq));
G = simplify(gradient(V, q));
C = sym(zeros(3));
for i = 1:3
  for j = 1:3
    for k = 1:3
      C(i,j) = C(i,j) + 0.5*( diff(M(i,j),q(k)) + diff(M(i,k),q(j)) ...
                            - diff(M(j,k),q(i)) )*dq(k);
    end
  end
end
C = simplify(C);

%% Checks
disp(simplify(M - M.'))                                  % should be all zeros
disp(subs(G, [q1 q2 q3 mp], [0 0 0 0]))                  % tau2 ~ 14.87 N·m, tau1 = 0
% Passivity check: dM/dt - 2C should be skew-symmetric
dM = sym(zeros(3));
for k = 1:3, dM = dM + diff(M, q(k))*dq(k); end
disp(simplify((dM - 2*C) + (dM - 2*C).'))                % should be all zeros

%% Export numeric function for the optimiser
syms ddq1 ddq2 ddq3 real
ddq = [ddq1; ddq2; ddq3];
tau = M*ddq + C*dq + G;

tau = simplify(tau);

tau1 = tau(1);
tau2 = tau(2);
tau3 = tau(3);

pretty(tau1)      % text version
pretty(tau2)      % text version
pretty(tau3)      % text version


tau_fun = matlabFunction(tau, 'Vars', {q, dq, ddq, mp});
% usage: tau_fun([0.1;0.5;-0.8], [0;0.2;0.1], [0;0;0], 0.5)
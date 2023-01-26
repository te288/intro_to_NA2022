%% Input Parmeters
L = 2*pi; % Length of Reservoir(m)
N = 100; % Number of Control Volume(-)
k = 0.2*ones(1,N); % Permiability(m^2)
phi = ones(1,N); % Porosity(m^2)
c = ones(1,N); % Compressibility(Pa^-1)
mu = 1; % Viscosity of Fluid(Pa^-1)
dx = L/N; % Size of Control Volume(m)
x = zeros(1,N); % x coordinate(m)
x(1) = dx/2;
for i = 2:N
  x(i) = x(i-1) + dx;
end

%% Parameters for Simulation
tmax = 20; % Time to stop simlation (s)
dt   = 0.005;   % dt (s)
nout = 500;  % output result ever nout step(s)
% Variable to decide Boundary condition
% 1 -> Neumann, 0 -> Dirichlet
B_right  = 0;% Boundary at right (x = L)
B_left   = 0;% Boundary at left (x = 0)
Pb_right = 0;% Pressure Value at right (x = L) 
Pb_left  = 0;% Pressure Value at left (x = 0)

%% Initial Conditions
% P_init = ones(N)      % Initial Pressure
P_init = sin(x);

%% Simulation
P_old  = P_init; % Pressure at n-th step
P_new  = P_init; % Pressure at n+1-th step
t = 0;
n = 0;

ax1 = subplot(1,1,1);
plot(ax1, x, P_init, 'DisplayName',[num2str(t,'%05.2f'),'[s]'])
xlim(ax1, [0, L]);
ylim(ax1, [-1,1]);
title(ax1, 'Pressure Diffusion 1D');
hold(ax1, 'on');
while true
    for i = 2:N-1
        alpha = dt / (phi(i)*c(i));
        lam_w = harmmean([k(i-1), k(i)])/mu;
        lam_e = harmmean([k(i+1), k(i)])/mu;
        % A = %-% Write Your Code Here %-%
        % B = %-% Write Your Code Here %-%
        % C = %-% Write Your Code Here %-%
        A = lam_e * alpha / (dx^2);C = lam_w * alpha / (dx^2);
        B = 1 - A - C;
        P_new(i) = A*P_old(i+1) + B*P_old(i) + C*P_old(i-1);
    end
    
    % P(1)
    if B_left == 1 % Neumann Condition
        alpha = dt / (phi(1)*c(1));
        lam_w = harmmean([k(1), k(1)])/mu;
        lam_e = harmmean([k(1+1), k(1)])/mu;
        A = lam_e * alpha / (dx^2);C = lam_w * alpha / (dx^2);
        B = 1 - A - C;
        P_new(1) = A*P_old(2) + B*P_old(1) + C*P_old(1);
    else % Dirichlet Condition
        alpha = dt / (phi(1)*c(1));
        lam_w = harmmean([k(1), k(1)])/mu;
        lam_e = harmmean([k(1+1), k(1)])/mu;
        A = lam_e * alpha / (dx^2);C = lam_w * alpha / (dx^2);
        B = 1 - A - C;
        P_new(1) = A*P_old(2) + B*P_old(1) + C*(2*Pb_right - P_old(1));
    end
    % P(N);
    if B_right == 1 % Neumann Condition;
        alpha = dt / (phi(N)*c(N));
        lam_w = harmmean([k(N-1), k(N)])/mu;
        lam_e = harmmean([k(N), k(N)])/mu;
        A = lam_e * alpha / (dx^2);C = lam_w * alpha / (dx^2);
        B = 1 - A - C;
        P_new(N) = A*P_old(N) + B*P_old(N) + C*P_old(N-1);
    else % Dirichlet Condition ;
        alpha = dt / (phi(N)*c(N));
        lam_w = harmmean([k(N-1), k(N)])/mu;
        lam_e = harmmean([k(N), k(N)])/mu;
        A = lam_e * alpha / (dx^2);C = lam_w * alpha / (dx^2);
        B = 1 - A - C;
        P_new(N) = A*(2*Pb_right- P_old(N)) + B*P_old(N) + C*P_old(N-1);
    end

    % Update Values, time step and Loop counter
    P_old = P_new;
    t = t + dt;
    n = n + 1;
    
    if mod(n,nout) == 0
        plot(ax1, x, P_new, 'DisplayName',[num2str(t,'%05.2f'),'[s]'])
        hold(ax1, 'on');
        disp([num2str(n, '%04d'),'th time step',num2str(t, '%05.2f')])
    end

    if t >= tmax
      break
    end
end

legend(ax1)
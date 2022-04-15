 clear all; close all;
%% Input Parameters
c        = 0.25; % Speed [m/s]
x_left   = -3; % left side of x coordinate[m]
x_right  = 10; % right side of x coordinate[m]
x_num    = 150; % Number of Lattice points[-]
t        = 0; %simulation time[s]
dt       = 0.2; % Simulation time step[s]
t_max    = 40; % time to finish simulation[s]
num_loop = 0; % Loop counter for simulation[-]
num_out  = 25; % plot ever num_out time step[-]


%% Initial Condition
x = linspace(x_left, x_right, x_num); % x coordinate system[m]
phi_init = exp(-x.^2); % initial distribution of phi


%% Simulation
phi_old = phi_init; % phi for nth time step
phi_new = phi_init; % phi for n+1th time step

% plot
ax_exact = subplot(2,1,1);
ax_num   = subplot(2,1,2);
xlim(ax_exact, [x_left, x_right]);
xlim(ax_num, [x_left, x_right]);
ylim(ax_exact, [-0.1, 1.1]);
ylim(ax_num, [-0.1, 1.1]);
% initial plot
plot(ax_exact, x, phi_init, "DisplayName", [num2str(t,'%04.1f'),'s']);
hold(ax_exact,"on");
plot(ax_num, x, phi_init, "DisplayName", [num2str(t,'%04.1f'),'s']);
hold(ax_num, "on");

while true
    % exact solution
    phi_exact = exp(-(x-c*t).^2);

    % numerical solution
    for i = 2:x_num
        phi_new(i) = phi_old(i) - c*dt*(phi_old(i)-phi_old(i-1))/(x(i)-x(i-1));
    end

    % updata
    phi_old = phi_new;
    num_loop = num_loop + 1;
    t = t + dt;


    if mod(num_loop, num_out) == 0
        disp([num2str(num_loop, '%04d'),'th iteration : ', num2str(t,'%04.2f'), 's'])
        plot(ax_exact, x, phi_exact, "DisplayName", [num2str(t,'%04.1f'),'s']);
        hold(ax_exact,"on");
        plot(ax_num, x, phi_new, "DisplayName", [num2str(t,'%04.1f'),'s']);
        hold(ax_exact,"on");
    end

    % judge time step and output
    if (t > t_max) || ( t== t_max)
        break
    end
end

title(ax_exact, 'Exact Solution');
title(ax_num, 'Numerical Solution');
legend(ax_exact, "Location", 'northwest')
legend(ax_num, "Location", 'northwest')
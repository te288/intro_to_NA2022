%% Input Parameters
m    = 1.0; % mass of a small ball [kg]
g    = 9.8; % gravitational acceleration [m/s^2]
k    = 1.0; % coefficient of air resistance [N*s/m]
v0   = 0.0; % Initial velocity [m/s]
dt   = 0.1; % time step for simulation [s]
t_max =  10; % time which simulation is stopped [s]

%% Calculation
t = 0;
i = 1; % counter
v_hist = [v0]; % array to hold velocity value
v_old = v0;
while true
    v = v_old + dt*(g - (k/m)*v_hist(i));
    t = t + dt;
    i = i+1;
    v_old = v;
    
    if t >= t_max
        break
    end    
    v_hist = [v_hist v]; % this code is not actually good.
end

%% Plot result
t_hist = 0:dt:t_max;
plot(t_hist, v_hist);
xlabel('t [s]');
ylabel('v [m/s]');
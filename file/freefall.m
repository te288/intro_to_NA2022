%% Input Parameters
m    = 1.0; % mass of a small ball [kg]
g    = 9.8; % gravitational acceleration [m/s^2]
k    = 1.0; % coefficient of air resistance [N*s/m]
v0   = 0.0; % Initial velocity [m/s]
dt   = 0.1; % time step for simulation [s]
t_max =  10; % time which simulation is stopped [s]

%% simulation
t = 0; % time
i = 1; % counter
v_hist = [v0]; % array to hold velocity value

v_old = v0;
while true
    % Update V, t, i
    v_new = %-% Write Your Code %-%
    t = %-% Write Your Code %-%
    i = %-% Write Your Code %-%
    v_old = %-% Write Your Code %-%
    if t >= t_max
        break
    end
    % Add Data
    v_hist = [v_hist v_new]; 
end

%% Plot result
t_hist = 0:dt:t_max;
scatter(t_hist(1:3:end), v_hist(1:3:end),'DisplayName','Numerical Ans.');
hold on
plot(t_hist, m*g/k*(1-exp(-k*t_hist./m)), 'DisplayName', 'Exact Sol.');
xlabel('t [s]');
ylabel('v [m/s]');
legend()
hold off
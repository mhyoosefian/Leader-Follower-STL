% Mohammad Hussein Yoosefian Nooshabadi
% Formal Methods.
% The code is written and tested on a MacOS machine.

clc
clear
close

%% Problem definition

r = 2; % leader circular path radius
omega_l = pi/18; % angular velocity of the leader
dt = 1; % discretization time
T = 70; % total simulation time
finish = 60; % horizon of the stl: F[start,finish](...)
start = 0;   % horizon of the stl: F[start,finish](...)
x0_l = [0; -2]; % initial location of the leader
x0_f = [-1; -3]; % initial location of the follower
theta0 = -pi/2; % initial angle of the leader
M = 10e2; % for "The Big M" method
eps = 1/M; % for "The Big M" method
uxMax = 0.4; % maximum follower velocity in x direction
uyMax = 0.4; % maximum follower velocity in y direction
uxMin = -0.4; % minimum follower velocity in x direction
uyMin = -0.4; % minimum follower velocity in y direction
a = 0.5; % half-width of unsafe square
b = sqrt(2); % half-width of safe square
robotsWidth = 70; % for visualization purposes
xLeft = -4.5; % coordinates of the vertices of the desired region
xRight = -3; % coordinates of the vertices of the desired region
yBottom = -0.5; % coordinates of the vertices of the desired region
yTop = 0.5; % coordinates of the vertices of the desired region


% Computing discrete locations the leader
N = T/dt; % number of locations excluding the initial location
loc_l = zeros(2, N+1);
loc_l(:, 1) = x0_l;
theta = theta0;
for i=2:N+1
    theta = theta + omega_l*dt;
    loc_l(:, i) = [r*cos(theta); r*sin(theta)];
end


%% Variables of the optimization problem

yalmip('clear');
x = sdpvar(N+1, 1);
y = sdpvar(N+1, 1);
ux = sdpvar(N, 1);
uy = sdpvar(N, 1);
zMu1 = binvar(N+1, 1);
zMu2 = binvar(N+1, 1);
zMu3 = binvar(N+1, 1);
zMu4 = binvar(N+1, 1);
zPsi1 = binvar(N+1, 1);
zPsi2 = binvar(1, 1);
zMu5 = binvar(N+1, 1);
zMu6 = binvar(N+1, 1);
zMu7 = binvar(N+1, 1);
zMu8 = binvar(N+1, 1);
zPsi3 = binvar(N+1, 1);
zPsi4 = binvar(1, 1);

%% Constraints of the optimization problem

mu1 = xRight - x;
mu2 = x - xLeft;
mu3 = yTop - y;
mu4 = y - yBottom;

% The specification related constraints (for visiting the desired reigon) 
Constraints = [];
% For mu1
for i=1:N+1
    Constraints = [Constraints, mu1(i) <= M*zMu1(i)-eps, -mu1(i) <= M*(1-zMu1(i))-eps]; %#ok<*AGROW> 
end
% For mu2
for i=1:N+1
    Constraints = [Constraints, mu2(i) <= M*zMu2(i)-eps, -mu2(i) <= M*(1-zMu2(i))-eps];
end
% For mu3
for i=1:N+1
    Constraints = [Constraints, mu3(i) <= M*zMu3(i)-eps, -mu3(i) <= M*(1-zMu3(i))-eps];
end
% For mu4
for i=1:N+1
    Constraints = [Constraints, mu4(i) <= M*zMu4(i)-eps, -mu4(i) <= M*(1-zMu4(i))-eps];
end

for i=1:N+1 % for conjunction
    Constraints = [Constraints, zPsi1(i) <= zMu1(i), zPsi1(i) <= zMu2(i), zPsi1(i) <= zMu3(i), zPsi1(i) <= zMu4(i), zPsi1(i) >= (1-4+zMu1(i)+zMu2(i)+zMu3(i)+zMu4(i))];
end


if finish-start > T
    error('The simulation time is not enough to check whether the STL is satisfied or not!');
end
for i=start/dt+1:finish/dt+1 % for finally
    Constraints = [Constraints, zPsi2 >= zPsi1(i)];
end
Constraints = [Constraints, zPsi2 <= sum(zPsi1(start/dt+1:finish/dt+1))];

Constraints = [Constraints, zPsi2 == 1]; % for specification satisfaction




% Safety related constraints
% for safety
for i=1:N+1
    Constraints = [Constraints, x(i) <= loc_l(1, i) + b, x(i) >= loc_l(1, i) - b, ... 
        y(i) <= loc_l(2, i) + b, y(i) >= loc_l(2, i) - b];
end
% for collision avoidance
mu5 = x - loc_l(1, :)' - a;
mu6 = -x + loc_l(1, :)' - a;
mu7 = y - loc_l(2, :)' - a;
mu8 = -y + loc_l(2, :)' - a;

% For mu5
for i=1:N+1
    Constraints = [Constraints, mu5(i) <= M*zMu5(i)-eps, -mu5(i) <= M*(1-zMu5(i))-eps]; %#ok<*AGROW> 
end
% For mu6
for i=1:N+1
    Constraints = [Constraints, mu6(i) <= M*zMu6(i)-eps, -mu6(i) <= M*(1-zMu6(i))-eps];
end
% For mu7
for i=1:N+1
    Constraints = [Constraints, mu7(i) <= M*zMu7(i)-eps, -mu7(i) <= M*(1-zMu7(i))-eps];
end
% For mu8
for i=1:N+1
    Constraints = [Constraints, mu8(i) <= M*zMu8(i)-eps, -mu8(i) <= M*(1-zMu8(i))-eps];
end
for i=1:N+1 % for disjunction
    Constraints = [Constraints, zPsi3(i) >= zMu5(i), zPsi3(i) >= zMu6(i), zPsi3(i) >= zMu7(i), zPsi3(i) >= zMu8(i), zPsi3(i) <= zMu5(i)+zMu6(i)+zMu7(i)+zMu8(i)];
end
for i=1:N+1 % for globally
    Constraints = [Constraints, zPsi4 <= zPsi3(i)];
end
Constraints = [Constraints, zPsi4 >= 1 - (N+1) + sum(zPsi3), zPsi4 == 1];




% Dynamics related constraints
Constraints = [Constraints, x(1) == x0_f(1), y(1) == x0_f(2)];
for i=1:N
    Constraints = [Constraints, x(i+1) == x(i) + ux(i)*dt, y(i+1) == y(i) + uy(i)*dt];
end
for i=1:N
    Constraints = [Constraints, ux(i) <= uxMax, uy(i) <= uyMax, ux(i) >= uxMin, uy(i) >= uyMin];
end

%% The objective function

Objective = sum(ux) + sum(uy);

%% Solving the optimization problem

options = sdpsettings('solver', 'mosek', 'verbose', 0);
sol = optimize(Constraints, Objective, options);

%% Plot the result

% Plot position
figure;
desRegion = polyshape([xLeft xLeft xRight xRight], [yBottom yTop yTop yBottom]);
for i=1:N+1
    leader = scatter(loc_l(1, i), loc_l(2, i), robotsWidth, 'k', 'filled', 'square');
    xlim([-5 3]);
    ylim([-3 3]);
    hold on;
    follower = scatter(x(i), y(i), robotsWidth, 'b', 'filled', 'square');
    pgon_b = polyshape([loc_l(1, i)-b loc_l(1, i)+b loc_l(1, i)+b loc_l(1, i)-b], [loc_l(2, i)+b loc_l(2, i)+b loc_l(2, i)-b loc_l(2, i)-b]);
    p_b = plot(pgon_b, 'EdgeColor', 'r', 'FaceAlpha', 0, 'LineWidth', 2);
    pgon_a = polyshape([loc_l(1, i)-a loc_l(1, i)+a loc_l(1, i)+a loc_l(1, i)-a], [loc_l(2, i)+a loc_l(2, i)+a loc_l(2, i)-a loc_l(2, i)-a]);
    p_a = plot(pgon_a, 'EdgeColor', 'r', 'FaceAlpha', 0, 'LineWidth', 2);
    plot(desRegion, "FaceColor", 'green', 'EdgeColor', 'green');
    pause(0.1);
    xlim([-5 3]);
    ylim([-3 3]);
    [in, on] = inpolygon(double(x(i)), double(y(i)), [xLeft xLeft xRight xRight], [yBottom yTop yTop yBottom]);
    if in ~= 0 || on ~= 0
        txt = text(double(x(i)), double(y(i)), 'Visited! $\rightarrow$', 'HorizontalAlignment', 'right', 'FontSize', 16, 'Interpreter','latex');
        pause(1.5);
        delete(txt);
    end
    if i~=N+1
        delete(follower);
        delete(p_b);
        delete(p_a);
    end
end
pause(1.5);

% Plot control actions
figure;
plot(0:dt:T-1, double(ux), 'k-*', 'LineWidth', 2, 'DisplayName', '$u_x$');
hold on
box on
grid on
plot(0:dt:T-1, double(uy), 'k--', 'LineWidth', 2, 'DisplayName', '$u_y$');
plot(0:dt:T-1, ones(N, 1)*uxMax, 'r', 'LineWidth', 2, 'DisplayName', '$u_{max}$');
plot(0:dt:T-1, ones(N, 1)*uxMin, 'b', 'LineWidth', 2, 'DisplayName', '$u_{min}$');
xlabel('Time (s)', 'Interpreter', 'latex');
legend('Interpreter','latex', 'FontSize', 14);

% Plot distance
figure;
hold on;
box on; 
grid on;
distX = abs(double(x) - loc_l(1,:)');
distY = abs(double(y) - loc_l(2,:)');
dist = max(distX, distY);
plot(0:dt:T, dist, 'k', 'LineWidth', 2);
plot(0:dt:T, ones(N+1, 1)*b, 'r', 'LineWidth', 2);
plot(0:dt:T, ones(N+1, 1)*a, 'b', 'LineWidth', 2);
legend('$||X_{leader} - X_{follower}||_{\infty}$', 'Maximum allowable distance', 'Minimum allowable distance', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Time (s)', 'Interpreter', 'latex');
ylim([0.4 1.5]);
ylabel('Distance (m)', 'Interpreter', 'latex');
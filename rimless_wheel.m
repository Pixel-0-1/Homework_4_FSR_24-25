clc
close all
clear all

% Apply default settings for all plots
set(0, 'DefaultTextInterpreter', 'latex')
set(0, 'DefaultLegendInterpreter', 'latex')
set(0, 'DefaultAxesTickLabelInterpreter', 'latex')
set(0, 'DefaultFigureColor', 'w')
lw = 2; % linewidth

% Parameters
g = 9.81;          % gravity (m/s^2)
l = 1.0;           % leg length (m) %% MODIFY HERE, default 1.0 %% 
alpha = pi/8;      % half inter-leg angle (rad) %% MODIFY HERE, default pi/8 %% 
gamma = 0.08;      % slope angle (rad) %%MODIFY HERE, default 0.08 %% 

% Initial conditions
thetadot0 = 0.95; %% MODIFY HERE, default 0.95 %% 
if (thetadot0 >= 0)
    theta0 = gamma-alpha;
else
    theta0 = gamma+alpha;
end

double_support = 0;

y0 = [theta0; thetadot0];

% Simulation settings
t0 = 0; %initial time
tf = 25; %final time
dt = 0.01; %max step time

% Time/state storage
T = [];
Y = [];

while t0 < tf
    options = odeset('Events', @(t, y) impact_event(t, y, alpha,gamma), 'MaxStep', dt);
    [t, y, te, ye, ie] = ode45(@(t, y) dynamics(t, y, g, l, double_support), [t0 tf], y0, options);
    
    T = [T; t];
    Y = [Y; y];

    if ~isempty(te)
        [y0,double_support] = impact_map(ye, alpha,g,l); % apply impact map
        t0 = te;
    else
        break;
    end
end

% Plot results
figure('Renderer', 'painters', 'Position', [10 10 900 700]);
subplot(1,2,1);
plot(T, Y(:,1), 'b', 'DisplayName', '\theta (rad)');
hold on;
plot(T, Y(:,2), 'r', 'DisplayName', '\theta dot (rad/s)');
xlabel('Time (s)');
ylabel('State');
title('Rimless Wheel Dynamics');
legend show;
grid on;

subplot(1,2,2);
plot(Y(:,1), Y(:,2), 'b', 'DisplayName', '\theta (rad)');
hold on
plot(Y(1,1), Y(1,2), 'r', 'Marker','*','LineWidth',5,'DisplayName','Initial point');
xlabel('\theta (rad)');
ylabel('\theta dot (rad/s)');
title('Rimless Wheel Limit Cycle');
legend show;
grid on;

% --- WHEEL ANIMATION ---
% Number of legs (computied form alpha - inter-leg angle)
num_legs = round(pi / alpha);  

figure(2); 
clf; % clean the figure
set(gcf, 'Name', 'Rimless Wheel Animation', 'Position', [100 100 800 600]);

% Visualization range based on leg length
view_range = 1.5 * l;  

% Animation Cycle
hub_x = 0;  % Initial position of the center of the wheel
for k = 1:10:length(T)
    figure(2); 
    cla; 
    
    theta = Y(k, 1); % Central angle of the stance leg
    
    if k > 1
        % Computing the displacement
        delta_theta = Y(k,1) - Y(k-1,1);
        hub_x = hub_x + l * delta_theta * cos(gamma);  % moving along the inclinated plane
    end
    
    % Height of the center of the wheel
    hub_y = hub_x * tan(gamma) + l;
    
    % Draw the inclinated plane
    x_ground = [-view_range:0.1:view_range] + hub_x;
    y_ground = (x_ground - hub_x) * tan(gamma);  % Depends on the slope inclination (gamma)
    plot(x_ground, y_ground, 'k-', 'LineWidth', 3, 'DisplayName', 'Ground'); 
    hold on;
    
    % Draw all the legs
    for i = 0:(num_legs-1)
        % Angle of each leg wrt the local vertical (perpendicular to the plane)
        leg_angle = theta + (2*i - (num_legs - 1)) * alpha;  % Depends on the inter-leg angle (alpha)
        
        % Legs coordinates in the reference frame of the wheel
        x_leg_local = l * sin(leg_angle);
        y_leg_local = -l * cos(leg_angle);
        
        % Transform the coordinates based on gamma
        x_leg_global = hub_x + x_leg_local * cos(gamma) + y_leg_local * sin(gamma);
        y_leg_global = hub_y - x_leg_local * sin(gamma) + y_leg_local * cos(gamma);
        
        % Draw the leg
        plot([hub_x, x_leg_global], [hub_y, y_leg_global], 'b-', 'LineWidth', 2); 
        
        % Highlight the leg that touches the ground
        if abs(leg_angle - gamma) < alpha/2  
            plot([hub_x, x_leg_global], [hub_y, y_leg_global], 'r-', 'LineWidth', 3);
        end
    end
    
    % Central pin of the wheel
    plot(hub_x, hub_y, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
    
    % Infos on parameters and time
    info_text = sprintf(['t = %.2f s\n' ...
                        'Leg length: %.2f m\n' ...
                        'Inter-leg angle: %.1f°\n' ...
                        'Slope angle: %.1f°'], ...
                        T(k), l, alpha*180/pi, gamma*180/pi);
    text(hub_x - view_range*0.9, hub_y + view_range*0.8, info_text, ...
         'FontSize', 10, 'Interpreter', 'latex', 'BackgroundColor', 'white', ...
         'EdgeColor', 'black');
    
    % Axis
    axis equal;
    xlim([hub_x - view_range, hub_x + view_range]);
    ylim([hub_y - view_range*0.7, hub_y + view_range*0.7]);
    xlabel('x (m)', 'Interpreter', 'latex');
    ylabel('y (m)', 'Interpreter', 'latex');
    title('Rimless Wheel Animation on Inclined Plane', 'Interpreter', 'latex');
    grid on;
    
    pause(0.005); 
end

function dydt = dynamics(~, y, g, l, ds)
    theta = y(1);
    thetadot = y(2);
    if (~ds)
        dtheta = thetadot;
        dthetadot = (g/l) * sin(theta);
    else
        dtheta = 0;
        dthetadot = 0;
    end
    dydt = [dtheta; dthetadot];
end

function [value, isterminal, direction] = impact_event(~, y, alpha,gamma)
    
    value = [y(1)-alpha-gamma; y(1)-gamma+alpha];% Trigger when theta = gamma+alpha
                                     %Trigger when theta = gamma-alpha
    isterminal = [1;1];         % Stop the integration
    direction = [1;-1];          % Detect only when increasing
end

function [yplus,ds] = impact_map(y_minus, alpha,g,l)%minus: before impact time; plus: after impact time
    if (y_minus(2)>=0)
        theta_plus = y_minus(1)-2*alpha;
    else
        theta_plus = y_minus(1)+2*alpha;
    end
    thetadot_plus = cos(2*alpha) * y_minus(2);
    if (thetadot_plus < 0.01*sqrt(g/l) && thetadot_plus >-0.01*sqrt(g/l)) 
        thetadot_plus = 0;
        ds = 1;
    else
        ds = 0;
    end
    yplus = [theta_plus; thetadot_plus];
end
function Parameter(Q2)

% Extended Kalman filter for parameter estimation.
% Estimate the natural frequency of a second order system.
% Input Q2 = artificial noise variance for parameter estimation

tf = 100; % simulation length
dt = 0.01; % simulation step size
N = round(tf / dt); % number of simulation time steps
wn = 2; % natural frequency
zeta = 0.1; % damping ratio
b = -2 * zeta * wn;
if ~exist('Q2', 'var')
    Q2 = 0.1; % artificial noise used for parameter estimation
end
Q = [1000 0; 0 Q2]; % covariance of process noise
R = [10 0; 0 10]; % covariance of measurement noise
H = [1 0 0; 0 1 0]; % measurement matrix
P = [0 0 0; 0 0 0; 0 0 20]; % covariance of estimation error
x = [0; 0; -wn*wn]; % initial state
xhat = 2 * x; % initial state estimate
% Initialize arrays for later plotting
xArray = zeros(3, N+1);
xArray(:, 1) = x;
xhatArray = zeros(3, N+1);
xhatArray(:, 1) = xhat;
for index = 1 : N
    % Simulate the system.
    w = sqrt(Q(1,1)) * randn;
    xdot = [x(2) ; x(3)*x(1) + b*x(2) - x(3)*w ; 0];
    x = x + xdot * dt;
    z = H * x + sqrt(R) * [randn ; randn];
    % Simulate the Kalman filter.
    F = [0 1 0 ; xhat(3) b xhat(1) ; 0 0 0];
    L = [0 0 ; -xhat(3) 0 ; 0 1];
    Pdot = F * P + P * F' + L * Q * L' - P * H' / R * H * P;
    P = P + Pdot * dt;
    K = P * H' / R;
    xhatdot = [xhat(2) ; xhat(3)*xhat(1) + b*xhat(2) ; 0];
    xhatdot = xhatdot + K * (z - H * xhat);
    xhat = xhat + xhatdot * dt;
    % Save data for plotting
    xArray(:, index+1) = x;
    xhatArray(:, index+1) = xhat;
end

% Plot results
close all
t = 0 : dt : tf;

figure
plot(t, xArray(3,:) - xhatArray(3,:));
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Seconds'); ylabel('w_n^2 Estimation Error');
title(['Q2 = ', num2str(Q2)]);

disp(['Final estimation error = ', num2str(xArray(3,end)-xhatArray(3,end))]);

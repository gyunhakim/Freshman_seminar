function [AltErr, VelErr, BallErr] = ExtendedBody(dt)

% Continuous time etended Kalman filter example
% Track a body falling through the atmosphere
% Input: 
%   dt = simulation/measurement step size
% Outputs:
%   AltErr = RMS altitude estimation error
%   VelErr = RMS velocity estimation error
%   BallErr = RMS ballistic coefficient estimation error

rho0 = 0.0034; % lb-sec^2/ft^4
g = 32.2; % ft/sec^2
k = 22000; % ft
R = 100; % measurement variance (ft^2)
x = [100000; -6000; 2000]; % initial state
xhat = [100010; -6100; 2500]; % initial state estimate
H = [1 0 0]; % measurement matrix
P = [500 0 0; 0 20000 0; 0 0 250000]; % initial estimation error covariance
tf = 16; % simulation length
if ~exist('dt', 'var')
    dt = 0.0004;
end
N = round(tf / dt);
i = 1;
xArray = zeros(3, N+1);
xArray(:, 1) = x;
xhatArray = zeros(3, N+1);
xhatArray(:, 1) = xhat;
for t = dt : dt : tf
   % Simulate the system (rectangular integration).
   xdot(1,1) = x(2);
   xdot(2,1) = rho0 * exp(-x(1)/k) * x(2)^2 / 2 / x(3) - g;
   xdot(3,1) = 0;
   x = x + xdot * dt;
   % Simulate the measurement.
   z = H * x + sqrt(R) * randn;
   % Simulate the filter.
   temp = rho0 * exp(-xhat(1)/k) * xhat(2) / xhat(3);
   F = [0 1 0; -temp * xhat(2) / 2 / k temp ...
      -temp * xhat(2) / 2 / xhat(3); ...
         0 0 0];
   Pdot = F * P + P * F' - P * H' / R * H * P;
   P = P + Pdot * dt;
   K = P * H' / R;
   xhatdot(1,1) = xhat(2);
   xhatdot(2,1) = temp * xhat(2) / 2 - g;
   xhatdot(3,1) = 0;
   xhatdot = xhatdot + K * (z - H * xhat);
   xhat = xhat + xhatdot * dt;
   % Save data for plotting
   i = i + 1;
   xArray(:, i) = x;
   xhatArray(:, i) = xhat;

end

% Plot data.
close all;
t = 0 : dt : tf;

figure;
subplot(3,1,1);
plot(t, xArray(1,:) - xhatArray(1,:));
title('Estimation Errors');
set(gca,'FontSize',12); set(gcf,'Color','White');
ylabel('Pos. (feet)');

AltErr = std(xArray(1,:) - xhatArray(1,:));
disp(['Continuous EKF RMS altitude estimation error = ', num2str(AltErr)]);

subplot(3,1,2);
plot(t, xArray(2,:) - xhatArray(2,:));
set(gca,'FontSize',12); set(gcf,'Color','White');
ylabel('Vel. (feet/sec)');

VelErr = std(xArray(2,:) - xhatArray(2,:));
disp(['Continuous EKF RMS velocity estimation error = ', num2str(VelErr)]);

subplot(3,1,3);
plot(t, xArray(3,:) - xhatArray(3,:));
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Seconds'); ylabel('Ball. Coeff.');

BallErr = std(xArray(3,:) - xhatArray(3,:));
disp(['Continuous EKF RMS ballistic coefficient estimation error = ', num2str(BallErr)]);

figure;
subplot(2,1,1);
plot(t, xArray(1,:));
title('Falling Body Simulation', 'FontSize', 12);
set(gca,'FontSize',12); set(gcf,'Color','White');
ylabel('True Position (feet)');

subplot(2,1,2);
plot(t, xArray(2,:));
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Seconds'); ylabel('True Velocity (feet/sec)');

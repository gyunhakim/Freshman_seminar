function [x1Array, xhatArray, KArray] = DiscreteKFEx2(Q)

% Simulate a discrete-time Kalman filter with a mismodel condition.
% We design the filter under the assumption that we are trying to estimate
% a constant bias state. In reality the state is a ramp.
% Performance can be improved by progressively increasing Q from 0 to 10.
% This illustrates the effectiveness of adding fictitious process noise 
% to the assumed model. Fictitious process noise can compensate for modeling errors.
% INPUTS:
%   Q = fictitious process noise

if ~exist('Q', 'var')
    Q = 0;
end
tf = 50; % final time
F = 1; % Assumed state matrix
H = 1; % Assumed measurement matrix
x = [0; 10]; % initial true state. x(2) is a constant.
xhat = 0; % initial estimate of x(1)
Pplus = 1; % initial estimation error covariance
R = 1; % covariance of measurement noise

x1Array = zeros(1, tf);
xhatArray =zeros(1, tf);
KArray = zeros(1, tf);
PArray = zeros(1, tf);
ResArray = zeros(1, tf);
for t = 1 : tf
   % Simulate the system
   x(1) = x(1) + x(2);
   y = x(1) + sqrt(R) * randn;
   % Kalman filter
   Pminus = F * Pplus * F' + Q;
   K = Pminus * H' / (H * Pminus * H' + R);
   xhat = F * xhat;
   ResArray(t) = y - xhat;
   xhat = xhat + K * (y - H * xhat);
   Pplus = (1 - K * H) * Pminus * (1 - K * H)' + K * R * K';
   % Save data for later
   x1Array(t) = x(1);
   xhatArray(t) = xhat(1);
   KArray(t) = K;
   PArray(t) = Pplus;
end

% Plot results
close all; t = 1 : tf;

figure
plot(t, x1Array, 'r-', t, xhatArray, 'b--')
title(['Q = ',num2str(Q)], 'FontSize', 14)
set(gca,'FontSize',12); set(gcf,'Color','White')
xlabel('time')
legend('true state', 'estimated state')

figure
plot(t, KArray)
title(['Kalman Gain; Q = ', num2str(Q)], 'FontSize', 14)
set(gca,'FontSize',12); set(gcf,'Color','White')
xlabel('time')
axis([0 tf 0 1])

figure
plot(t, PArray)
title(['Variance; Q = ', num2str(Q)], 'FontSize', 14)
set(gca,'FontSize',12); set(gcf,'Color','White')
xlabel('time')
axis([0 tf 0 1])

figure
plot(t, ResArray)
title(['Residual; Q = ', num2str(Q)], 'FontSize', 14)
set(gca,'FontSize',12); set(gcf,'Color','White')
xlabel('time')

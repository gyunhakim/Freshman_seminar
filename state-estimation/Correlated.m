function Correlated(M, MFilter)

% Kalman filter simulation using correlated process and measurement noise.
% This illustrates the improvement in filter results that can be attained
% when the correlation is taken into account.
% M = correlation between process noise and measurement noise.
% MFilter = value of M used in Kalman filter.

if ~exist('M', 'var')
    M = 0;
end
if ~exist('MFilter', 'var')
    MFilter = 0;
end

kf = 500; % number of time steps in simulation
phi = 0.8; % system matrix
H = 1; % measurement matrix
Q = 1; % process noise covariance
R = 0.1; % measurement noise covariance

% Compute the eigendata of the covariance matrix so that we can simulate correlated noise.
Q1 = [Q M; M' R];
[d, lambda] = eig(Q1);
if (lambda(1) < 0) || (lambda(2) < 0)
    error('Q1 is not positive semidefinite');
end
ddT = d * d';
if (norm(eye(size(ddT)) - ddT) > eps)
   error(['d is not orthogonal. d = ',d]);
end
   
x = 0; % initial state
xhatplus = x; % initial state estimate
Pplus = 0; % initial uncertainty in state estimate
xArray = zeros(1, kf);
xhatArray = zeros(1, kf);
KArray = zeros(1, kf);
PArray = zeros(1, kf);
zArray = zeros(1, kf);

rng(0); % initialize random number generator

for k = 1 : kf
   % Generate correlated process noise (w) and measurement noise (n)
   v = [sqrt(lambda(1,1))*randn; sqrt(lambda(2,2))*randn];
   Dv = d * v;
   w = Dv(1);
   n = Dv(2);
   % Simulate the system dynamics and the measurement
   x = phi * x + w;
   z = H * x + n;
   % Simulate the Kalman filter
   Pminus = phi * Pplus * phi' + Q;
   K = (Pminus * H' + MFilter) / (H * Pminus * H' + H * MFilter + MFilter' * H' + R);
   xhatminus = phi * xhatplus;
   xhatplus = xhatminus + K * (z - H * xhatminus);
   Pplus = Pminus - K * (H * Pminus + MFilter');
   % Save data for plotting
   xArray(k) = x;
   xhatArray(k) = xhatplus;
   KArray(k) = K;
   PArray(k) = Pplus;
   zArray(k) = z;
end

% Plot the data
k = 1 : kf;
close all;

figure;
plot(k, zArray - xArray, 'r-', k, xhatArray - xArray, 'b:');
set(gca,'FontSize',12); set(gcf,'Color','White');
title(['M = ',num2str(M),', MFilter = ',num2str(MFilter)]);
xlabel('Time');
legend('Measurement Error', 'Estimation Error');

figure;
plot(k, KArray);
set(gca,'FontSize',12); set(gcf,'Color','White');
title(['Kalman Gain for M = ',num2str(M),', MFilter = ',num2str(MFilter)]);
xlabel('Time');

figure;
plot(k, PArray);
set(gca,'FontSize',12); set(gcf,'Color','White');
title(['Estimation Error Covariance for M = ',num2str(M),', MFilter = ',num2str(MFilter)]);
xlabel('Time');

% Compute statistics
err = zArray - xArray;
err = norm(err)^2 / kf;
disp(['Measurement Error Variance = ',num2str(err)]);
err = xhatArray - xArray;
err = norm(err)^2 / kf;
disp(['Estimation Error Variance = ',num2str(err)]);
disp(['Analytical Variance = ',num2str(Pplus)]);

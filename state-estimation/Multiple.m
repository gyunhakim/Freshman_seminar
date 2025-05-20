function Multiple

% Multiple-model Kalman filtering
% Second-order system

zeta = 0.1; % damping ratio
wn = [sqrt(4); sqrt(4.4); sqrt(4.8)]; % possible wn values
T = 0.1; % sample period
tf = 60; % Length of simulation
pr = [0.1; 0.6; 0.3]; % a priori probabilities
Q = diag([0, 1600]); % discrete time process noise variance
R = diag([10 10]); % discrete time measurement noise covariance
H = eye(2); % measurement matrix
x = [0; 0]; % initial state

N = size(wn, 1); % number of parameter sets
wnhat = sum(wn .* pr); % Compute the initial estimate of wn
Nsteps = round(tf / T);
q = size(H, 1); % number of measurements
phi = zeros(2, 2, N);
Pplus = zeros(2, 2, N);
Pminus = zeros(2, 2, N);
xhat = zeros(2, N);
pdf = zeros(N, 1);
% Compute the alternative phi matrices
for i = 1 : N
   Fi = [0 1; -wn(i)^2 -2*zeta*wn(i)];
   phii = expm(Fi*T);
   phi(:,:,i) = phii;
   xhat(:,i) = x;
end

% Create arrays for later plotting
wnhatArray = zeros(1, Nsteps+1);
wnhatArray(1) = wnhat;
prArray = zeros(3, Nsteps+1);
prArray(:, 1) = pr;
index = 1;
for t = T : T : tf
   % Simulate the system
   % The first parameter set is the true parameter set
   w = sqrt(Q) * randn(2, 1);
   x = phi(:,:,1) * x + w;
   z = H * x + sqrt(R) * [randn; randn];
   % Run a separate Kalman filter for each parameter set
   for i = 1 : N
      Pminus(:,:,i) = phi(:,:,i) * Pplus(:,:,i) * phi(:,:,i)' + Q;
      K = Pminus(:,:,i) * H' / (H * Pminus(:,:,i) * H' + R);
      xhat(:,i) = phi(:,:,i) * xhat(:,i);
      r = z - H * xhat(:,i); % measurement residual
      S = H * Pminus(:,:,i) * H' + R; % covariance of measurement residual
      pdf(i) = exp(-r' / S * r/2) / (2*pi)^(q/2) / sqrt(det(S));
      xhat(:,i) = xhat(:,i) + K * (z - H * xhat(:,i));
      Pplus(:,:,i) = (eye(2) - K * H) * Pminus(:,:,i) * (eye(2) - K * H)' + K * R * K';
   end
   % Update the probability of each parameter set
   pr = pdf .* pr / sum(pdf .* pr);
   % Compute the best state estimate and the best parameter estimate
   xhatbest = zeros(2, 1);
   wnhat = 0;
   for i = 1 : N
      xhatbest = xhatbest + pr(i) * xhat(:,i);
      wnhat = wnhat + pr(i) * wn(i);
   end
   % Save data for plotting
   index = index + 1;
   wnhatArray(index) = wnhat;
   prArray(:, index) = pr;
end

close all
t = 0 : T : tf;
figure
plot(t, wnhatArray.^2);
title('Estimate of square of natural frequency', 'FontSize', 12);
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Seconds');

figure
plot(t, prArray(1,:), 'b-', t, prArray(2,:), 'k--', t, prArray(3,:), 'r:');
title('Probabilities of square of natural frequency', 'FontSize', 12);
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Seconds');
legend('Probability that \omega_n^2 = 4', 'Probability that \omega_n^2 = 4.4', 'Probability that \omega_n^2 = 4.8');

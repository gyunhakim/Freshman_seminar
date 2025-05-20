function HybridUKF(T)

% Hybrid unscented Kalman filter example
% Track a body falling through the atmosphere
% This example is taken from [Jul00], which was based on [Ath68].
% INPUT: T = sample period

rho0 = 2; % lb-sec^2/ft^4
g = 32.2; % ft/sec^2
k = 2e4; % ft
R = 10^4; % measurement noise variance (ft^2)
Q = diag([0 0 0]); % process noise covariance
M = 10^5; % horizontal range of position sensor
a = 10^5; % altitude of position sensor
x = [3e5; -2e4; 1e-3]; % initial state
xhat = [3e5; -2e4; 1e-3]; % initial EKF state estimate
xhatukf = xhat; % initial UKF state estimation
P0 = diag([1e6 4e6 10]); % initial EKF covariance
P = P0;
Pukf = P; % initial UKF covariance
if ~exist('T', 'var')
    T = 0.5; % measurement time step
end
tf = 30; % simulation length (seconds)
dt = 0.001; % time step for system dynamics integration (seconds)
N = round(tf / T);
xArray = zeros(3, N+1);
xArray(:, 1) = x;
xhatArray = zeros(3, N+1);
xhatArray(:, 1) = xhat;
xhatukfArray = zeros(3, N+1);
xhatukfArray(:, 1) = xhatukf;
Parray = zeros(3, N+1);
Parray(:, 1) = diag(P);
Pukfarray = zeros(3, N+1);
Pukfarray(:, 1) = diag(Pukf);
W = ones(6,1) / 6; % UKF weights
for index = 1 : N
    % Simulate the system
    for tau = dt : dt : T
        xdot(1,1) = x(2);
        xdot(2,1) = rho0 * exp(-x(1)/k) * x(2)^2 / 2 * x(3) - g;
        xdot(3,1) = 0;
        xdot = xdot + sqrt(dt * Q) * [randn; randn; randn];
        x = x + xdot * dt;
    end
    % Simulate the noisy measurement
    z = sqrt(M^2 + (x(1)-a)^2) + sqrt(R) * randn;
    % Simulate the continuous-time part of the Kalman filter (time update)
    for tau = dt : dt : T
        xhatdot(1,1) = xhat(2);
        xhatdot(2,1) = rho0 * exp(-xhat(1)/k) * xhat(2)^2 / 2 * xhat(3) - g;
        xhatdot(3,1) = 0;
        xhat = xhat + xhatdot * dt;
        F = [0 1 0; -rho0 * exp(-xhat(1)/k) * xhat(2)^2 / 2 / k * xhat(3) ...
            rho0 * exp(-xhat(1)/k) * xhat(2) * xhat(3) ...
            rho0 * exp(-xhat(1)/k) * xhat(2)^2 / 2; ...
            0 0 0];
        %H = [(xhat(1) - a) / sqrt(M^2 + (xhat(1)-a)^2) 0 0];
        Pdot = F * P + P * F' + Q * dt;% - P * H' / R * H * P;
        P = P + Pdot * dt;
        if sum(sum(isinf(P))) + sum(sum(isnan(P))) > 0
            P = P0; % reinitialize the EKF covariance
        end
    end
    % Simulate the discrete-time part of the Kalman filter (measurement update)
    H = [(xhat(1) - a) / sqrt(M^2 + (xhat(1)-a)^2) 0 0];
    K = P * H' / (H * P * H' + R);
    zhat = sqrt(M^2 + (xhat(1)-a)^2);
    xhat = xhat + K * (z - zhat);
    P = (eye(3) - K * H) * P * (eye(3) - K * H)' + K * R * K';
    % Generate the UKF sigma points
    [root, ~] = chol(3*Pukf);
    xbreve(:, 1:3) = xhatukf + root(1:3,:)'; % Eq. 14.56
    xbreve(:, 4:6) = xhatukf - root(1:3,:)';
    % UKF time update
    for i = 1 : 6
        for tau = dt : dt : T
            xbrevedot(1,1) = xbreve(2,i);
            xbrevedot(2,1) = rho0 * exp(-xbreve(1,i)/k) * xbreve(2,i)^2 / 2 * xbreve(3,i) - g;
            xbrevedot(3,1) = 0;
            xbreve(:, i) = xbreve(:, i) + xbrevedot * dt; % Eq. 14.57
        end
    end
    xhatukf = zeros(3,1);
    for i = 1 : 6
        xhatukf = xhatukf + W(i) * xbreve(:,i); % Eq. 14.58
    end
    Pukf = zeros(3,3);
    for i = 1 : 6
        Pukf = Pukf + W(i) * (xbreve(:,i) - xhatukf) * (xbreve(:,i) - xhatukf)';
    end
    Pukf = Pukf + Q; % Eq. 14.58
    if sum(sum(isinf(Pukf))) + sum(sum(isnan(Pukf))) > 0
        Pukf = P0; % reinitialize the UKF covariance
    end
    % UKF measurement update
    zukf = sqrt(M^2 + (xbreve(1,:)-a).^2); % Eq. 14.61
    zhat = 0;
    for i = 1 : 6
        zhat = zhat + W(i) * zukf(:,i); % Eq. 14.62
    end
    Py = 0;
    Pxy = zeros(3,1);
    for i = 1 : 6
        Py = Py + W(i) * (zukf(:,i) - zhat) * (zukf(:,i) - zhat)'; % Eq. 14.63
        Pxy = Pxy + W(i) * (xbreve(:,i) - xhatukf) * (zukf(:,i) - zhat)'; % Eq. 14.64
    end
    Py = Py + R;
    Kukf = Pxy / Py;
    xhatukf = xhatukf + Kukf * (z - zhat);
    Pukf = Pukf - Kukf * Py * Kukf';
    if sum(sum(isinf(Pukf))) + sum(sum(isnan(Pukf))) > 0
        Pukf = P0; % reinitialize the UKF covariance
    end
    % Save data for plotting
    xArray(:, index+1) = x;
    xhatArray(:, index+1) = xhat;
    xhatukfArray(:, index+1) = xhatukf;
    Parray(:, index+1) = diag(P);
    Pukfarray(:, index+1) = diag(Pukf);
end

close all
t = 0 : T : tf;
figure; subplot(3,1,1)
semilogy(t, abs(xArray(1,:) - xhatArray(1,:)), 'b'); hold;
semilogy(t, abs(xArray(1,:) - xhatukfArray(1,:)), 'r:');
set(gca,'FontSize',12); set(gcf,'Color','White');
ylabel('Pos Est Err');
legend('Kalman filter', 'Unscented filter');

subplot(3,1,2)
semilogy(t, abs(xArray(2,:) - xhatArray(2,:)), 'b'); hold;
semilogy(t, abs(xArray(2,:) - xhatukfArray(2,:)), 'r:');
set(gca,'FontSize',12); set(gcf,'Color','White');
ylabel('Vel Est Err');

subplot(3,1,3)
semilogy(t, abs(xArray(3,:) - xhatArray(3,:)), 'b'); hold;
semilogy(t, abs(xArray(3,:) - xhatukfArray(3,:)), 'r:');
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Seconds');
ylabel('Ball Coeff Est Err');

figure; subplot(2,1,1)
plot(t, xArray(1,:));
set(gca,'FontSize',12); set(gcf,'Color','White');
ylabel('True Position');

subplot(2,1,2)
plot(t, xArray(2,:));
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Seconds');
ylabel('True Velocity');

% Calculate and output RMS estimation errors
KalPos = std(xArray(1,:) - xhatArray(1,:));
UKFPos = std(xArray(1,:) - xhatukfArray(1,:));
KalVel = std(xArray(2,:) - xhatArray(2,:));
UKFVel = std(xArray(2,:) - xhatukfArray(2,:));
KalBal = std(xArray(3,:) - xhatArray(3,:));
UKFBal = std(xArray(3,:) - xhatukfArray(3,:));
disp('RMS Estimation Errors:');
disp(['  Position = ', num2str(KalPos), ' (Kalman), ', num2str(UKFPos), ' (UKF)']);
disp(['  Velocity = ', num2str(KalVel), ' (Kalman), ', num2str(UKFVel), ' (UKF)']);
disp(['  Ball. Coeff. = ', num2str(KalBal), ' (Kalman), ', num2str(UKFBal), ' (UKF)']);

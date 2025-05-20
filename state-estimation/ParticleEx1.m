function ParticleEx1(N)

% Particle filter example, adapted from Gordon, Salmond, and Smith paper.
% INPUT: N = number of particles

x = 0.1; % initial state
Q = 1; % process noise covariance
R = 1; % measurement noise covariance
tf = 50; % simulation length (time steps)
if ~exist('N', 'var')
    N = 100; % number of particles in the particle filter
end
xhat = x;
P = 2;
xhatPart = x;
% Initialize the particle filter
xpart = x + sqrt(P) * randn(N, 1);
xpartminus = zeros(N, 1);
q = zeros(N, 1);
% Initialize arrays
xArr = zeros(1, tf);
xArr(1) = x;
yArr = zeros(1, tf);
yArr(1) = x^2 / 20 + sqrt(R) * randn;
xhatArr = zeros(1, tf);
xhatArr(1) = x;
PArr = zeros(1, tf);
PArr(1) = P;
xhatPartArr = zeros(1, tf);
xhatPartArr(1) = xhatPart;
close all
for k = 1 : tf
    % System simulation
    x = 0.5 * x + 25 * x / (1 + x^2) + 8 * cos(1.2*(k-1)) + sqrt(Q) * randn;
    y = x^2 / 20 + sqrt(R) * randn;
    % Extended Kalman filter
    F = 0.5 + 25 * (1 - xhat^2) / (1 + xhat^2)^2;
    P = F * P * F' + Q;
    H = xhat / 10;
    K = P * H' * (H * P * H' + R)^(-1);
    xhat = 0.5 * xhat + 25 * xhat / (1 + xhat^2) + 8 * cos(1.2*(k-1));
    xhat = xhat + K * (y - xhat^2 / 20);
    P = (1 - K * H) * P;
    % Particle filter
    for i = 1 : N
        xpartminus(i) = 0.5 * xpart(i) + 25 * xpart(i) / (1 + xpart(i)^2) + 8 * cos(1.2*(k-1)) + sqrt(Q) * randn;
        ypart = xpartminus(i)^2 / 20;
        vhat = y - ypart;
        q(i) = (1 / sqrt(R) / sqrt(2*pi)) * exp(-vhat^2 / 2 / R);
    end
    % Normalize the likelihood of each a priori estimate
    q = q / sum(q);
    % Resample
    for i = 1 : N
        u = rand; % uniform random number between 0 and 1
        qtempsum = 0;
        for j = 1 : N
            qtempsum = qtempsum + q(j);
            if qtempsum >= u
                xpart(i) = xpartminus(j);
                break
            end
        end
    end
    % The particle filter estimate is the mean of the particles
    xhatPart = mean(xpart);
    % Plot the estimated pdf's at a specific time
    if k == 20
        % Particle filter pdf
        pdf = zeros(81, 1);
        for m = -40 : 40
            for i = 1 : N
                if (m <= xpart(i)) && (xpart(i) < m+1)
                    pdf(m+41) = pdf(m+41) + 1;
                end
            end
        end
        figure
        m = -40 : 40;
        plot(m, pdf / N, 'r')
        hold
        title('Estimated pdf at k=20')
        % Kalman filter pdf
        pdf = (1 / sqrt(P) / sqrt(2*pi)) .* exp(-(m - xhat).^2 / 2 / P);
        plot(m, pdf, 'b');
        plot(x, 0, 'k*');
        legend('Particle filter', 'Kalman filter', 'True state');
    end
    % Save data in arrays for later plotting
    xArr(k+1) = x;
    yArr(k+1) = y;
    xhatArr(k+1) = xhat;
    PArr(k+1) = P;
    xhatPartArr(k+1) = xhatPart;
end
t = 0 : tf;

figure
plot(t, xArr, 'b.', t, xhatArr, 'k-', t, xhatArr-2*sqrt(PArr), 'r:', t, xhatArr+2*sqrt(PArr), 'r:');
axis([0 tf -40 40]);
set(gca,'FontSize',12); set(gcf,'Color','White'); 
xlabel('time step'); ylabel('state');
legend('True state', 'EKF estimate', '95% confidence region'); 

figure
plot(t, xArr, 'b.', t, xhatPartArr, 'k-');
set(gca,'FontSize',12); set(gcf,'Color','White'); 
xlabel('time step'); ylabel('state');
legend('True state', 'Particle filter estimate'); 

xhatRMS = sqrt((norm(xArr - xhatArr))^2 / tf);
xhatPartRMS = sqrt((norm(xArr - xhatPartArr))^2 / tf);
disp(['Kalman filter RMS error = ', num2str(xhatRMS)]);
disp(['Particle filter RMS error = ', num2str(xhatPartRMS)]);

function HinfEx1b(SSFlag, BiasFlag, theta)

% Compare the Kalman and H-infinity estimators for a scalar system.
% SSFlag = true if the steady state estimator is used.
% BiasFlag = true if the process noise is biased
% theta = H-infinity performance bound

if ~exist('SSFlag', 'var')
    SSFlag = true;
end
if ~exist('BiasFlag', 'var')
    BiasFlag = true;
end
if ~exist('theta', 'var')
    theta = 1/10;
end
    
if SSFlag
    KK = (1+sqrt(5))/(3+sqrt(5)); % steady-state Kalman gain
    PH = (1-theta+sqrt((theta-1)*(theta-5)))/2/(1-theta);
    KH = PH / (1 - theta * PH + PH); % steady-state H-infinity gain
end
    
kf = 20;
PK = 1;
PH = PK;
Q = 1;
R = 1;
x = 0;
xhatK = x + sqrt(PK); % initial Kalman filter estimate (in error by one sigma)
xhatH = xhatK; % initial H-infinity estimate
xArr = zeros(1, kf+1);
xArr(1) = x;
xhatKArr = zeros(1, kf+1);
xhatKArr(1) = xhatK;
xhatHArr = zeros(1, kf+1);
xhatHArr(1) = xhatH;
for k = 1 : kf
    y = x + sqrt(R) * randn;
    if ~SSFlag
        KK = PK / (1 + PK/R) / R;
        PK = PK / (1 + PK/R) + Q;
        KH = PH / (1 - theta*PH + PH/R) / R;
        PH = PH / (1 - theta*PH + PH/R) + Q;
        if ((1/PH - theta + 1/R) <= 0) || (PH <= 0)
            error('theta is too large')
        end
    end
    x = x + sqrt(Q) * randn;
    if BiasFlag 
        x = x + 10;
    end
    xhatK = xhatK + KK * (y - xhatK);
    xhatH = xhatH + KH * (y - xhatH);
    xArr(k+1) = x;
    xhatKArr(k+1) = xhatK;
    xhatHArr(k+1) = xhatH;
end

k = 0 : kf;
close all; figure; hold on;
plot(k, xArr, 'k', 'LineWidth', 2.5);
plot(k, xhatKArr, 'r--', k, xhatHArr, 'b:');
set(gca,'FontSize',12); set(gcf,'Color','White');
legend('true state', 'Kalman estimate', 'H_{\infty} estimate');
xlabel('time'); ylabel('state value');
set(gca,'box','on');

RMSK = sqrt(norm(xArr-xhatKArr)^2/(kf+1));
RMSH = sqrt(norm(xArr-xhatHArr)^2/(kf+1));
disp(['Kalman RMS estimation error = ', num2str(RMSK)]);
disp(['H-infinity RMS estimation error = ', num2str(RMSH)]);

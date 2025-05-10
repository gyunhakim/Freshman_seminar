clear all
close all
clc
% 중요: 비교할게 뭔지, 내가 뭘보여주려고 하는건지 확실히하기
% load('noise_d_sample.mat'); % noise의 영향이 아닌 다른 것의 영향을 보기 위해 랜덤값을 고정시킬때 사용 
% inital position (5.3, 3.6)
pos_ini = [5.3 3.6];
% initial state
x = [pos_ini(1) 0 pos_ini(2) 0]'; % 2,4에 속도를 추가가능, 속도값으로 필터의 성능을 알 수 있음.-> 속도 성능이 지정값에 가까우면 필터성능이 좋다고 볼 수 있음

%% plant motion
dt = 0.5; % sampling time
vx = 0.2; % x축 속도
vy = 0.1; % y축 속도
% saved measurement (sample data)
saved_pos = [pos_ini(1) 0:vx*dt:2; pos_ini(2) 5:vy*dt:6];
saved_pos_noise = saved_pos + 0.1*randn(size(saved_pos)); % 노이즈 크기 조절할 여지 존재
% matrix A
A = [1 dt 0 0; 0 1 0 0; 0 0 1 dt; 0 0 0 1]; % 다룰 값에 따라 열의 값 수정
% matrix C
C = [1 0 0 0; 0 0 1 0]; % 다룰 값에 따라 열의 값 수정

%% observability matrix
% rank([C; C*A; C*A^2; C*A^3])

%% Kalman filter initialization
P = cov(x);% initial covariance % x의 초기값이 존재하지 않을 경우엔 일반적인 값을 넣음

for k = 1:size(saved_pos,2)
    % noise
    %Q=0;
    %R=0;  
    Q = 0.1*eye(4); %A*P*A에 4x4여서 사이즈를 맞춰야 한다. 4x1이면 자동으로 열이 복사되어 4x4가 된다. 웬만하면 단위행렬 형태로.
    R = 0.1*eye(2); %C*Pk_predict*C'가 2x2행렬이라서 사이즈를 맞춰야 한다. 2x1이면 자동으로 열하나가 복사되어 2x2가 된다. 이게 작아질수록 update 성능이 올라간다.
    % Q,R 비율중요.
    % measurement
    %y = saved_pos(:,k);
    y = saved_pos_noise(:,k); % measurement값에 noise를 섞어 현실성을 부여
    %% Kalman filter
    xk_predict = A*x; % 매 루프마다 값이 갱신되어서 값 따로 볼거 아니면 사실 구별지을 필요는 없다.
    Pk_predict = A*P*A' + Q;
    Kk = Pk_predict*C'/(C*Pk_predict*C' + R);
    x = xk_predict + Kk*(y - C*xk_predict);
    P = Pk_predict - Kk*C*Pk_predict;
    
    % trajectory
    x_predict_trajectory(:,k) = xk_predict;
    x_update_trajectory(:,k) = x;
end


figure(1)
plot(saved_pos(1,:), saved_pos(2,:), '-ok', 'MarkerSize', 5, 'LineWidth', 3); hold on;
plot(saved_pos_noise(1,:), saved_pos_noise(2,:), '-*r', 'MarkerSize', 5); hold on;
plot(x_predict_trajectory(1,:), x_predict_trajectory(3,:), '-ob', 'MarkerSize', 5); hold on;
plot(x_update_trajectory(1,:), x_update_trajectory(3,:), '-+b', 'MarkerSize', 5); hold on;
grid on;
%xlim([-1 7]); ylim([3 7]);
xlabel('x(m)'); ylabel('y(m)');
legend({'measurement', 'noisy measurement','predict','update'}, 'Fontsize',20);

%% velocity
figure(2)
subplot(1,2,1)
plot(1:k, x_predict_trajectory(2,:), '-or'); hold on;
plot(1:k, vx*ones(1,k), '-*k'); hold on;
grid on;
xlabel('x(sample)'); ylabel('y(velocity)');
legend({'x velocity', 'desired x velocity'}, 'Fontsize',20);

subplot(1,2,2)
plot(1:k, x_predict_trajectory(4,:), '-or'); hold on;
plot(1:k, vy*ones(1,k), '-*k'); hold on;
grid on;
xlabel('x(sample)'); ylabel('y(velocity)');
legend({'y velocity', 'desired y velocity'}, 'Fontsize',20);

%% position
figure(3)
subplot(1,2,1)
plot(1:k, x_predict_trajectory(1,:), '-or'); hold on;
plot(1:k, saved_pos(1,:), '-*k'); hold on;
plot(1:k, saved_pos_noise(1,:), '->b'); hold on;
grid on;
xlabel('x(sample)'); ylabel('y(position)');
legend({'x position', 'desired x position', 'desired x position + noise'}, 'Fontsize',20);

subplot(1,2,2)
plot(1:k, x_predict_trajectory(3,:), '-or'); hold on;
plot(1:k, saved_pos(2,:), '-*k'); hold on;
plot(1:k, saved_pos_noise(2,:), '->b'); hold on;
grid on;
xlabel('x(sample)'); ylabel('y(position)');
legend({'y position', 'desired y position', 'desired y position + noise'}, 'Fontsize',20);

%% absolute error
%논문 쓸 때 필수적으로 넣을 파라미터
figure,
plot(1:k, abs(saved_pos(1,:)-x_predict_trajectory(1,:)), '-or'); hold on;
% plot(1:k, saved_pos(1,:), '-*k'); hold on;
% plot(1:k, saved_pos_noise(1,:), '->b'); hold on;
grid on;
xlabel('x(sample)'); ylabel('y(error)');
legend({'absolute error'}, 'Fontsize',20);
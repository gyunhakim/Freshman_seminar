clear all
close all
clc
% �߿�: ���Ұ� ����, ���� �������ַ��� �ϴ°��� Ȯ�����ϱ�
% load('noise_d_sample.mat'); % noise�� ������ �ƴ� �ٸ� ���� ������ ���� ���� �������� ������ų�� ��� 
% inital position (5.3, 3.6)
pos_ini = [5.3 3.6];
% initial state
x = [pos_ini(1) 0 pos_ini(2) 0]'; % 2,4�� �ӵ��� �߰�����, �ӵ������� ������ ������ �� �� ����.-> �ӵ� ������ �������� ������ ���ͼ����� ���ٰ� �� �� ����

%% plant motion
dt = 0.5; % sampling time
vx = 0.2; % x�� �ӵ�
vy = 0.1; % y�� �ӵ�
% saved measurement (sample data)
saved_pos = [pos_ini(1) 0:vx*dt:2; pos_ini(2) 5:vy*dt:6];
saved_pos_noise = saved_pos + 0.1*randn(size(saved_pos)); % ������ ũ�� ������ ���� ����
% matrix A
A = [1 dt 0 0; 0 1 0 0; 0 0 1 dt; 0 0 0 1]; % �ٷ� ���� ���� ���� �� ����
% matrix C
C = [1 0 0 0; 0 0 1 0]; % �ٷ� ���� ���� ���� �� ����

%% observability matrix
% rank([C; C*A; C*A^2; C*A^3])

%% Kalman filter initialization
P = cov(x);% initial covariance % x�� �ʱⰪ�� �������� ���� ��쿣 �Ϲ����� ���� ����

for k = 1:size(saved_pos,2)
    % noise
    %Q=0;
    %R=0;  
    Q = 0.1*eye(4); %A*P*A�� 4x4���� ����� ����� �Ѵ�. 4x1�̸� �ڵ����� ���� ����Ǿ� 4x4�� �ȴ�. �����ϸ� ������� ���·�.
    R = 0.1*eye(2); %C*Pk_predict*C'�� 2x2����̶� ����� ����� �Ѵ�. 2x1�̸� �ڵ����� ���ϳ��� ����Ǿ� 2x2�� �ȴ�. �̰� �۾������� update ������ �ö󰣴�.
    % Q,R �����߿�.
    % measurement
    %y = saved_pos(:,k);
    y = saved_pos_noise(:,k); % measurement���� noise�� ���� ���Ǽ��� �ο�
    %% Kalman filter
    xk_predict = A*x; % �� �������� ���� ���ŵǾ �� ���� ���� �ƴϸ� ��� �������� �ʿ�� ����.
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
%�� �� �� �ʼ������� ���� �Ķ����
figure,
plot(1:k, abs(saved_pos(1,:)-x_predict_trajectory(1,:)), '-or'); hold on;
% plot(1:k, saved_pos(1,:), '-*k'); hold on;
% plot(1:k, saved_pos_noise(1,:), '->b'); hold on;
grid on;
xlabel('x(sample)'); ylabel('y(error)');
legend({'absolute error'}, 'Fontsize',20);
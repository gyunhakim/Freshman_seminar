clear all
close all
clc


%% plant motion
dt = 0.1; % sampling time
syms x1 y1 seta v w A B  % ��ȣ��������
% initial state
% x = [x1 y1 seta]';
% control input
% u = [v w]';
% matrix A
A = [1 0 0; 0 1 0; 0 0 1];
B = [2*cos(seta)*dt 0; sin(seta)*dt 0; 0 1*dt];
% matrix C
C = eye(3);

%% Jacobian matrix
Jx = @(a1,a2) [ 1 0 -(conj(a2)*sin(a1)*dt); ...
                0 1 (cos(a1)*conj(a2)*dt); ...
                0 0 1];
% jacobian(A*x+B*u, [x1, y1, seta]); %xhat = f(x,u,w,t)
% Jy = @(x) C*x; %y = h(x,v,t)
% Jy = @(x) jacobian(C*x, [x1, y1, seta]); %y = h(x,v,t)
%% observability matrix
% rank([C; C*A; C*A^2; C*A^3])

%% Kalman filter initialization
% inital position
x_temp = [0 0 0]';
% x_temp = subs(x,[x1, y1, seta], [pos_ini(1), pos_ini(2), pos_ini(3)]); % x�� ���ں���� ���� ���ڿ�, x_temp�� ���� ���꿡 ����� ���ڿ�, coveriance ������ ���� ����
P = cov(x_temp); % initial covariance % x�� �ʱⰪ�� �������� ���� ��쿣 �Ϲ����� ���� ����
%u_temp = subs(u,[v w],[2, 0.25]); % u�� ���ں���� ���� ���ڿ�, u_temp�� ���� ���꿡 ����� ���ڿ�
sim_time = 60;

for k = 1:sim_time % �� �������� block diagram�� ���ư��ٰ� �����ϱ�    
%     Jx1 = Jx(x,u); 
%     Jy1 = Jy(x);
%     x_temp = subs(x,[x1, y1, seta], [pos_ini(1), pos_ini(2), (k-1)*pi/180]); % x�� ���ں���� ���� ���ڿ�, x_temp�� ���� ���꿡 ����� ���ڿ�
%     u_temp = subs(u,[v w],[sin(k/5), 0.01]); % u�� ���ں���� ���� ���ڿ�, u_temp�� ���� ���꿡 ����� ���ڿ�
%     Jx1_temp = subs(Jx1,[v, sin(seta), cos(seta)],[u_temp(1), sin(x_temp(3)), cos(x_temp(3))]); %seta���ڿ��� �Լ��ȿ� �������� �Լ���°�� �����ؾ� ��ȯ��
    
    %% Car Model
    if k == 1
        car_pos = x_temp;
    end
    
    car_input = [1 0.01]'; % �տ��� ������ ���ڿ� u��
    car_pos = Car_Model(car_pos, car_input, dt); %�ý����� ���ĳ��� ������ ����
    y = car_pos + randn(size(car_pos))*0.01; % x,y,seta�� x,y
    y_real = car_pos;
    saved_pos_real(:,k) = y_real; % x,y,seta ����
    saved_pos(:,k) = y;
    Jx1 = Jx(car_pos(3),car_input(1)); 
    Jy1 = C;        
    % noise,     % Q,R �����߿�.
    Q = 0.4*eye(3); %A*P*A�� 4x4���� ����� ����� �Ѵ�. 4x1�̸� �ڵ����� ���� ����Ǿ� 4x4�� �ȴ�. �����ϸ� ������� ���·�.
    R = 0.1*eye(3); %C*Pk_predict*C'�� 2x2����̶� ����� ����� �Ѵ�. 2x1�̸� �ڵ����� ���ϳ��� ����Ǿ� 2x2�� �ȴ�. �̰� �۾������� update ������ �ö󰣴�.
    
    %% Extended Kalman filter
    xk_predict = Jx1*x_temp; % �� �������� ���� ���ŵǾ �� ���� ���� �ƴϸ� ��� �������� �ʿ�� ����.
    Pk_predict = Jx1*P*Jx1' + Q;
    Kk = Pk_predict*Jy1'*(Jy1*Pk_predict*Jy1' + R)^(-1);
    x_temp = xk_predict + Kk*(y - Jy1*(xk_predict));
    P = Pk_predict - Kk*Jy1*Pk_predict;
    
    % trajectory
    x_predict_trajectory(:,k) = xk_predict;
    x_update_trajectory(:,k) = x_temp;
end

figure(1)
plot(saved_pos(1,:), saved_pos(2,:), '-ob', 'MarkerSize', 5); hold on;
plot(saved_pos_real(1,:), saved_pos_real(2,:), '-og', 'MarkerSize', 5); hold on;
plot(x_update_trajectory(1,:), x_update_trajectory(2,:), '-+r', 'MarkerSize', 5); hold on;
grid on;
xlabel('x(m)'); ylabel('y(m)');
legend({'measurement','real','update'}, 'Fontsize',20);
%% velocity
%{
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
%}
%% position
figure(3)
subplot(1,2,1)
plot(1:k, x_update_trajectory(1,:), '-or'); hold on;
plot(1:k, saved_pos(1,:), '-*k'); hold on;
grid on;
xlabel('x(sample)'); ylabel('y(position)');
legend({'x position', 'desired x position'}, 'Fontsize',20);
subplot(1,2,2)
plot(1:k, x_update_trajectory(2,:), '-or'); hold on;
plot(1:k, saved_pos(2,:), '-*k'); hold on;
grid on;
xlabel('x(sample)'); ylabel('y(position)');
legend({'y position', 'desired y position'}, 'Fontsize',20);
%% absolute error
%�� �� �� �ʼ������� ���� �Ķ����
figure,
plot(1:k, abs(saved_pos(1,:)-x_update_trajectory(1,:)), '-or'); hold on;
grid on;
xlabel('x(sample)'); ylabel('y(error)');
legend({'absolute error'}, 'Fontsize',20);
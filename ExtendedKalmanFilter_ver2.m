clear all
close all
clc


%% plant motion
dt = 0.1; % sampling time
syms x1 y1 seta v w A B  % 기호변수선언
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
% x_temp = subs(x,[x1, y1, seta], [pos_ini(1), pos_ini(2), pos_ini(3)]); % x는 자코비안을 위한 문자열, x_temp는 실제 연산에 사용할 숫자열, coveriance 연산을 위한 정의
P = cov(x_temp); % initial covariance % x의 초기값이 존재하지 않을 경우엔 일반적인 값을 넣음
%u_temp = subs(u,[v w],[2, 0.25]); % u는 자코비안을 위한 문자열, u_temp는 실제 연산에 사용할 숫자열
sim_time = 60;

for k = 1:sim_time % 이 루프동안 block diagram이 돌아간다고 생각하기    
%     Jx1 = Jx(x,u); 
%     Jy1 = Jy(x);
%     x_temp = subs(x,[x1, y1, seta], [pos_ini(1), pos_ini(2), (k-1)*pi/180]); % x는 자코비안을 위한 문자열, x_temp는 실제 연산에 사용할 숫자열
%     u_temp = subs(u,[v w],[sin(k/5), 0.01]); % u는 자코비안을 위한 문자열, u_temp는 실제 연산에 사용할 숫자열
%     Jx1_temp = subs(Jx1,[v, sin(seta), cos(seta)],[u_temp(1), sin(x_temp(3)), cos(x_temp(3))]); %seta문자열이 함수안에 들어가있으면 함수통째로 선언해야 변환됨
    
    %% Car Model
    if k == 1
        car_pos = x_temp;
    end
    
    car_input = [1 0.01]'; % 앞에서 선언한 문자열 u값
    car_pos = Car_Model(car_pos, car_input, dt); %시스템을 거쳐나온 측정값 생성
    y = car_pos + randn(size(car_pos))*0.01; % x,y,seta중 x,y
    y_real = car_pos;
    saved_pos_real(:,k) = y_real; % x,y,seta 축적
    saved_pos(:,k) = y;
    Jx1 = Jx(car_pos(3),car_input(1)); 
    Jy1 = C;        
    % noise,     % Q,R 비율중요.
    Q = 0.4*eye(3); %A*P*A에 4x4여서 사이즈를 맞춰야 한다. 4x1이면 자동으로 열이 복사되어 4x4가 된다. 웬만하면 단위행렬 형태로.
    R = 0.1*eye(3); %C*Pk_predict*C'가 2x2행렬이라서 사이즈를 맞춰야 한다. 2x1이면 자동으로 열하나가 복사되어 2x2가 된다. 이게 작아질수록 update 성능이 올라간다.
    
    %% Extended Kalman filter
    xk_predict = Jx1*x_temp; % 매 루프마다 값이 갱신되어서 값 따로 볼거 아니면 사실 구별지을 필요는 없다.
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
%논문 쓸 때 필수적으로 넣을 파라미터
figure,
plot(1:k, abs(saved_pos(1,:)-x_update_trajectory(1,:)), '-or'); hold on;
grid on;
xlabel('x(sample)'); ylabel('y(error)');
legend({'absolute error'}, 'Fontsize',20);
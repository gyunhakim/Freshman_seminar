clear all
close all
clc

%% 
% x position y position width height
Tunnel = [5,10,5,20]; % 터널 사이즈 지정


car = figure('name','car');
ax_sim = axes(car);
% 異? limit
xlim(ax_sim,[-10 50]); %구역 지정
ylim(ax_sim,[-10 50]);
TunnelArea = rectangle(ax_sim, 'position', Tunnel,'facecolor',[0.5 0.5 0.5]);
% car definition
[car1_front1, car1_front2, car1_rear1, car1_rear2, car1_chassis, car1_center, car1_measurement] = Set_func(ax_sim);
[car2_front1, car2_front2, car2_rear1, car2_rear2, car2_chassis, car2_center, car2_measurement] = Set_func(ax_sim);
[car3_front1, car3_front2, car3_rear1, car3_rear2, car3_chassis, car3_center, car3_measurement] = Set_func(ax_sim);

%% initial pose : x position y position steering angle(deg ? rad?) angular velocity
car1_pose = [0 20 0 0];      car1_trajectory = car1_pose; car1_tra_mea = car1_pose;
car2_pose = [0 8 0 0];      car2_trajectory = car2_pose; car2_tra_mea = car2_pose;
car3_pose = [-1 10 0 0];      car3_trajectory = car3_pose; car3_tra_mea = car3_pose; % 초기 지점

dt = 0.1; % 샘플링 간격
car1_length = 2;
car2_length = 2;
car3_length = 2;
%% Kalman filter initialization
car1_est = zeros(4,1); P_car1 = cov(car1_est);
car2_est = zeros(4,1); P_car2 = cov(car2_est);
car3_est = zeros(4,1); P_car3 = cov(car3_est);
car1_est2 = zeros(4,1); P_car12 = cov(car1_est);
car2_est2 = zeros(4,1); P_car22 = cov(car2_est);
car3_est2 = zeros(4,1); P_car32 = cov(car3_est);

sim_t = 60;

for k = 1:sim_t % 샘플링 횟수
    %% measurement 
    
    car1_mea = GetMeasurement(Tunnel, car1_pose);
    car2_mea = GetMeasurement(Tunnel, car2_pose);
    car3_mea = GetMeasurement(Tunnel, car3_pose);

    
    %% filter Kalman
    [car1_est, P_car1] = KalmanFilter(car1_mea,dt,car1_est, P_car1);
    [car2_est, P_car2] = KalmanFilter(car2_mea,dt,car2_est, P_car2);
    [car3_est, P_car3] = KalmanFilter(car3_mea,dt,car3_est, P_car3);
    
    %% FIR filter
    [car1_est2] = FIRFilter(car1_mea,dt,car1_est2,k);
    [car2_est2] = FIRFilter2(car2_mea,dt,car2_est2,k);
    [car3_est2] = FIRFilter3(car3_mea,dt,car3_est2,k);
    
    %% path planning
    
    car1_ref = [20 0];
    car1_error = car1_pose(1:2) - car1_ref;
    car2_ref = [50 car2_pose(2)]; 
    car2_error = car2_ref - car2_pose(1:2);
    %% controller
    % velocity % angular velocity, 속도제어기
    car1_u = [4 0.05];
%     car1_u(1) = car1_ref(1) - car1_mea(1);
%     car1_u(2) = car1_ref(2) - car1_mea(2);
    car2_u = [car2_error(1)*0.5 0];
    car3_u = [5 0.01*sin(k/100)];
    
    %% Car movement from control input
    [car1_trajectory, car1_pose, car1_tra_mea] = Car_model(car1_front1,car1_front2,car1_rear1,car1_rear2,car1_chassis,car1_center,car1_measurement,car1_pose,car1_trajectory,car1_mea,car1_tra_mea,dt,car1_u,car1_length,k,ax_sim);
    [car2_trajectory, car2_pose, car2_tra_mea] = Car_model(car2_front1,car2_front2,car2_rear1,car2_rear2,car2_chassis,car2_center,car2_measurement,car2_pose,car2_trajectory,car2_mea,car2_tra_mea,dt,car2_u,car2_length,k,ax_sim);
    [car3_trajectory, car3_pose, car3_tra_mea] = Car_model(car3_front1,car3_front2,car3_rear1,car3_rear2,car3_chassis,car3_center,car3_measurement,car3_pose,car3_trajectory,car3_mea,car3_tra_mea,dt,car3_u,car3_length,k,ax_sim);
    
    title(ax_sim,['t = ', num2str(k)]);
    pause(dt);
    
    save_car1_measurement(k,:) = car1_mea;
    save_car2_measurement(k,:) = car2_mea;
    save_car3_measurement(k,:) = car3_mea;
    
    save_car1_kf(k,:) = car1_est;
    save_car2_kf(k,:) = car2_est;
    save_car3_kf(k,:) = car3_est;
    save_car1_fir(k,:) = car1_est2;
    save_car2_fir(k,:) = car2_est2;
    save_car3_fir(k,:) = car3_est2;
    
    if k == 1
        car1_velocity_trajectory(k,:) = zeros(4,1);
        car2_velocity_trajectory(k,:) = zeros(4,1);
        car3_velocity_trajectory(k,:) = zeros(4,1);
    else
        car1_velocity_trajectory(k,:) = (car1_trajectory(k,:) - car1_trajectory(k-1,:))/dt;
        car2_velocity_trajectory(k,:) = (car2_trajectory(k,:) - car2_trajectory(k-1,:))/dt;
        car3_velocity_trajectory(k,:) = (car3_trajectory(k,:) - car3_trajectory(k-1,:))/dt;
    end
end

%% absolute error
%논문 쓸 때 필수적으로 넣을 파라미터
% figure,
% plot(1:k, abs(car1_tra_mea(:,1)-car1_trajectory(:,1)), '-or'); hold on;
% plot(1:k, abs(car2_tra_mea(:,1)-car2_trajectory(:,1)), '-or'); hold on;
% plot(1:k, abs(car3_tra_mea(:,1)-car3_trajectory(:,1)), '-or'); hold on;
%  plot(1:k, saved_pos(1,:), '-*k'); hold on;
%  plot(1:k, saved_pos_noise(1,:), '->b'); hold on;
% grid on;
% xlabel('x(sample)'); ylabel('y(error)');
% legend({'car1absolute error','car2absolute error','car3absolute error'}, 'Fontsize',20);


figure,
plot(car1_trajectory(:,1), car1_trajectory(:,2), '-or'); hold on;
plot(save_car1_measurement(:,1), save_car1_measurement(:,2), '-ob'); hold on;
plot(save_car1_kf(:,1), save_car1_kf(:,3), '-og'); hold on;
plot(save_car1_fir(:,1), save_car1_fir(:,3), '-om'); hold on;
grid on;
xlabel('x position'); ylabel('y position');
legend('real', 'measurement','Kalman','FIR');
title('car1');
figure,
plot(car2_trajectory(:,1), car2_trajectory(:,2), '-or'); hold on;
plot(save_car2_measurement(:,1), save_car2_measurement(:,2), '-ob'); hold on;
plot(save_car2_kf(:,1), save_car2_kf(:,3), '-og'); hold on;
plot(save_car2_fir(:,1), save_car2_fir(:,3), '-om'); hold on;
grid on;
xlabel('x position'); ylabel('y position');
legend('real', 'measurement','Kalman','FIR');
title('car2');
figure,
plot(car3_trajectory(:,1), car3_trajectory(:,2), '-or'); hold on;
plot(save_car3_measurement(:,1), save_car3_measurement(:,2), '-ob'); hold on;
plot(save_car3_kf(:,1), save_car3_kf(:,3), '-og'); hold on;
plot(save_car3_fir(:,1), save_car3_fir(:,3), '-om'); hold on;
grid on;
xlabel('x position'); ylabel('y position');
legend('real', 'measurement','Kalman','FIR');
title('car3');


figure,
subplot(3,2,1)
plot(1:sim_t, car1_velocity_trajectory(:,1), '-or'); hold on;
plot(1:sim_t, save_car1_kf(:,2), '-og'); hold on;
plot(1:sim_t, save_car1_fir(:,2), '-ob'); hold on;
grid on;
xlabel('simulation time'); ylabel('car1 x velocity(:,1)');
legend('real','Kalman','FIR');

subplot(3,2,2)
plot(1:sim_t, car1_velocity_trajectory(:,2), '-or'); hold on;
plot(1:sim_t, save_car1_kf(:,4), '-og'); hold on;
plot(1:sim_t, save_car1_fir(:,4), '-ob'); hold on;
grid on;
xlabel('simulation time'); ylabel('car1 y velocity(:,2)');
legend('real','Kalman','FIR');

subplot(3,2,3)
plot(1:sim_t, car2_velocity_trajectory(:,1), '-or'); hold on;
plot(1:sim_t, save_car2_kf(:,2), '-og'); hold on;
plot(1:sim_t, save_car2_fir(:,2), '-ob'); hold on;
grid on;
xlabel('simulation time'); ylabel('car2 x velocity(:,1)');
legend('real','Kalman','FIR');

subplot(3,2,4)
plot(1:sim_t, car2_velocity_trajectory(:,2), '-or'); hold on;
plot(1:sim_t, save_car2_kf(:,4), '-og'); hold on;
plot(1:sim_t, save_car2_fir(:,4), '-ob'); hold on;
grid on;
xlabel('simulation time'); ylabel('car2 x velocity(:,2)');
legend('real','Kalman','FIR');

subplot(3,2,5)
plot(1:sim_t, car3_velocity_trajectory(:,1), '-or'); hold on;
plot(1:sim_t, save_car3_kf(:,2), '-og'); hold on;
plot(1:sim_t, save_car3_fir(:,2), '-ob'); hold on;
grid on;
xlabel('simulation time'); ylabel('car3 x velocity(:,1)');
legend('real','Kalman','FIR');

subplot(3,2,6)
plot(1:sim_t, car3_velocity_trajectory(:,2), '-or'); hold on;
plot(1:sim_t, save_car3_kf(:,4), '-og'); hold on;
plot(1:sim_t, save_car3_fir(:,4), '-ob'); hold on;
grid on;
xlabel('simulation time'); ylabel('car3 x velocity(:,2)');
legend('real','Kalman','FIR');
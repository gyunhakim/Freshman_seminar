clear; clc; close all;

% [1] 시스템/측정 행렬 정의
A = [0.8  0.1;
     0    0.9];
C = [1 0];
n = size(A,1);
m = size(C,1);

% [2] 파라미터 설정
N = 4;                            % N 샘플
Omega = eye(m*N);                 % (m*N x m*N) ; 여기서는 4x4 단위행렬

[G, ThetaN] = FIR_gain(A, C, Omega, N);

disp('Filter gain G =');
disp(G);

% [3] 시뮬레이션 설정
T = 20;                  % 시뮬레이션 총 스텝
x_true = zeros(n, T+1);  % 실제 상태 (k=0 ~ k=T)
z_meas = zeros(m, T+1);  % 실제 측정
x_hat = zeros(n, T+1);   % 추정값 (DTFWNUF)

% 초기값
x0 = [1; -1];
x_true(:,1) = x0;
z_meas(:,1) = C*x0;

% [4] 메인 루프
for k = 1:T
    % -- 실제 시스템 진행
    w_k = [0; 0];       % 프로세스 잡음 (여기선 0으로)
    x_true(:,k+1) = A * x_true(:,k) + w_k;

    % -- 측정
    v_k = 0;            % 측정 잡음
    z_meas(:,k+1) = C * x_true(:,k+1) + v_k;

    % -- FIR 필터 적용 (k >= N일 때)
    if k >= N
        % Z_{k-1} = [z_meas(k); z_meas(k-1); ...; z_meas(k-N+1)]^T
        %           (m*N x 1)  = (1 x 4)^T 여기서는 4x1
        Z_k_1 = [];
        for j = 0:N-1
            Z_k_1 = [ Z_k_1; z_meas(:, k-j) ];  
        end

        % x_hat(k) = G * Z_{k-1}
        x_hat(:, k+1) = G * Z_k_1;  
    else
        x_hat(:, k+1) = x_hat(:, k);
    end
end

% [5] 결과
t = 0:T;
figure; 
subplot(2,1,1);
plot(t, x_true(1,:), 'o-k','LineWidth',1.5); hold on;
plot(t, x_hat(1,:),  'x--r','LineWidth',1.5);
legend('true x_1','hat x_1','Location','best');
grid on; xlabel('time step k'); ylabel('x(1)');

subplot(2,1,2);
plot(t, x_true(2,:), 'o-k','LineWidth',1.5); hold on;
plot(t, x_hat(2,:),  'x--r','LineWidth',1.5);
legend('true x_2','hat x_2','Location','best');
grid on; xlabel('time step k'); ylabel('x(2)');

sgtitle('FIR Filter Example');

%% FIR Gain
function [G, ThetaN] = FIR_gain(A, C, Omega, N)
    [n, ~] = size(A);
    [m, ~] = size(C);

    ThetaN = zeros(m*N, n);
    for i = 1:N
        Ai = A^(N - i);      % A^(N-i)
        ThetaN( (i-1)*m+1 : i*m, : ) = C * Ai;
    end

    temp = ThetaN' * (Omega^2) * ThetaN;   % (n x n)
    G = A^N * ( temp \ ( ThetaN' * (Omega^2) ) );  % inv(temp)*(...) 꼴

end

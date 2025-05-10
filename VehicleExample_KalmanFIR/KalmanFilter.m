function [car_hat, P] = KalmanFilter(car_pose,dt,x,P) % 내가 해야할 부분
   %% plant motion
    % matrix A
    A1 = [1 dt; 0 1]; 
    A = blkdiag(A1,A1);
    % matrix C
    C1 = [1 0];
    C = blkdiag(C1, C1);
    % noise
    Q = 1.5*eye(4); %A*P*A에 4x4여서 사이즈를 맞춰야 한다. 4x1이면 자동으로 열이 복사되어 4x4가 된다. 웬만하면 단위행렬 형태로.
    R = 0.1*eye(2); %C*Pk_predict*C'가 2x2행렬이라서 사이즈를 맞춰야 한다. 2x1이면 자동으로 열하나가 복사되어 2x2가 된다. 이게 작아질수록 update 성능이 올라간다.
    % Q,R 비율중요.

    % saved measurement (sample data)
    y = [car_pose(1) car_pose(2)]';

    %% Kalman filter
    if car_pose(1) == 0 && car_pose(2) == 0
        x = A*x;
    else
        xk_predict = A*x; % 매 루프마다 값이 갱신되어서 값 따로 볼거 아니면 사실 구별지을 필요는 없다.
        Pk_predict = A*P*A' + Q;
        Kk = Pk_predict*C'/(C*Pk_predict*C' + R);
        x = xk_predict + Kk*(y - C*xk_predict);
        P = Pk_predict - Kk*C*Pk_predict;
    end

    car_hat=x;
    
end
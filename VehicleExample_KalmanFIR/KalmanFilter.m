function [car_hat, P] = KalmanFilter(car_pose,dt,x,P) % ���� �ؾ��� �κ�
   %% plant motion
    % matrix A
    A1 = [1 dt; 0 1]; 
    A = blkdiag(A1,A1);
    % matrix C
    C1 = [1 0];
    C = blkdiag(C1, C1);
    % noise
    Q = 1.5*eye(4); %A*P*A�� 4x4���� ����� ����� �Ѵ�. 4x1�̸� �ڵ����� ���� ����Ǿ� 4x4�� �ȴ�. �����ϸ� ������� ���·�.
    R = 0.1*eye(2); %C*Pk_predict*C'�� 2x2����̶� ����� ����� �Ѵ�. 2x1�̸� �ڵ����� ���ϳ��� ����Ǿ� 2x2�� �ȴ�. �̰� �۾������� update ������ �ö󰣴�.
    % Q,R �����߿�.

    % saved measurement (sample data)
    y = [car_pose(1) car_pose(2)]';

    %% Kalman filter
    if car_pose(1) == 0 && car_pose(2) == 0
        x = A*x;
    else
        xk_predict = A*x; % �� �������� ���� ���ŵǾ �� ���� ���� �ƴϸ� ��� �������� �ʿ�� ����.
        Pk_predict = A*P*A' + Q;
        Kk = Pk_predict*C'/(C*Pk_predict*C' + R);
        x = xk_predict + Kk*(y - C*xk_predict);
        P = Pk_predict - Kk*C*Pk_predict;
    end

    car_hat=x;
    
end
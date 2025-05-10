function [car_hat] = FIRFilter(car_pose,dt,x,k) % ���� �ؾ��� �κ�
%% global ������ ��ȯ
persistent F H Y
persistent firstrun
if isempty(firstrun)
    F = [];
    H = [];
    Y = [];
    firstrun = 0;
end
%���� iteraation����...(��Ȳ���� iteration����)
%% plant motion
% matrix A
A1 = [1 dt; 0 1];
A = blkdiag(A1,A1);
% matrix C
C1 = [1 0];
C = blkdiag(C1, C1);
y = [car_pose(1) car_pose(2)]';

%% Extended state space equation
N=6; % horizon size
K=2; % Batch size
m = k-N+1; % Batch start
s = k-N+K; % Batch end

if K>=k % ó�� N sample���� batch�����ŭ data stack
    F(1+size(A,1)*(k-1):4+size(A,1)*(k-1),:) = A*A^(k-1);
    H(1+size(y,1)*(k-1):2+size(y,1)*(k-1),:) = C*A^(k-1);
end
    
%H3 = [C; C*A];
% Y3(1+2*(k-1):2+2*(k-1),:) = y;

%% FIR filter
if car_pose(1) == 0 && car_pose(2) == 0 % measurement���� �������ʴ� ��Ȳ�� �� predict���� �־��ֱ� ����.
    x = A*x;
else
    if m > 0    % ó�� N sample�� ������ batch ���� ��
        if m == 1
            Y(1+size(y,1)*(k-1):2+size(y,1)*(k-1),:) = y; % if�� else�� �ڵ带 ���������ʾƼ� ���� �����������.
        else
            Y = [Y(1+size(y,1):N*size(y,1),:); y]; % ���� Y���� ù ������ ������ ���õ�� ���ο� ���÷� �籸��
        end
        % Batch
        G = inv(H'*H);
        x = G*H'*Y(1:K*2,:);
        if m > K
            %Iteration
            x = A*x;
            G = inv(C'*C+inv(A*G*A'));
            ki = G*C';
            x = x + ki*(y-C*x);

            %Pk_predict = A*P*A' + Q;
            %P = (eye() - Kk*C)*Pk_predict*(eye() - Kk*C)' + ;
        end
    else % batch ���� ������
        Y(1+size(y,1)*(k-1):2+size(y,1)*(k-1),:) = y;
    end
end

disp(k), disp(Y), disp(H), disp(F)
car_hat = x;

end
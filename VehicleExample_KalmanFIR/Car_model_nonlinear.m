% function [Front_wheel_1,Front_wheel_2,Rear_wheel_1,Rear_wheel_2,chassis,car_pose] ...
    function [car_trajectory, car_pose, car_tra_mea] ...
    = Car_model_nonlinear(car_front1,car_front2,car_rear1,car_rear2,car_chassis,car_center,car_measurement,car_pose,car_trajectory,car_mea,car_tra_mea,dt,u,car_length,k,ax)
    car_x_dot = u(1) * cos(car_pose(3));
    car_y_dot = u(1) * sin(car_pose(3));
    car_theta_dot = u(1) / car_length *tan(car_pose(4));
    car_phi_dot = u(2);
    
    car_pose_dot = [car_x_dot car_y_dot car_theta_dot car_phi_dot];
    
    car_pose = car_pose + dt*car_pose_dot;
    
    lim_steering = [-pi/3 pi/3];
    if car_pose(4) > lim_steering(2)
        car_pose(4) = lim_steering(2);
    elseif car_pose(4) < lim_steering(1)
        car_pose(4) = lim_steering(1);
    end
    
    
    r = car_length / 2;
    chassis = [
        car_pose(1) + car_length/2*cos(car_pose(3)-pi/2)                                , car_pose(2) + car_length/2*sin(car_pose(3)-pi/2)                                  ;
        car_pose(1) + car_length/2*cos(car_pose(3)+pi/2)                                , car_pose(2) + car_length/2*sin(car_pose(3)+pi/2)                                  ;
        car_pose(1)                                                                     , car_pose(2)                                                                       ;
        car_pose(1) + car_length*cos(car_pose(3))                                       , car_pose(2) + car_length*sin(car_pose(3))                                         ;
        car_pose(1) + car_length*cos(car_pose(3)) + car_length/2*cos(car_pose(3)-pi/2)  , car_pose(2) + car_length*sin(car_pose(3)) + car_length/2*sin(car_pose(3)-pi/2)    ;
        car_pose(1) + car_length*cos(car_pose(3)) + car_length/2*cos(car_pose(3)+pi/2)  , car_pose(2) + car_length*sin(car_pose(3)) + car_length/2*sin(car_pose(3)+pi/2)    ;
        car_pose(1) + car_length*cos(car_pose(3))                                       , car_pose(2) + car_length*sin(car_pose(3))
        ];
    
    Front_wheel_1 = [
        car_pose(1) + car_length*cos(car_pose(3)) + car_length/2*cos(car_pose(3)-pi/2) + r*cos(car_pose(3)+car_pose(4))  , car_pose(2) + car_length*sin(car_pose(3)) + car_length/2*sin(car_pose(3)-pi/2) + r*sin(car_pose(3)+car_pose(4)) ;
        car_pose(1) + car_length*cos(car_pose(3)) + car_length/2*cos(car_pose(3)-pi/2) - r*cos(car_pose(3)+car_pose(4))  , car_pose(2) + car_length*sin(car_pose(3)) + car_length/2*sin(car_pose(3)-pi/2) - r*sin(car_pose(3)+car_pose(4))
        ];
    
    Front_wheel_2 = [
        car_pose(1) + car_length*cos(car_pose(3)) + car_length/2*cos(car_pose(3)+pi/2) + r*cos(car_pose(3)+car_pose(4))  , car_pose(2) + car_length*sin(car_pose(3)) + car_length/2*sin(car_pose(3)+pi/2) + r*sin(car_pose(3)+car_pose(4)) ;
        car_pose(1) + car_length*cos(car_pose(3)) + car_length/2*cos(car_pose(3)+pi/2) - r*cos(car_pose(3)+car_pose(4))  , car_pose(2) + car_length*sin(car_pose(3)) + car_length/2*sin(car_pose(3)+pi/2) - r*sin(car_pose(3)+car_pose(4))
        ];
    
    Rear_wheel_1 = [
        car_pose(1) + car_length/2*cos(car_pose(3)-pi/2) + r*cos(car_pose(3))  , car_pose(2) + car_length/2*sin(car_pose(3)-pi/2) + r*sin(car_pose(3)) ;
        car_pose(1) + car_length/2*cos(car_pose(3)-pi/2) - r*cos(car_pose(3))  , car_pose(2) + car_length/2*sin(car_pose(3)-pi/2) - r*sin(car_pose(3))
        ];
    
    Rear_wheel_2 = [
        car_pose(1) + car_length/2*cos(car_pose(3)+pi/2) + r*cos(car_pose(3))  , car_pose(2) + car_length/2*sin(car_pose(3)+pi/2) + r*sin(car_pose(3)) ;
        car_pose(1) + car_length/2*cos(car_pose(3)+pi/2) - r*cos(car_pose(3))  , car_pose(2) + car_length/2*sin(car_pose(3)+pi/2) - r*sin(car_pose(3))
        ];
    
    car_trajectory(k,:) = car_pose;
    car_tra_mea(k,:) = car_mea;
    car_front1.XData = Front_wheel_1(:,1);
    car_front1.YData = Front_wheel_1(:,2);
    car_front2.XData = Front_wheel_2(:,1);
    car_front2.YData = Front_wheel_2(:,2);
    car_rear1.XData = Rear_wheel_1(:,1);
    car_rear1.YData = Rear_wheel_1(:,2);
    car_rear2.XData = Rear_wheel_2(:,1);
    car_rear2.YData = Rear_wheel_2(:,2);
    car_chassis.XData = chassis(:,1);
    car_chassis.YData = chassis(:,2);
    car_center.XData = car_trajectory(:,1);
    car_center.YData = car_trajectory(:,2);
    car_measurement.XData = car_tra_mea(:,1);
    car_measurement.YData = car_tra_mea(:,2);
    title(ax,['t = ', num2str(k)]);
    % pause(dt);
    end


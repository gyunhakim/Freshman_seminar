function [car_mea] = GetMeasurement(Tunnel,car_pose) 

area_x_min = Tunnel(1);
area_x_max = Tunnel(1) + Tunnel(3);
area_y_min = Tunnel(2);
area_y_max = Tunnel(2) + Tunnel(4);
car_mea = car_pose;

if car_pose(1) >= area_x_min && car_pose(1) <= area_x_max && car_pose(2) >= area_y_min && car_pose(2) <= area_y_max
    car_mea(1) = 0;
    car_mea(2) = 0;
else
    %% noise 
    car_mea = car_pose + 0.1*rand(size(car_pose));
end

end

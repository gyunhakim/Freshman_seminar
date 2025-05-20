function car_measurement = Car_Model(car_pos, car_input, dt)
car_measurement = [car_pos(1) + dt*car_input(1)*cos(car_pos(3)); ...    %car matrix(car_pos에 위치값, car_input에 속도값, dt에 샘플타임)
                    car_pos(2) + dt*car_input(1)*sin(car_pos(3)); ...
                    car_pos(3) + dt*car_input(2)  ];
end

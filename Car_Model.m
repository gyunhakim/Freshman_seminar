function car_measurement = Car_Model(car_pos, car_input, dt)
car_measurement = [car_pos(1) + dt*car_input(1)*cos(car_pos(3)); ...    %car matrix(car_pos�� ��ġ��, car_input�� �ӵ���, dt�� ����Ÿ��)
                    car_pos(2) + dt*car_input(1)*sin(car_pos(3)); ...
                    car_pos(3) + dt*car_input(2)  ];
end

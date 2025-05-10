function [car_front1, car_front2, car_rear1, car_rear2, car_chassis, car_center, car_measurement]= Set_func(ax, Tunnel)
hold(ax,'on');
car_front1 = plot(ax,0,0,'linewidth',3, 'color','red');
car_front2 = plot(ax,0,0,'linewidth',3, 'color','red');
car_rear1 = plot(ax,0,0,'linewidth',3, 'color',[1, 0.84 0]);
car_rear2 = plot(ax,0,0,'linewidth',3, 'color',[1, 0.84 0]);
car_chassis = plot(ax,0,0,'linewidth',3, 'color','blue');
car_center = plot(ax,0,0,'linewidth',1, 'color','black','Marker','+');
car_measurement = plot(ax,0,0,'linewidth',1, 'color','red','Marker','o');
% TunnelArea = rectangle(ax, 'position', Tunnel,'facecolor',[0.5 0.5 0.5]);



legend(ax,[car_center car_measurement], 'Actual pose', 'Measured pose');
hold(ax,'off');
end


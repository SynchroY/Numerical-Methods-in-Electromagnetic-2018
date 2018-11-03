clearvars;
a_m = 0.1;
h_m = 0.3;
phih_V = 100;
phi0_V = 0;

step_a = 0.01;
step_h = 0.01;

length_a = a_m/step_a+1;
length_h = h_m/step_h+1;

phi = zeros(length_a, length_h, 2);

for i=1:length_a
    phi(i,length_h,1) = phih_V;
end

times = 0;
maxerror_V = 1000;
tolerance_V = 0.0001;
this_time = 1;

while maxerror_V > tolerance_V
    times = times + 1;
    last_time = this_time;
    this_time = mod(times,2)+1;
    maxerror_V = 0;
    for i = 1:length_a
        phi(i,1,this_time) = phi0_V;
        phi(i,length_h,this_time) = phih_V;
    end
    for k = 2:length_h-1
        phi(1,k,this_time) = (2*phi(2,k,last_time)/step_a^2+(phi(1,k+1,last_time)+phi(1,k-1,this_time))/step_h^2)/(2/step_a^2+2/step_h^2);
        if (phi(1,k,this_time)-phi(1,k,last_time)) > maxerror_V
            maxerror_V = phi(1,k,this_time)-phi(1,k,last_time);
        end
        phi(length_a,k,this_time) = phi0_V;
    end
    for i = 2:length_a-1
        for k = 2:length_h-1
            phi(i,k,this_time) = ((phi(i+1,k,last_time)+phi(i-1,k,this_time))/step_a^2+(phi(i+1,k,last_time)-phi(i-1,k,this_time))/(2*(i-1)*step_a^2)+(phi(i,k+1,last_time)+phi(i,k-1,this_time))/step_h^2)/(2/step_a^2+2/step_h^2);
            if (phi(i,k,this_time)-phi(i,k,last_time)) > maxerror_V
                maxerror_V = phi(1,k,this_time)-phi(1,k,last_time);
            end
        end
    end
    clc;
    disp(maxerror_V);
end

phi_final = zeros(2*length_a-1,length_h);
for i = 1:2*length_a-1
    phi_final(i,:) = phi(abs(i-length_a)+1,:,this_time);
end

R = -a_m:step_a:a_m;
Z = 0:step_h:h_m;
figure
set(gcf,'Position',[300 100 600 900]);
contourf(R,Z,phi_final(:,:)',200,'Linestyle','none');
title('金属圆筒内部电场分布','FontSize',15);
xlabel('\it r /\rm m');
ylabel('\it z /\rm m');

hold on

[Er,Ez] = gradient(-phi_final(:,:)',step_a,step_h);

quiver(R(2:2*length_a-2),Z(1:length_h),Er(1:length_h,2:2*length_a-2),Ez(1:length_h,2:2*length_a-2),'Color','k','LineWidth',1.2);


r_inf = 5;
z_inf = 15;
h_gap = 0.02;
metal_a = a_m/step_a+1;
metal_h = h_m/step_h+1;
gap_h = (h_m-h_gap)/step_h+1;
length_a = r_inf/step_a+1;
length_h = z_inf/step_h+1;

phi = zeros(length_a, length_h, 2);

for i=1:metal_a
    phi(i,metal_h,1) = phih_V;
end

times = 0;
maxerror_V = 1000;
tolerance_V = 0.0001;
this_time = 1;

while maxerror_V > tolerance_V
    times = times + 1;
    last_time = this_time;
    this_time = mod(times,2)+1;
    maxerror_V = 0;
    for i = 1:length_a
        phi(i,1,this_time) = phi0_V;
        phi(i,length_h,this_time) = phi0_V;
    end
    for i = 1:metal_a
        phi(i,metal_h,this_time) = phih_V;
    end
    for k = 2:length_h-1
        if k~=metal_h
            phi(1,k,this_time) = (2*phi(2,k,last_time)/step_a^2+(phi(1,k+1,last_time)+phi(1,k-1,this_time))/step_h^2)/(2/step_a^2+2/step_h^2);
            if (phi(1,k,this_time)-phi(1,k,last_time)) > maxerror_V
                maxerror_V = phi(1,k,this_time)-phi(1,k,last_time);
            end
        end
        phi(length_a,k,this_time) = phi0_V;
    end
    for k = 2:gap_h
        phi(metal_a,k,this_time) = phi0_V;
    end
    for i = 2:length_a-1
        for k = 2:length_h-1
            if ((k~=metal_h)&&(i~=metal_a))||((i==metal_a)&&(k<=gap_h))
                phi(i,k,this_time) = ((phi(i+1,k,last_time)+phi(i-1,k,this_time))/step_a^2+(phi(i+1,k,last_time)-phi(i-1,k,this_time))/(2*(i-1)*step_a^2)+(phi(i,k+1,last_time)+phi(i,k-1,this_time))/step_h^2)/(2/step_a^2+2/step_h^2);
                if (phi(i,k,this_time)-phi(i,k,last_time)) > maxerror_V
                    maxerror_V = phi(1,k,this_time)-phi(1,k,last_time);
                end
            end
        end
    end
    clc;
    disp(maxerror_V);
end

draw_a = metal_a*2;
draw_h = metal_h*2;

phi_final = zeros(2*draw_a-1,draw_h);
for i = 1:2*draw_a-1
    phi_final(i,:) = phi(abs(i-draw_a)+1,1:draw_h,this_time);
end

R = -(draw_a-1)*step_a:step_a:(draw_a-1)*step_a;
Z = 0:step_h:(draw_h-1)*step_h;
figure
set(gcf,'Position',[300 100 600 900]);
contourf(R,Z,phi_final(:,:)',200,'Linestyle','none');
title('金属圆筒内外电场分布','FontSize',15);
xlabel('\it r /\rm m');
ylabel('\it z /\rm m');

hold on

[Er,Ez] = gradient(-phi_final(:,:)',step_a,step_h);

quiver(R(2:2*draw_a-2),Z(1:draw_h),Er(1:draw_h,2:2*draw_a-2),Ez(1:draw_h,2:2*draw_a-2),'LineWidth',1);
plot([a_m a_m],[0 h_m-h_gap],'k-','LineWidth',2);
plot([-a_m -a_m],[0 h_m-h_gap],'k-','LineWidth',2);
plot([-a_m a_m],[h_m h_m],'k-','LineWidth',2);



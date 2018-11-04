clearvars;

%设置金属圆筒尺寸
a_m = 0.1;
h_m = 0.3;
%设置边界电势值
phih_V = 100;
phi0_V = 0;

%设置划分场域的步长
step_a = 0.005;
step_h = 0.005;
%根据步长得到需要的数组长度
length_a = a_m/step_a+1;
length_h = h_m/step_h+1;

%phi[i,k,t] -- 存储电势值
%      i    -- r方向坐标
%      k    -- z方向坐标
%      t    -- 存储当前和上一次迭代值
%初始化电势均为零
phi = zeros(length_a, length_h, 2);
%设定顶盖电势
for i=1:length_a-1
    phi(i,length_h,1) = phih_V;
end

%记录迭代次数
times = 0;
%记录两次迭代间的最大误差
maxerror_V = 1000;
%设定容许值
tolerance_V = 0.0001;
%当前电势结果存储的位置指针
this_time = 1;

%迭代开始
while maxerror_V > tolerance_V
    times = times + 1;
    %确定当前和上一次迭代数据存储的位置指针
    last_time = this_time;
    this_time = mod(times,2)+1;
    maxerror_V = 0;
    %设定顶盖和底部的边界条件
    for i = 1:length_a
        phi(i,1,this_time) = phi0_V;
        phi(i,length_h,this_time) = phih_V;
    end
    phi(length_a,length_h,this_time) = phi0_V;
    %设定中心轴线和侧壁的边界条件
    for k = 2:length_h-1
        phi(1,k,this_time) = (2*phi(2,k,last_time)/step_a^2+(phi(1,k+1,last_time)+phi(1,k-1,this_time))/step_h^2)/(2/step_a^2+2/step_h^2);
        if (phi(1,k,this_time)-phi(1,k,last_time)) > maxerror_V
            maxerror_V = phi(1,k,this_time)-phi(1,k,last_time);
        end
        phi(length_a,k,this_time) = phi0_V;
    end
    %对场域内各节点进行迭代计算
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

%取出迭代结果，展开到整个对称截面，便于观察
phi_final = zeros(2*length_a-1,length_h);
for i = 1:2*length_a-1
    phi_final(i,:) = phi(abs(i-length_a)+1,:,this_time);
end
%给出实际尺寸坐标
R = -a_m:step_a:a_m;
Z = 0:step_h:h_m;
figure
set(gcf,'Position',[300 100 600 900]);
%绘制等势线，由于i对应横向半径，k对应纵向高度，为了美观，将矩阵转置
contour(R,Z,phi_final(:,:)',20);
title('金属圆筒内部电场分布','FontSize',15);
xlabel('\it r /\rm m');
ylabel('\it z /\rm m');
hold on
%求得电场强度
[Er,Ez] = gradient(-phi_final(:,:)',step_a,step_h);
%绘制电场强度矢量
quiver(R(2:2*length_a-2),Z(1:length_h),Er(1:length_h,2:2*length_a-2),Ez(1:length_h,2:2*length_a-2),'Color','k','LineWidth',1.2);

%另一种风格的等势线
figure
set(gcf,'Position',[300 100 600 900]);
contourf(R,Z,phi_final(:,:)',200,'LineStyle','None');
title('金属圆筒内部电场分布','FontSize',15);
xlabel('\it r /\rm m');
ylabel('\it z /\rm m');
hold on
[Er,Ez] = gradient(-phi_final(:,:)',step_a,step_h);
quiver(R(2:2*length_a-2),Z(1:length_h),Er(1:length_h,2:2*length_a-2),Ez(1:length_h,2:2*length_a-2),'Color','k','LineWidth',1.2);

%计算开放场域电场分布
%设定认为的无穷远位置
r_inf = 5;
z_inf = 15;
%设定缝隙尺寸
a_gap = 0.01;
%设定划分场域步长
step_a = a_gap/5;
step_h = a_gap/5;
%设定金属圆筒边界对应的坐标
metal_a = a_m/step_a+1;
metal_h = h_m/step_h+1;
gap_a = round((a_m-a_gap)/step_a+1);
%设定需要的数组长度
length_a = r_inf/step_a+1;
length_h = z_inf/step_h+1;

%初始化电势值
phi = zeros(length_a, length_h, 2);
for i=1:gap_a
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
    %设定无穷高处和地面边界条件
    for i = 1:length_a
        phi(i,1,this_time) = phi0_V;
        phi(i,length_h,this_time) = phi0_V;
    end
    %设定顶盖边界条件
    for i = 1:gap_a
        phi(i,metal_h,this_time) = phih_V;
    end
    %设定中心轴线和侧面无穷远处边界条件
    for k = 2:length_h-1
        %剔除顶盖上的节点
        if k~=metal_h
            phi(1,k,this_time) = (2*phi(2,k,last_time)/step_a^2+(phi(1,k+1,last_time)+phi(1,k-1,this_time))/step_h^2)/(2/step_a^2+2/step_h^2);
            if (phi(1,k,this_time)-phi(1,k,last_time)) > maxerror_V
                maxerror_V = phi(1,k,this_time)-phi(1,k,last_time);
            end
        end
        phi(length_a,k,this_time) = phi0_V;
    end
    %设定圆筒侧壁边界条件
    for k = 2:metal_h
        phi(metal_a,k,this_time) = phi0_V;
    end
    %进行一轮迭代计算
    for i = 2:length_a-1
        for k = 2:length_h-1
            %剔除金属部分节点
            if ((k~=metal_h)&&(i~=metal_a))||((k==metal_h)&&(i>gap_a)&&(i<metal_a))
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

%设定绘制场域大小
draw_a = metal_a*2;
draw_h = round(metal_h*1.6667);

%去除电势结果
phi_final = zeros(2*draw_a-1,draw_h);
for i = 1:2*draw_a-1
    phi_final(i,:) = phi(abs(i-draw_a)+1,1:draw_h,this_time);
end

R = -(draw_a-1)*step_a:step_a:(draw_a-1)*step_a;
Z = 0:step_h:(draw_h-1)*step_h;
%绘制等势线
figure
set(gcf,'Position',[300 100 600 900]);
contour(R,Z,phi_final(:,:)',20);
title('金属圆筒内外电场分布','FontSize',15);
xlabel('\it r /\rm m');
ylabel('\it z /\rm m');
hold on
%绘制电场强度
[Er,Ez] = gradient(-phi_final(:,:)',step_a,step_h);
quiver(R(2:2*draw_a-2),Z(1:draw_h),Er(1:draw_h,2:2*draw_a-2),Ez(1:draw_h,2:2*draw_a-2),'Color','k');
%画出金属圆筒
plot([a_m a_m],[0 h_m],'k-','LineWidth',3);
plot([-a_m -a_m],[0 h_m],'k-','LineWidth',3);
plot([-(a_m-a_gap) a_m-a_gap],[h_m h_m],'k-','LineWidth',3);

%绘制另一种风格等势线
figure
set(gcf,'Position',[300 100 600 900]);
contourf(R,Z,phi_final(:,:)',200,'LineStyle','None');
title('金属圆筒内外电场分布','FontSize',15);
xlabel('\it r /\rm m');
ylabel('\it z /\rm m');
hold on
[Er,Ez] = gradient(-phi_final(:,:)',step_a,step_h);
quiver(R(2:2*draw_a-2),Z(1:draw_h),Er(1:draw_h,2:2*draw_a-2),Ez(1:draw_h,2:2*draw_a-2),'Color','k');
plot([a_m a_m],[0 h_m],'k-','LineWidth',3);
plot([-a_m -a_m],[0 h_m],'k-','LineWidth',3);
plot([-(a_m-a_gap) a_m-a_gap],[h_m h_m],'k-','LineWidth',3);


%计算电容
%设定介电常数
e = 8.854187817e-12;
%初始化电荷面密度数组
q = zeros(gap_a,1);
%根据顶盖表面Z轴方向电场强度，计算电荷面密度
for i=1:gap_a
    q(i) = e*(abs(Ez(metal_h-1,draw_a+i-1))+abs(Ez(metal_h+1,draw_a+i-1)));
end
%积分得到总电荷量
Q = q(1)*pi*step_a^2 + q(gap_a)*pi*((a_m-a_gap)^2-(a_m-a_gap-step_a)^2);
for i=2:gap_a-1
    Q = Q + q(i)*pi*(((i-0.5)*step_a)^2-((i-1.5)*step_a)^2);
end
%由定义式计算得到电容值
C = Q/(phih_V-phi0_V);
disp(['电容大小为',num2str(C*1E12),'pF']);

%绘制电荷分布
figure
set(gcf,'Position',[300 100 800 800]);
%生成柱坐标计算网格
theta_n = 200;
theta = linspace(0,2*pi,theta_n);
r = 0:step_a:a_m-a_gap;
[Theta,R] = meshgrid(theta,r);
%填入电荷面密度
z = zeros(gap_a,theta_n);
for i=1:gap_a
    z(i,:) = q(i).*ones(1,theta_n);
end
%转化到直角坐标
[X,Y,Z] = pol2cart(Theta,R,z);
%绘制
contourf(X,Y,Z,200,'LineStyle','None');
axis([-a_m*1.1 a_m*1.1 -a_m*1.1 a_m*1.1]);
title('金属圆筒顶盖电荷分布','FontSize',15);



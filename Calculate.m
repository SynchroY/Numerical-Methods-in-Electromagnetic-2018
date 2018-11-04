clearvars;

%���ý���ԲͲ�ߴ�
a_m = 0.1;
h_m = 0.3;
%���ñ߽����ֵ
phih_V = 100;
phi0_V = 0;

%���û��ֳ���Ĳ���
step_a = 0.005;
step_h = 0.005;
%���ݲ����õ���Ҫ�����鳤��
length_a = a_m/step_a+1;
length_h = h_m/step_h+1;

%phi[i,k,t] -- �洢����ֵ
%      i    -- r��������
%      k    -- z��������
%      t    -- �洢��ǰ����һ�ε���ֵ
%��ʼ�����ƾ�Ϊ��
phi = zeros(length_a, length_h, 2);
%�趨���ǵ���
for i=1:length_a-1
    phi(i,length_h,1) = phih_V;
end

%��¼��������
times = 0;
%��¼���ε������������
maxerror_V = 1000;
%�趨����ֵ
tolerance_V = 0.0001;
%��ǰ���ƽ���洢��λ��ָ��
this_time = 1;

%������ʼ
while maxerror_V > tolerance_V
    times = times + 1;
    %ȷ����ǰ����һ�ε������ݴ洢��λ��ָ��
    last_time = this_time;
    this_time = mod(times,2)+1;
    maxerror_V = 0;
    %�趨���Ǻ͵ײ��ı߽�����
    for i = 1:length_a
        phi(i,1,this_time) = phi0_V;
        phi(i,length_h,this_time) = phih_V;
    end
    phi(length_a,length_h,this_time) = phi0_V;
    %�趨�������ߺͲ�ڵı߽�����
    for k = 2:length_h-1
        phi(1,k,this_time) = (2*phi(2,k,last_time)/step_a^2+(phi(1,k+1,last_time)+phi(1,k-1,this_time))/step_h^2)/(2/step_a^2+2/step_h^2);
        if (phi(1,k,this_time)-phi(1,k,last_time)) > maxerror_V
            maxerror_V = phi(1,k,this_time)-phi(1,k,last_time);
        end
        phi(length_a,k,this_time) = phi0_V;
    end
    %�Գ����ڸ��ڵ���е�������
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

%ȡ�����������չ���������Գƽ��棬���ڹ۲�
phi_final = zeros(2*length_a-1,length_h);
for i = 1:2*length_a-1
    phi_final(i,:) = phi(abs(i-length_a)+1,:,this_time);
end
%����ʵ�ʳߴ�����
R = -a_m:step_a:a_m;
Z = 0:step_h:h_m;
figure
set(gcf,'Position',[300 100 600 900]);
%���Ƶ����ߣ�����i��Ӧ����뾶��k��Ӧ����߶ȣ�Ϊ�����ۣ�������ת��
contour(R,Z,phi_final(:,:)',20);
title('����ԲͲ�ڲ��糡�ֲ�','FontSize',15);
xlabel('\it r /\rm m');
ylabel('\it z /\rm m');
hold on
%��õ糡ǿ��
[Er,Ez] = gradient(-phi_final(:,:)',step_a,step_h);
%���Ƶ糡ǿ��ʸ��
quiver(R(2:2*length_a-2),Z(1:length_h),Er(1:length_h,2:2*length_a-2),Ez(1:length_h,2:2*length_a-2),'Color','k','LineWidth',1.2);

%��һ�ַ��ĵ�����
figure
set(gcf,'Position',[300 100 600 900]);
contourf(R,Z,phi_final(:,:)',200,'LineStyle','None');
title('����ԲͲ�ڲ��糡�ֲ�','FontSize',15);
xlabel('\it r /\rm m');
ylabel('\it z /\rm m');
hold on
[Er,Ez] = gradient(-phi_final(:,:)',step_a,step_h);
quiver(R(2:2*length_a-2),Z(1:length_h),Er(1:length_h,2:2*length_a-2),Ez(1:length_h,2:2*length_a-2),'Color','k','LineWidth',1.2);

%���㿪�ų���糡�ֲ�
%�趨��Ϊ������Զλ��
r_inf = 5;
z_inf = 15;
%�趨��϶�ߴ�
a_gap = 0.01;
%�趨���ֳ��򲽳�
step_a = a_gap/5;
step_h = a_gap/5;
%�趨����ԲͲ�߽��Ӧ������
metal_a = a_m/step_a+1;
metal_h = h_m/step_h+1;
gap_a = round((a_m-a_gap)/step_a+1);
%�趨��Ҫ�����鳤��
length_a = r_inf/step_a+1;
length_h = z_inf/step_h+1;

%��ʼ������ֵ
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
    %�趨����ߴ��͵���߽�����
    for i = 1:length_a
        phi(i,1,this_time) = phi0_V;
        phi(i,length_h,this_time) = phi0_V;
    end
    %�趨���Ǳ߽�����
    for i = 1:gap_a
        phi(i,metal_h,this_time) = phih_V;
    end
    %�趨�������ߺͲ�������Զ���߽�����
    for k = 2:length_h-1
        %�޳������ϵĽڵ�
        if k~=metal_h
            phi(1,k,this_time) = (2*phi(2,k,last_time)/step_a^2+(phi(1,k+1,last_time)+phi(1,k-1,this_time))/step_h^2)/(2/step_a^2+2/step_h^2);
            if (phi(1,k,this_time)-phi(1,k,last_time)) > maxerror_V
                maxerror_V = phi(1,k,this_time)-phi(1,k,last_time);
            end
        end
        phi(length_a,k,this_time) = phi0_V;
    end
    %�趨ԲͲ��ڱ߽�����
    for k = 2:metal_h
        phi(metal_a,k,this_time) = phi0_V;
    end
    %����һ�ֵ�������
    for i = 2:length_a-1
        for k = 2:length_h-1
            %�޳��������ֽڵ�
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

%�趨���Ƴ����С
draw_a = metal_a*2;
draw_h = round(metal_h*1.6667);

%ȥ�����ƽ��
phi_final = zeros(2*draw_a-1,draw_h);
for i = 1:2*draw_a-1
    phi_final(i,:) = phi(abs(i-draw_a)+1,1:draw_h,this_time);
end

R = -(draw_a-1)*step_a:step_a:(draw_a-1)*step_a;
Z = 0:step_h:(draw_h-1)*step_h;
%���Ƶ�����
figure
set(gcf,'Position',[300 100 600 900]);
contour(R,Z,phi_final(:,:)',20);
title('����ԲͲ����糡�ֲ�','FontSize',15);
xlabel('\it r /\rm m');
ylabel('\it z /\rm m');
hold on
%���Ƶ糡ǿ��
[Er,Ez] = gradient(-phi_final(:,:)',step_a,step_h);
quiver(R(2:2*draw_a-2),Z(1:draw_h),Er(1:draw_h,2:2*draw_a-2),Ez(1:draw_h,2:2*draw_a-2),'Color','k');
%��������ԲͲ
plot([a_m a_m],[0 h_m],'k-','LineWidth',3);
plot([-a_m -a_m],[0 h_m],'k-','LineWidth',3);
plot([-(a_m-a_gap) a_m-a_gap],[h_m h_m],'k-','LineWidth',3);

%������һ�ַ�������
figure
set(gcf,'Position',[300 100 600 900]);
contourf(R,Z,phi_final(:,:)',200,'LineStyle','None');
title('����ԲͲ����糡�ֲ�','FontSize',15);
xlabel('\it r /\rm m');
ylabel('\it z /\rm m');
hold on
[Er,Ez] = gradient(-phi_final(:,:)',step_a,step_h);
quiver(R(2:2*draw_a-2),Z(1:draw_h),Er(1:draw_h,2:2*draw_a-2),Ez(1:draw_h,2:2*draw_a-2),'Color','k');
plot([a_m a_m],[0 h_m],'k-','LineWidth',3);
plot([-a_m -a_m],[0 h_m],'k-','LineWidth',3);
plot([-(a_m-a_gap) a_m-a_gap],[h_m h_m],'k-','LineWidth',3);


%�������
%�趨��糣��
e = 8.854187817e-12;
%��ʼ��������ܶ�����
q = zeros(gap_a,1);
%���ݶ��Ǳ���Z�᷽��糡ǿ�ȣ����������ܶ�
for i=1:gap_a
    q(i) = e*(abs(Ez(metal_h-1,draw_a+i-1))+abs(Ez(metal_h+1,draw_a+i-1)));
end
%���ֵõ��ܵ����
Q = q(1)*pi*step_a^2 + q(gap_a)*pi*((a_m-a_gap)^2-(a_m-a_gap-step_a)^2);
for i=2:gap_a-1
    Q = Q + q(i)*pi*(((i-0.5)*step_a)^2-((i-1.5)*step_a)^2);
end
%�ɶ���ʽ����õ�����ֵ
C = Q/(phih_V-phi0_V);
disp(['���ݴ�СΪ',num2str(C*1E12),'pF']);

%���Ƶ�ɷֲ�
figure
set(gcf,'Position',[300 100 800 800]);
%�����������������
theta_n = 200;
theta = linspace(0,2*pi,theta_n);
r = 0:step_a:a_m-a_gap;
[Theta,R] = meshgrid(theta,r);
%���������ܶ�
z = zeros(gap_a,theta_n);
for i=1:gap_a
    z(i,:) = q(i).*ones(1,theta_n);
end
%ת����ֱ������
[X,Y,Z] = pol2cart(Theta,R,z);
%����
contourf(X,Y,Z,200,'LineStyle','None');
axis([-a_m*1.1 a_m*1.1 -a_m*1.1 a_m*1.1]);
title('����ԲͲ���ǵ�ɷֲ�','FontSize',15);



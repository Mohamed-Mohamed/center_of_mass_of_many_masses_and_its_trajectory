%% Coded by
% Mohamed Mohamed El-Sayed Atyya
% mohamed.atyya94@eng-st.cu.edu.eg


% this program is used for geting center of mass of two masses and its trajoctry and velocity at each point using RK4
% numerical integration method under Newton's gravitaional force
close all; clear all; clc;
%% constants
G=6.67428E-11; % gravitational constant
%% Intial condition
M=[5.974e24,73.48e21];       % [M1,M2]
R1=[0;0;0];                 % position of M(1)
R2=[384.4e6;0;0];                 % position of M(2)
I=[0,0,0];              % location of initial axis
V1=[0;0;0];                 % velocity of M(1)
V2=[0;1;0]*1.02306e3;                 % velocity of M(2)
%% RK4 parameter
tf=365*24*360*5;   % final time of soution
dt=3600;            % time step
X0=[R1;R2;V1;V2];
B=[0;0;0;0;0;0;0;0;0;0;0;0];
sol(1:12,1)=X0;
order=12;
%% solution by RK4
for n=1:length(0:dt:tf)
    b=G*M(2)/(norm(sol(1:3,n)-sol(4:6,n)))^3;
    c=-G*M(1)/(norm(sol(1:3,n)-sol(4:6,n)))^3;
    A=[0,0,0,0,0,0,1,0,0,0,0,0; ...
        0,0,0,0,0,0,0,1,0,0,0,0; ...
        0,0,0,0,0,0,0,0,1,0,0,0; ...
        0,0,0,0,0,0,0,0,0,1,0,0; ...
        0,0,0,0,0,0,0,0,0,0,1,0; ...
        0,0,0,0,0,0,0,0,0,0,0,1;...
        -b,0,0,b,0,0,0,0,0,0,0,0; ...
        0,-b,0,0,b,0,0,0,0,0,0,0; ...
        0,0,-b,0,0,b,0,0,0,0,0,0; ...
        -c,0,0,c,0,0,0,0,0,0,0,0; ...
        0,-c,0,0,c,0,0,0,0,0,0,0; ...
        0,0,-c,0,0,c,0,0,0,0,0,0 ];
    [ XX ] = RK4( A,B,sol(1:12,n),dt,n*dt,(n+1)*dt,order );
    sol(1:12,n+1)=XX(1:12,2);
end
R1_x=sol(1,:);
R1_y=sol(2,:);
R1_z=sol(3,:);
R2_x=sol(4,:);
R2_y=sol(5,:);
R2_z=sol(6,:);
V1_x=sol(7,:);
V1_y=sol(8,:);
V1_z=sol(9,:);
V2_x=sol(10,:);
V2_y=sol(11,:);
V2_z=sol(12,:);
%% center of masses parameters
r=[R2_x;R2_y;R2_z]-[R1_x;R1_y;R1_z];                                            % the distance betweem M1 & M2
Rc=(M(1)*[R1_x;R1_y;R1_z]+M(2)*[R2_x;R2_y;R2_z])/sum(M);                        % location of center of masses
Vc=(M(1)*[V1_x;V1_y;V1_z]+M(2)*[V2_x;V2_y;V2_z])/sum(M);                        % Vc = constant
Ac=(M(1)*G*M(2)*r.^3.*r-M(1)*G*M(2)*r.^3.*r)/sum(M);    % for check Ac = 0
for H=1:length(Rc(1,:))
    acc1(1:3,H)=(-1)^(1)*G*M(1)/(norm([R2_x(H);R2_y(H);R2_z(H)]'-[R1_x(H);R1_y(H);R1_z(H)]'))^3*([R2_x(H);R2_y(H);R2_z(H)]-[R1_x(H);R1_y(H);R1_z(H)]); % acceleration of M1
    acc2(1:3,H)=(-1)^(2)*G*M(2)/(norm([R2_x(H);R2_y(H);R2_z(H)]'-[R1_x(H);R1_y(H);R1_z(H)]'))^3*([R2_x(H);R2_y(H);R2_z(H)]-[R1_x(H);R1_y(H);R1_z(H)]); % acceleration of M2
end
%% velocity and acceleration magnituides
for h=1:length(Rc(1,:))
        MagV1(1,h)=norm([V1_x;V1_y;V1_z]);  % velocity magnituide of M1
        MagV2(1,h)=norm([V2_x;V2_y;V2_z]);  % velocity magnituide of M2
        MagA1(1,h)=norm(acc1(1:3,h));  % acceleration magnituide of M1
        MagA2(1,h)=norm(acc2(1:3,h));  % acceleration magnituide of M2
end
%% plotting
%--------------------------------------------------------------------------------------------------------------------------------------------------------
% axes at postion I
figure(1);
view(3);
set(gcf,'Color','w');
hold all;
% M1 trajectory from I
plot3(I(1)*ones(1,length(R1_x))+R1_x,I(2)*ones(1,length(R1_y))+R1_y,I(3)*ones(1,length(R1_z))+R1_z,'cyan','LineWidth',2);
% M1 start
plot3(I(1)+R1_x(1),I(2)+R1_y(1),I(3)+R1_z(1),'o','color',[.5,0.2,0.9],'LineWidth',5);
% M2 trajectory from I
plot3(I(1)*ones(1,length(R2_x))+R2_x,I(2)*ones(1,length(R2_y))+R2_y,I(3)*ones(1,length(R2_z))+R2_z,'g','LineWidth',2);
% M2 start
plot3(I(1)+R2_x(1),I(2)+R2_y(1),I(3)+R2_z(1),'o','color',[.9,0.2,0.5],'LineWidth',5);
% Rc trajectory from I
plot3(I(1)*ones(1,length(Rc(1,:)))+Rc(1,:),I(2)*ones(1,length(Rc(2,:)))+Rc(2,:),I(3)*ones(1,length(Rc(3,:)))+Rc(3,:),'red','LineWidth',2);
grid on;
xlabel('X','Fontsize',18);
ylabel('Y','Fontsize',18);
zlabel('Z','Fontsize',18);
title('Solution at Position I','Fontsize',18);
xlim auto;
ylim auto;
zlim auto;
legend('M1 trajectory fom I','M1 start','M2 trajectory from I','M2 start','Rc trajectory from I');
%--------------------------------------------------------------------------------------------------------------------------------------------------------
% axes at M1
figure(2);
view(3);
set(gcf,'Color','w');
hold all;
% M1 position
plot3(0,0,0,'o','color','cyan','LineWidth',5);
% M2 position relative to M1
plot3(R2_x-R1_x,R2_y-R1_y,R2_z-R1_z,'color','g','LineWidth',2);
% M2 start
plot3(R2_x(1)-R1_x(1),R2_y(1)-R1_y(1),R2_z(1)-R1_z(1),'o','color',[.9,0.2,0.5],'LineWidth',5);
% Rc position relative to M1
plot3(Rc(1,:)-R1_x,Rc(2,:)-R1_y,Rc(3,:)-R1_z,'red','LineWidth',2);
grid on;
xlabel('X','Fontsize',18);
ylabel('Y','Fontsize',18);
zlabel('Z','Fontsize',18);
title('Solution at M1 Center of Mass','Fontsize',18);
xlim auto;
ylim auto;
zlim auto;
legend('M1 position','M2 position relative to M1','M2 start','Rc position relative to M1');
%--------------------------------------------------------------------------------------------------------------------------------------------------------
% axes at M2
figure(3);
view(3);
set(gcf,'Color','w');
hold all;
% M1 position relative to M2 curve
plot3(R1_x-R2_x,R1_y-R2_y,R1_z-R2_z,'color','cyan','LineWidth',2);
% M1 start
plot3(R1_x(1)-R2_x(1),R1_y(1)-R2_y(1),R1_z(1)-R2_z(1),'o','color',[.5,0.2,0.9],'LineWidth',5);
%M2 position
plot3(0,0,0,'o','color','g','LineWidth',5);
%Rc position relative to M2 curve
plot3(Rc(1,:)-R2_x,Rc(2,:)-R2_y,Rc(3,:)-R2_z,'red','LineWidth',2);
grid on;
xlabel('X','Fontsize',18);
ylabel('Y','Fontsize',18);
zlabel('Z','Fontsize',18);
title('Solution at M2 Center of Mass','Fontsize',18);
xlim auto;
ylim auto;
zlim auto;
legend('M1 position relative to M1','M1 start','M2 position','Rc position relative to M2');
%--------------------------------------------------------------------------------------------------------------------------------------------------------
% axes at Rc
figure(4);
view(3);
set(gcf,'Color','w');
hold all;
% change in M1 position relative to Rc
plot3(R1_x-Rc(1,:),R1_y-Rc(2,:),R1_z-Rc(3,:),'color','cyan','LineWidth',2);
% M1 start
plot3(R1_x(1)-Rc(1,1),R1_y(1)-Rc(2,1),R1_z(1)-Rc(3,1),'o','color',[.5,0.2,0.9],'LineWidth',5);
% change in M2 position relative to Rc
plot3(R2_x-Rc(1,:),R2_y-Rc(2,:),R2_z-Rc(3,:),'color','g','LineWidth',2);
% M2 start
plot3(R2_x(1)-Rc(1,1),R2_y(1)-Rc(2,1),R2_z(1)-Rc(3,1),'o','color',[.9,0.2,0.5],'LineWidth',5);
% Rc position
plot3(0,0,0,'o','color','red','LineWidth',5);
grid on;
xlabel('X','Fontsize',18);
ylabel('Y','Fontsize',18);
zlabel('Z','Fontsize',18);
title('Solution at Rc','Fontsize',18);
xlim auto;
ylim auto;
zlim auto;
legend('M1 position relative to Rc','M1 start','M2 position relative to Rc','M2 start','Rc position');
%--------------------------------------------------------------------------------------------------------------------------------------------------------
% All
% axes at postion I
figure(5);
subplot(2,2,1);
view(3);
set(gcf,'Color','w');
hold all;
% M1 trajectory from I
plot3(I(1)*ones(1,length(R1_x))+R1_x,I(2)*ones(1,length(R1_y))+R1_y,I(3)*ones(1,length(R1_z))+R1_z,'cyan','LineWidth',2);
% M2 trajectory from I
plot3(I(1)*ones(1,length(R2_x))+R2_x,I(2)*ones(1,length(R2_y))+R2_y,I(3)*ones(1,length(R2_z))+R2_z,'g','LineWidth',2);
% Rc trajectory from I
plot3(I(1)*ones(1,length(Rc(1,:)))+Rc(1,:),I(2)*ones(1,length(Rc(2,:)))+Rc(2,:),I(3)*ones(1,length(Rc(3,:)))+Rc(3,:),'red','LineWidth',2);
grid on;
xlabel('X','Fontsize',18);
ylabel('Y','Fontsize',18);
zlabel('Z','Fontsize',18);
title('Solution at Position I','Fontsize',18);
xlim auto;
ylim auto;
zlim auto;
legend('M1 trajectory fom I','M2 trajectory from I','Rc trajectory from I','Location','NorthWest');
% axes at M1
subplot(2,2,2);
view(3);
set(gcf,'Color','w');
hold all;
% M1 position
plot3(0,0,0,'o','color','cyan','LineWidth',5);
% M2 position relative to M1
plot3(R2_x-R1_x,R2_y-R1_y,R2_z-R1_z,'color','g','LineWidth',2);
% Rc position relative to M1
plot3(Rc(1,:)-R1_x,Rc(2,:)-R1_y,Rc(3,:)-R1_z,'red','LineWidth',2);
grid on;
xlabel('X','Fontsize',18);
ylabel('Y','Fontsize',18);
zlabel('Z','Fontsize',18);
title('Solution at M1 Center of Mass','Fontsize',18);
xlim auto;
ylim auto;
zlim auto;
legend('M1 position','M2 position relative to M1','Rc position relative to M1','Location','NorthEast');
% axes at M2
subplot(2,2,3);
view(3);
set(gcf,'Color','w');
hold all;
% M1 position relative to M2 curve
plot3(R1_x-R2_x,R1_y-R2_y,R1_z-R2_z,'color','cyan','LineWidth',2);
%M2 position
plot3(0,0,0,'o','color','g','LineWidth',5);
%Rc position relative to M2 curve
plot3(Rc(1,:)-R2_x,Rc(2,:)-R2_y,Rc(3,:)-R2_z,'red','LineWidth',2);
grid on;
xlabel('X','Fontsize',18);
ylabel('Y','Fontsize',18);
zlabel('Z','Fontsize',18);
title('Solution at M2 Center of Mass','Fontsize',18);
xlim auto;
ylim auto;
zlim auto;
legend('M1 position relative to M1','M2 position','Rc position relative to M2','Location','NorthWest');
% axes at Rc
subplot(2,2,4);
view(3);
set(gcf,'Color','w');
hold all;
% change in M1 position relative to Rc
plot3(R1_x-Rc(1,:),R1_y-Rc(2,:),R1_z-Rc(3,:),'color','cyan','LineWidth',2);
% change in M2 position relative to Rc
plot3(R2_x-Rc(1,:),R2_y-Rc(2,:),R2_z-Rc(3,:),'color','g','LineWidth',2);
% Rc position
plot3(0,0,0,'o','color','red','LineWidth',5);
grid on;
xlabel('X','Fontsize',18);
ylabel('Y','Fontsize',18);
zlabel('Z','Fontsize',18);
title('Solution at Rc','Fontsize',18);
xlim auto;
ylim auto;
zlim auto;
legend('M1 position relative to Rc','M2 position relative to Rc','Rc position','Location','NorthEast');
%--------------------------------------------------------------------------------------------------------------------------------------------------------
figure(6);
subplot(2,2,1);
view(3);
set(gcf,'Color','w');
hold all;
% M1 position
plot3(0,0,0,'o','color','cyan','LineWidth',5);
% M2 position relative to M1
plot3(R2_x-R1_x,R2_y-R1_y,R2_z-R1_z,'color','g','LineWidth',2);
% Rc position relative to M1
plot3(Rc(1,:)-R1_x,Rc(2,:)-R1_y,Rc(3,:)-R1_z,'red','LineWidth',2);
%Velocity distribution of M2
stem3(R2_x-R1_x,R2_y-R1_y,R2_z-R1_z+MagV2,':^m','color','blue');
grid on;
xlabel('X','Fontsize',18);
ylabel('Y','Fontsize',18);
zlabel('Z','Fontsize',18);
xlim auto;
ylim auto;
zlim auto;
title('Solution at M1 with velocity distribution of M2','Fontsize',12);
legend('M1 position','M2 position relative to M1','Rc position relative to M1','Velocity distribution of M2','Location','NorthEast');
%--------------------------------------------------------------------------------------------------------------------------------------------------------
subplot(2,2,2);
view(3);
set(gcf,'Color','w');
hold all;
% M1 position relative to M2 curve
plot3(R1_x-R2_x,R1_y-R2_y,R1_z-R2_z,'color','cyan','LineWidth',2);
%M2 position
plot3(0,0,0,'o','color','g','LineWidth',5);
%Rc position relative to M2 curve
plot3(Rc(1,:)-R2_x,Rc(2,:)-R2_y,Rc(3,:)-R2_z,'red','LineWidth',2);
% Velocity distribution of M1
stem3(R1_x-R2_x,R1_y-R2_y,R1_z-R2_z+MagV1,':^m','color','blue');
grid on;
xlabel('X','Fontsize',18);
ylabel('Y','Fontsize',18);
zlabel('Z','Fontsize',18);
xlim auto;
ylim auto;
zlim auto;
title('Solution at M2 with vlocity distribution of M1','Fontsize',12);
legend('M1 position relative to M1','M2 position','Rc position relative to M2','Velocity distribution of M1','Location','NorthEast');
%--------------------------------------------------------------------------------------------------------------------------------------------------------
subplot(2,2,3);
view(3);
set(gcf,'Color','w');
hold all;
% M1 position
plot3(0,0,0,'o','color','cyan','LineWidth',5);
% M2 position relative to M1
plot3(R2_x-R1_x,R2_y-R1_y,R2_z-R1_z,'color','g','LineWidth',2);
% Rc position relative to M1
plot3(Rc(1,:)-R1_x,Rc(2,:)-R1_y,Rc(3,:)-R1_z,'red','LineWidth',2);
%Acceleration distribution of M2
stem3(R2_x-R1_x,R2_y-R1_y,R2_z-R1_z+MagA2,':^m','color','blue');
grid on;
xlabel('X','Fontsize',18);
ylabel('Y','Fontsize',18);
zlabel('Z','Fontsize',18);
xlim auto;
ylim auto;
zlim auto;
title('Solution at M1 with acceleration distribution of M2','Fontsize',12);
legend('M1 position','M2 position relative to M1','Rc position relative to M1','Acceleration distribution of M2','Location','NorthEast');
%--------------------------------------------------------------------------------------------------------------------------------------------------------
subplot(2,2,4);
view(3);
set(gcf,'Color','w');
hold all;
% M1 position relative to M2 curve
plot3(R1_x-R2_x,R1_y-R2_y,R1_z-R2_z,'color','cyan','LineWidth',2);
%M2 position
plot3(0,0,0,'o','color','g','LineWidth',5);
%Rc position relative to M2 curve
plot3(Rc(1,:)-R2_x,Rc(2,:)-R2_y,Rc(3,:)-R2_z,'red','LineWidth',2);
% Acceleration distribution of M1
stem3(R1_x-R2_x,R1_y-R2_y,R1_z-R2_z+MagA1,':^m','color','blue');
alpha(0.2)
grid on;
xlabel('X','Fontsize',18);
ylabel('Y','Fontsize',18);
zlabel('Z','Fontsize',18);
xlim auto;
ylim auto;
zlim auto;
title('Solution at M2 with acceleration distribution of M1','Fontsize',12);
legend('M1 position relative to M1','M2 position','Rc position relative to M2','Acceleration distribution of M1','Location','NorthEast');
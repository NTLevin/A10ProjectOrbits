file1 = 'data\filter\data one.xlsx';
file2 = 'data\filter\data two.xlsx';
file3 = 'data\filter\data three.xlsx';

STM_file = 'data\filter\state-transition matrix.xlsx';

%t = readtable(file)
F = readmatrix(STM_file);
DATA_one = readmatrix(file1);
DATA_two = readmatrix(file2);
DATA_three = readmatrix(file3);

%Initialisation:

Q = eye(3).*0.01;
H = eye(3);

x_0_0 = DATA_three(:,1);
P_0_0 = eye(3).*100;
R = eye(3)*100;

x_ni_n = F*x_0_0;
P_ni_n = F*P_0_0*(F.') + Q;

N = size(DATA_three,2);

dataKalman = [x_0_0];

for i=2:N
zi = DATA_three(:,i);    

Ki = P_ni_n*(H.')*inv(H*P_ni_n*(H.')+R);

x_ni_ni = x_ni_n + Ki*(zi-H*x_ni_n);
P_ni_ni = (eye(3)-Ki*H)*P_ni_n*(eye(3)-Ki*H).'+Ki*R*(Ki.');


x_n2i_ni = F*x_ni_ni;
P_n2i_ni = F*P_ni_ni*(F.') + Q;

%next iteration assignment
dataKalman = [dataKalman x_n2i_ni];
x_ni_n = x_n2i_ni;
P_ni_n = P_n2i_ni;
end


timestep = 2*pi*sqrt((6378+500)^3/398600)/3600;
time = [];
for j = 1:N
time = [time timestep*j];
end

X1 = DATA_one(1,:);
X2 = dataKalman(1,:);
X3 = DATA_three(1,:);
figure(1)
hold on
plot(time, X3 ,'red','linewidth',2 )
plot(time, X1 ,'green','linewidth',2 )
plot(time, X2 ,'blue','linewidth',2 )
xlabel('Time (sec)');
ylabel('X');
hold off


Y1 = DATA_one(2,:);
Y2 = dataKalman(2,:);
Y3 = DATA_three(2,:);
figure(2)
hold on
plot(time, Y3 ,'red','linewidth',2 )
plot(time, Y1 ,'green','linewidth',2 )
plot(time, Y2 ,'blue','linewidth',2 )
xlabel('Time (sec)');
ylabel('Y');
hold off



Z1 = DATA_one(3,:);
Z2 = dataKalman(3,:);
Z3 = DATA_three(3,:);
figure(3)
hold on
plot(time, Z3 ,'red','linewidth',2 )
plot(time, Z1 ,'green','linewidth',2 )
plot(time, Z2 ,'blue','linewidth',2 )
xlabel('Time (sec)');
ylabel('Z');
hold off

%find slope 
ddx = zeros(length(X1), 1, 'double');
ddx(1) = (X1(1+1)-X1(1))/(time(1+1)-time(1));

ddx(length(X1)) = (X1(length(X1))-X1(length(X1)-1))/ (time(length(X1))-time(length(X1)-1));

for i = 2:(length(X1)-1)
    ddx(i,1)=(X1(i+1)-X1(i-1))/(time(i+1)-time(i-1));
end

%find slope based on average of 5 signal data derivative values 
ddx3 = zeros(length(X3), 1, 'double');
ddx3(1) = (X3(1+1)-X1(1))/(time(1+1)-time(1));

ddx3(length(X3)) = (X3(length(X3))-X3(length(X3)-1))/ (time(length(X3))-time(length(X3)-1));

for i = 2:(length(X3)-1)
    ddx3(i,1)=(X3(i+1)-X3(i-1))/(time(i+1)-time(i-1));
end

ddx3avg = zeros(length(X3), 1, 'double');
for i = 3:(length(X3)-2)
     ddx3avg(i,1)=(ddx3(i-2)+ ddx3(i-1)+ddx3(i)+ddx3(i+1)+ddx3(i+2))/5;
end
ddx3avg(2,1)= (ddx3(1)+ ddx3(2)+ddx3(3))/3;
ddx3avg(1,1)= (ddx3(1)+ ddx3(2))/2;
ddx3avg(length(X3),1)= (ddx3(length(X3))+ ddx3(length(X3)-1))/2;
ddx3avg(length(X3)-1,1)= (ddx3(length(X3))+ ddx3(length(X3)-1)+ddx3(length(X3)-2))/3;






%find slope y axis
ddy = zeros(length(Y1), 1, 'double');
ddy(1) = (Y1(1+1)-Y1(1))/(time(1+1)-time(1));

ddy(length(Y1)) = (Y1(length(Y1))-Y1(length(Y1)-1))/(time(length(Y1))-time(length(Y1)-1));

for i = 2:(length(Y1)-1)
    ddy(i,1)=(Y1(i+1)-Y1(i-1))/(time(i+1)-time(i-1));
end


%find slope based on average of 5 signal data derivative values 
ddy3 = zeros(length(Y3), 1, 'double');
ddy3(1) = (Y3(1+1)-Y1(1))/(time(1+1)-time(1));

ddy3(length(Y3)) = (Y3(length(Y3))-Y3(length(Y3)-1))/ (time(length(Y3))-time(length(Y3)-1));

for i = 2:(length(Y3)-1)
    ddy3(i,1)=(Y3(i+1)-Y3(i-1))/(time(i+1)-time(i-1));
end

ddy3avg = zeros(length(Y3), 1, 'double');
for i = 3:(length(Y3)-2)
     ddy3avg(i,1)=(ddy3(i-2)+ ddy3(i-1)+ddy3(i)+ddy(i+1)+ddy(i+2))/5;
end
ddy3avg(2,1)= (ddy3(1)+ ddy3(2)+ddy3(3))/3;
ddy3avg(1,1)= (ddy3(1)+ ddy3(2))/2;
ddy3avg(length(Y3),1)= (ddy3(length(Y3))+ ddy3(length(Y3)-1))/2;
ddy3avg(length(Y3)-1,1)= (ddy3(length(Y3))+ ddy3(length(Y3)-1)+ddy3(length(Y3)-2))/3;

%find slope y axis
ddz = zeros(length(Z1), 1, 'double');
ddz(1) = (Z1(1+1)-Z1(1))/(time(1+1)-time(1));
 
ddz(length(Z1)) = (Z1(length(Z1))-Z1(length(Z1)-1))/(time(length(Z1))-time(length(Z1)-1));
 
for i = 2:(length(Z1)-1)
    ddz(i,1)=(Z1(i+1)-Z1(i-1))/(time(i+1)-time(i-1));
end
 
 
%find z axis slope based on average of 5 signal data derivative values 
ddz3 = zeros(length(Z3), 1, 'double');
ddz3(1) = (Z3(1+1)-Z1(1))/(time(1+1)-time(1));
 
ddz3(length(Z3)) = (Z3(length(Z3))-Z3(length(Z3)-1))/ (time(length(Z3))-time(length(Z3)-1));
 
for i = 2:(length(Z3)-1)
    ddz3(i,1)=(Z3(i+1)-Z3(i-1))/(time(i+1)-time(i-1));
end
 
ddz3avg = zeros(length(Z3), 1, 'double');
for i = 3:(length(Z3)-2)
     ddz3avg(i,1)=(ddz3(i-2)+ ddz3(i-1)+ddz3(i)+ddz(i+1)+ddz(i+2))/5;
end
ddz3avg(2,1)= (ddz3(1)+ ddz3(2)+ddz3(3))/3;
ddz3avg(1,1)= (ddz3(1)+ ddz3(2))/2;
ddz3avg(length(Z3),1)= (ddz3(length(Z3))+ ddz3(length(Z3)-1))/2;
ddz3avg(length(Z3)-1,1)= (ddz3(length(Z3))+ ddz3(length(Z3)-1)+ddz3(length(Z3)-2))/3;

Z_error_one = dataKalman(3,:)-DATA_one(3,:);
Z_error_three = dataKalman(3,:)-DATA_three(3,:);
figure(12)
plot(time, ddz3avg,'green','linewidth',2 )
hold on
plot(time, Z_error_one ,'blue','linewidth',2 )
xlabel('Time (sec)');
ylabel('Z Error Kalman Filter wrt. true data (km)');
hold off
legend({'Sensor Data averaged derivative over 5 data points','KF Error'})
title('Z Axis Kalman Filter Error (km) vs Averaged Y Axis Signal Derivatives (km/s)')

Y_error_one = dataKalman(2,:)-DATA_one(2,:);
Y_error_three = dataKalman(2,:)-DATA_three(2,:);
figure(11)
plot(time, ddy3avg,'green','linewidth',2 )
hold on
plot(time, Y_error_one ,'blue','linewidth',2 )
xlabel('Time (sec)');
ylabel('Y Error Kalman Filter wrt. true data (km)');
hold off
legend({'Sensor Data averaged derivative over 5 data points','KF Error'})
title('Y Axis Kalman Filter Error (km) vs Averaged Y Axis Signal Derivatives (km/s)')

X_error_one = dataKalman(1,:)-DATA_one(1,:);
X_error_three = dataKalman(1,:)-DATA_three(1,:);
figure(10)
plot(time, ddx3avg,'green','linewidth',2 )
hold on
% plot(time, X_error_three ,'red','linewidth',2 )
plot(time, X_error_one ,'blue','linewidth',2 )
xlabel('Time (sec)');
ylabel('X Error Kalman wrt. true data');
hold off

Y_error_one = dataKalman(2,:)-DATA_one(2,:);
Y_error_three = dataKalman(2,:)-DATA_three(2,:);
figure(9)
plot(time, ddy ,'green','linewidth',2 )
hold on
% plot(time, X_error_three ,'red','linewidth',2 )
plot(time, Y_error_one ,'blue','linewidth',2 )
xlabel('Time (sec)');
ylabel('Y Error Kalman wrt. true data');
hold off

X_error_one = dataKalman(1,:)-DATA_one(1,:);
X_error_three = dataKalman(1,:)-DATA_three(1,:);
figure(8)
plot(time, ddx ,'green','linewidth',2 )
hold on
% plot(time, X_error_three ,'red','linewidth',2 )
plot(time, X_error_one ,'blue','linewidth',2 )
xlabel('Time (sec)');
ylabel('X Error Kalman wrt. true data');
hold off


figure(7)
plot(X_error_one, ddx ,'green','linewidth',2 )
hold on
xlabel('Error (sec)');
ylabel('d/dx orbit');
hold off

figure(6)
plot(dataKalman(2,:)-DATA_one(2,:), ddx2 ,'green','linewidth',2 )
hold on
xlabel('Error (sec)');
ylabel('d/dx orbit');
hold off




Y_error_one = dataKalman(2,:)-DATA_one(2,:);
Y_error_three = dataKalman(2,:)-DATA_three(2,:);
figure(5)
plot(time, Y1 ,'green','linewidth',2 )
hold on
plot(time, Y_error_three ,'red','linewidth',2 )
plot(time, Y_error_one ,'blue','linewidth',2 )
xlabel('Time (sec)');
ylabel('Y Error Kalman wrt. true data');
hold off


Z_error_one = dataKalman(3,:)-DATA_one(3,:);
Z_error_three = dataKalman(3,:)-DATA_three(3,:);
figure(4)
plot(time, Z1 ,'green','linewidth',2 )
hold on
plot(time, Z_error_three ,'red','linewidth',2 )
plot(time, Z_error_one ,'blue','linewidth',2 )
xlabel('Time (sec)');
ylabel('Z Error Kalman wrt. true data');
hold off



X1 = dataKalman(1,:);
X2 = DATA_one(1,:);
X3 = DATA_three(1,:);
figure(1)
plot(time, X3 ,'red','linewidth',2 )
hold all
plot(time, X2 ,'-^g','linewidth',2,'MarkerIndices',1:1:length(X2))
plot(time, X1 ,'-sb','linewidth',2,'MarkerIndices',1:1:length(X1))
xlabel('Time (sec)');
ylabel('J200 X coordinate (km)');
hold off
% legend({'Sensor Data','True Orbit','KF Prediction'},'Location','northwest','Orientation','horizontal')
% title('Orbital Plot of Sensor, True, and Filtered Data of X Axis')

Y1 = dataKalman(2,:);
Y2 = DATA_one(2,:);
Y3 = DATA_three(2,:);
figure(2)
plot(time, Y3 ,'red','linewidth',2 )
hold all
plot(time, Y2 ,'-^g','linewidth',2,'MarkerIndices',1:1:length(Y2))
plot(time, Y1 ,'-sb','linewidth',2,'MarkerIndices',1:1:length(Y1))
xlabel('Time (sec)');
ylabel('J200 Y coordinate (km)');
hold off
% legend({'Sensor Data','True Orbit','KF Prediction'},'Location','northwest','Orientation','horizontal')
% title('Orbital Plot of Sensor, True, and Filtered Data of Y Axis')

Z1 = dataKalman(3,:);
Z2 = DATA_one(3,:);
Z3 = DATA_three(3,:);
figure(3)
plot(time, Z3,'red','linewidth',2);
hold on
plot(time, Z2 ,'-^g','linewidth',2,'MarkerIndices',1:600:length(Z2))
plot(time, Z1 ,'-sb','linewidth',2,'MarkerIndices',1:600:length(Z1))
xlabel('Time (sec)');
ylabel('J2000 Z coordinate (km)');
hold off
% legend({'Sensor Data','True Orbit','KF Prediction'},'Location','northwest','Orientation','horizontal')
% title('Orbital Plot of Sensor, True, and Filtered Data of Z Axis')


SE_X = 0;
SE_X3=0;
E_X=0;
E_X3=0
maxErrorX = 0;
maxErrorX3=0;
i_maxErrorX=0;
i_maxErrorX3=0;
for i=1:length(X1)
    SE_X = SE_X + (X2(i)-X1(i))^2;
    SE_X3=SE_X3 + (X3(i)-X1(i))^2;
    E_X = E_X + (X2(i)-X1(i));
    E_X3 = E_X3 + (X3(i)-X1(i));
    if abs(X1(i)-X2(i)) > abs(maxErrorX)
        maxErrorX = X1(i) - X2(i);
        i_maxErrorX = i;
    end
    if abs(X1(i)-X3(i)) > abs(maxErrorX3)
        maxErrorX3 = X1(i) - X3(i);
        i_maxErrorX3 = i;
    end
end

MSE_X = 1/length(X1)*SE_X 
MSE_X3 = 1/length(X3)*SE_X3 
RMSE_X = sqrt(MSE_X)
RMSE_X3 = sqrt(MSE_X3)
MeanError_X = 1/length(X1)*E_X
MeanError_X3 = 1/length(X3)*E_X3
maxErrorX
maxErrorX3
i_maxErrorX
i_maxErrorX3


SE_Y = 0;
E_Y = 0;
SE_Y3 = 0;
E_Y3 = 0;
maxErrorY = 0;
maxErrorY3 = 0;
i_maxErrorY = 0;
i_maxErrorY3=0;
for i=1:length(Y1)
    SE_Y = SE_Y + (Y2(i)-Y1(i))^2;
    SE_Y3 = SE_Y3 + (Y3(i)-Y1(i))^2;
    
    E_Y = E_Y +(Y2(i)-Y1(i));
    E_Y3 = (Y3(i)-Y1(i))^2;
    if abs(Y1(i)-Y2(i)) > abs(maxErrorY);
        maxErrorY = Y1(i) - Y2(i);
        i_maxErrorY = i;
    end
    if abs(Y1(i)-Y3(i)) > abs(maxErrorY3);
       maxErrorY3 = Y1(i) - Y3(i);
       i_maxErrorY3 = i;
    end
end

%KF error metrics
MSE_Y = 1/length(Y1)*SE_Y 
RMSE_Y = sqrt(MSE_Y)
MeanError_Y = 1/length(Y1)*E_Y
maxErrorY
i_maxErrorY 

%Signal data error metrics
MSE_Y3 = 1/length(Y3)*SE_Y3
RMSE_Y3 = sqrt(MSE_Y3)
MeanError_Y3 = 1/length(Y3)*E_Y3
maxErrorY3
i_maxErrorY3



SE_Z = 0;
SE_Z3 = 0; 
E_Z = 0;
E_Z3 = 0;
maxErrorZ = 0;
maxErrorZ3 = 0;
for i=1:length(Z1)
    SE_Z = SE_Z + (Z2(i)-Z1(i))^2;
    SE_Z3 = SE_Z3 + (Z3(i)-Z1(i))^2;
    E_Z = E_Z +(Z2(i)-Z1(i));
    E_Z3 = E_Z3 + (Z3(i)-Z1(i));
  
    if abs(Z1(i)-Z2(i)) > abs(maxErrorZ)
        maxErrorZ = Z1(i) - Z2(i);
    end
    if abs(Z1(i)-Z3(i)) > abs(maxErrorZ3)
        maxErrorZ3 = Z1(i) - Z3(i);
    end
end

MSE_Z = 1/length(Z1)*SE_Z
MSE_Z3 = 1/length(Z3)*SE_Z3
RMSE_Z = sqrt(MSE_Z)
RMSE_Z3 = sqrt(MSE_Z3)
MeanError_Z = 1/length(Z1)*E_Z
MeanError_Z3 = 1/length(Z3)*E_Z3
maxErrorZ
maxErrorZ3




    
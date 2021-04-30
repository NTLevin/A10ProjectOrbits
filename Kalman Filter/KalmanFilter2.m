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



X_error_one = dataKalman(1,:)-DATA_one(1,:);
X_error_three = dataKalman(1,:)-DATA_three(1,:);
figure(4)
hold on
plot(time, X_error_three ,'red','linewidth',2 )
plot(time, X_error_one ,'blue','linewidth',2 )
xlabel('Time (sec)');
ylabel('Error Kalman wrt. true data');
hold off



X1 = dataKalman(1,:);
X2 = DATA_one(1,:);
figure(1)
plot(time, X1 ,'green','linewidth',2 )
hold all
plot(time, X2 ,'red','linewidth',2 )
xlabel('Time (sec)');
ylabel('J200 X coordinate (km)');
hold off
legend({'KF Prediction','True Orbit'},'Location','northwest','Orientation','horizontal')

Y1 = dataKalman(2,:);
Y2 = DATA_one(2,:);
figure(2)
plot(time, Y1 ,'green','linewidth',2 )
hold all
plot(time, Y2 ,'red','linewidth',2 )
xlabel('Time (sec)');
ylabel('J200 Y coordinate (km)');
hold off
legend({'KF Prediction','True Orbit'},'Location','northwest','Orientation','horizontal')

Z1 = dataKalman(3,:);
Z2 = DATA_one(3,:);
figure(3)
plot(time, Z1 ,'green','linewidth',2)
hold on
plot(time, Z2 ,'red','linewidth',2 )
xlabel('Time (sec)');
ylabel('J2000 Z coordinate (km)');
hold off
legend({'KF Prediction','True Orbit'},'Location','northwest','Orientation','horizontal')
dataSet1 = dataone;
dataSet2 = datatwo;
dataSet3 = datathree;

data = table2array(datathree);
F = table2array(statetransitionmatrix);
x_0_0 = data(:,1);
P_0_0 = [100 0 0; 0 100 0; 0 0 100];
n = width(data);
Q = [100 0 0;0 100 0;0 0 100];
R = Q;

x_n_n = x_0_0;
P_n_n = P_0_0;

datacorrected = [];

for i = 1:n
%Prediction Matrix
x_n1_n = F*x_n_n;
P_n1_n = F*P_n_n*transpose(F)+Q;
%Kalman Gain
K_n1 = P_n1_n*(P_n1_n+R)^(-1);
%Update Matrix
x_n1_n1 = x_n1_n + K_n1*(data(:,(n))-x_n1_n);
P_n1_n1 = (eye(3)- K_n1)*P_n1_n*transpose(eye(3)-K_n1)+K_n1*R*transpose(K_n1);
x_n_n = x_n1_n1;
P_n_n = P_n1_n1;
if i == 1
    x_n1_n1 = x_0_0;  
end 
datacorrected = [datacorrected, x_n1_n1];
end

timestep = 2*pi*sqrt((6378+500)^3/398600)/3600;
time = [];
for j = 1:n
time = [time timestep*j];
end

X1 = datacorrected(1,:);
X2 = data(1,:);
figure(1)
plot(time, X1 ,'green','linewidth',2 )
hold all
plot(time, X2 ,'red','linewidth',2 )
xlabel('Time (sec)');
ylabel('X');

Y1 = datacorrected(2,:);
Y2 = data(2,:);
figure(2)
plot(time, Y1 ,'green','linewidth',2 )
hold all
plot(time, Y2 ,'red','linewidth',2 )
xlabel('Time (sec)');
ylabel('Y');


Z1 = datacorrected(3,:);
Z2 = data(3,:);
figure(3)
plot(time, Z1 ,'green','linewidth',2 )
hold all
plot(time, Z2 ,'red','linewidth',2 )
xlabel('Time (sec)');
ylabel('Z');
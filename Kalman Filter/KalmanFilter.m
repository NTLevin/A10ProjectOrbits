dataSet1 = dataone
dataSet2 = datatwo
dataSet3 = datathree

data = table2array(dataone);
F = table2array(statetransitionmatrix);
x_0_0 = data(:,1);
P_0_0 = [100 0 0; 0 100 0; 0 0 100];
n = width(data);
Q = [100 0 0;0 100 0;0 0 100];
R = Q;



%Prediction Matrix
x_1_0 = F*x_0_0;
P_1_0 = F*P_0_0*rot90(F)+Q
%Kalman Gain
K_1 = P_1_0*(P_1_0+R)^-1;
%Update Matrix
x_1_1 = x_1_0 + K_1*(data(:,2)-x_1_0);
P_1_1 = (eye(3)-K_1)*P_1_0*rot90((eye(3)-K_1))+K_1*R*rot90(K_1);

x_n_n = x_0_0;
P_n_n = P_0_0;

for i = 1:n
%Prediction Matrix
x_n1_n = F*x_n_n;
P_n1_n = F*P_n_n*rot90(F)+1;
%Kalman Gain
K_n1 = P_n1_n*(P_n1_n+R)^(-1);
%Update Matrix
x_n1_n1 = x_n1_n + K_n1*(data(:,(n))-x_n1_n);
P_n1_n1 = (eye(3)- K_n1)*P_n1_n*rot90(eye(3)-K_n1)+K_n1*R*rot90(K_n1);
datacorrected = [datacorrected x_n1_n1];
end

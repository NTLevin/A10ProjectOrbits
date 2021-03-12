A = table2array(PositionVelocity)
% B = table2array(Airpertub);

colRow = size(A);
colLength = colRow(1);
rowLength = colRow(2);

x = A(1:rowLength,2:2)';
y = A(1:rowLength,3:3)';
z = A(1:rowLength,4:4)';

Earth_radius = 6371;
[X,Y,Z] = sphere(20);
X = X*Earth_radius;
Y = Y*Earth_radius;
Z = Z*Earth_radius;
surf(X,Y,Z,'FaceColor','b', 'FaceAlpha',0.2, 'EdgeColor','none');
hold on
figure; plot3(x,y,z,'.-'); 

vx = A(1:rowLength,5:5)';
vy = A(1:rowLength,6:6)';
vz = A(1:rowLength,7:7)';

v = sqrt(A(1:1441,5:5).^2+A(1:1441,6:6).^2+A(1:1441,7:7).^2);


inclination = atan(-1*z/sqrt(x.^2+y.^2));
psi = atan((y/x));
disp(psi);
inclination_v = atan(-1*vz/sqrt(vx.^2+vy.^2));
psi_v = atan(vy/vx);
radius = sqrt(x.^2+y.^2+z.^2);
figure; plot3(vx,vy,vz,'.-');


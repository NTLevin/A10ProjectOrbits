A = PositionVelocity;
%---------------------
B = table2array(A);
size = height(B);
time = A(1:size,1);
t = table2array(time);
x = A(1:size,2);
X = table2array(x)';
y = A(1:size,3);
Y = table2array(y)';
z = A(1:size,4);
Z = table2array(z)';
vx = A(1:size,5);
Vx = table2array(vx)';
vy = A(1:size,6);
Vy = table2array(vy)';
vz = A(1:size,7);
Vz = table2array(vz)';
%CORRECT-------------
absr = sqrt(X.^2 + Y.^2 + Z.^2)';
absv = sqrt(Vx.^2 + Vy.^2 + Vz.^2)';
r = A(1:size,2:4);
r = table2array(r)';
v = A(1:size,5:7);
v = table2array(v)';
h = cross(r(:,1),v(:,1));
for i =2:size
    htemp = cross(r(:,i),v(:,i));
    h = [h;htemp];
end
hsub = sqrt(dot(h,h));
K = [0 0 1]';
%N = cross(h,repmat(K,1,size(h,2))); %wrong as h is wrong
mu = 398600;
E =(absv.^2)/2 - mu * absr.^-1;
%a
a = -0.5*mu*E.^-1;
%e
e1 = hsub.^2/mu^2;
e2 = absv.^2 - 2 * mu * absr.^-1;
e1 = diag(e1);
e3 = e1*e2;
e4 = 1+e3;
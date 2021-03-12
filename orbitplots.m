Earth_radius = 6371;

A = table2array(Airpertub);
x1 = A(1:1441,2);
y1 = A(1:1441,3);
z1 = A(1:1441,4);
figure;plot3(x1,y1,z1,'.-');
axis equal
hold on

[X,Y,Z] = sphere(20);
X = X*Earth_radius;
Y = Y*Earth_radius;
Z = Z*Earth_radius;
surf(X,Y,Z,'FaceColor','b', 'FaceAlpha',0.2, 'EdgeColor','none');
hold on

B = table2array(J2pertub);
x2 = B(1:483,2);
y2 = B(1:483,3);
z2 = B(1:483,4);
figure;plot3(x2,y2,z2,'.-');
axis equal
hold on 

[X,Y,Z] = sphere(20);
X = X*Earth_radius;
Y = Y*Earth_radius;
Z = Z*Earth_radius;
surf(X,Y,Z,'FaceColor','b', 'FaceAlpha',0.2, 'EdgeColor','none');
hold on

C = table2array(Moniya);
x3 = C(1:12962,2);
y3 = C(1:12962,3);
z3 = C(1:12962,4);
figure;plot3(x3,y3,z3,'.-');
axis equal
hold on 

[X,Y,Z] = sphere(20);
X = X*Earth_radius;
Y = Y*Earth_radius;
Z = Z*Earth_radius;
surf(X,Y,Z,'FaceColor','b', 'FaceAlpha',0.2, 'EdgeColor','none');
hold on

global mu;
eps = 1.e-10;


r = norm(R);
v = norm;

vr = dot(R,V)/r;

H = cross(R,V);
h = norm(H);

incl = acos(H(3)/h);

N = cross([0 0 1],H);
n = norm(N);

if n




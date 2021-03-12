Earth_radius = 6371;

A = table2array(Airpertub);
x1 = A(:,2);
y1 = A(:,3);
z1 = A(:,4);
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
x2 = B(:,2);
y2 = B(:,3);
z2 = B(:,4);
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
x3 = C(:,2);
y3 = C(:,3);
z3 = C(:,4);
figure;plot3(x3,y3,z3,'.-');
axis equal
hold on 

[X,Y,Z] = sphere(20);
X = X*Earth_radius;
Y = Y*Earth_radius;
Z = Z*Earth_radius;
surf(X,Y,Z,'FaceColor','b', 'FaceAlpha',0.2, 'EdgeColor','none');
hold on

mu = 398600;
eps = 1.e-10;

R = [-2157.462049 2197.773282 -6150.115340];
V = [-7.167979 -7.036505 -0.000000];

r = norm(R);
v = norm(V);

vr = dot(R,V)./r;

H = cross(R,V);
h = norm(H);

incl = acos(H(3)./h).*180./pi;

N = cross([0 0 1],H);
n = norm(N);

if n ~= 0
    RAAN = acos(N(1)./n);
    if N(2) < 0
        RAAN = 2.*pi - RAAN;
    end
else 
    RAAN = 0;
end
RAAN = RAAN.*180./pi;

E = 1/mu*((v^2-mu/r)* R - r*vr*V);
e = norm(E);

if n ~= 0
    if e > eps
        w = acos(dot(N,E)./n./e);
        if E(3) < 0
            w = 2.*pi - w;
        end
    else 
        w = 0;
    end
else 
    w = 0;
end 
w = w.*180./pi;

if e > eps
    TA = acos(dot(E,R)./e./r);
    if vr < 0
        TA = 2.*pi - TA;
    end
else 
    cp = cross(N,R);
    if cp(3) >= 0
        TA = acos(dot(N,R)./n./r);
    else
        TA = 2.*pi - acos(dot(N,R)./n./r);
    end
end
TA = TA.*180./pi;
while TA >= 360
    TA = TA - 360;
end

a = h.^2./mu./(1-e.^2);

disp(a)
disp(e)
disp(incl)
disp(RAAN)
disp(w)
disp(TA)
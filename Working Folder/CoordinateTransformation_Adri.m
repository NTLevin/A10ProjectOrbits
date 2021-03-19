A = Moniya; %name of input table

titleText = ' Molniya Orbit';
% titleText = ' J2 Perturbation';
% titleText = ' Drag Perturbation';
% Uncomment correct title
    
titleText
%---------------------
B = table2array(A);
size = height(B);
time = A(1:size,1);
time = table2array(time);
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

%Coordinate System - Correct

absr = sqrt(X.^2 + Y.^2 + Z.^2)';
absv = sqrt(Vx.^2 + Vy.^2 + Vz.^2)';
r = A(1:size,2:4);
r = table2array(r)';
v = A(1:size,5:7);
v = table2array(v)';
mu = 398600;

%H 
h = cross(r(:,1),v(:,1))';
for i =2:size
    htemp = cross(r(:,i),v(:,i))';
    h = [h;htemp];
end


absh = sqrt(dot(h(1,:),h(1,:)));
for i = 2:size
    hsubtemp = sqrt(dot(h(i,:),h(i,:)));
    absh = [absh;hsubtemp];
end
K = [0 0 1]';

%N 
N = cross(K,h(1,:));
for i = 2:size
    Ntemp = cross(K,h(i,:));
    N = [N;Ntemp];
end
abs_N = sqrt(dot(N(1,:),N(1,:)));
for i = 2:size
    Nabstemp = sqrt(dot(N(i,:),N(i,:)));
    abs_N= [abs_N;Nabstemp];
end

E =(absv.^2)/2 - mu * absr.^-1;

%a 
a = -0.5*mu*E.^-1;

%e 
e1 = r(:,1).*(absv(1,1)^2-mu/absr(1,:))';
e2 = v(:,1).*dot(r(:,1),v(:,1))';
e = e1/mu - e2/mu;
e = e';
for i=2:size
   e1temp = r(:,i).*(absv(i,1)^2-mu/absr(i,1))';
   e2temp = v(:,i).*dot(r(:,i),v(:,i))';
   etemp = e1temp/mu - e2temp/mu;
   e = [e;etemp'];
end
abse = sqrt(dot(e(1,:),(e(1,:))));
for i=2:size
   absetemp =  sqrt(dot(e(i,:),(e(i,:))));
   abse = [abse;absetemp];
end

%i 
inc = acosd(h(1,3)/absh(1,1));
for i = 2:size
    itemp = acosd(h(i,3)/absh(i,1));
    inc= [inc;itemp];
end

%RAAN 
RAAN = acosd(N(1,1)/abs_N(1,1));
if N(1,2)<0
    RAAN = 360 - RAAN;
end
for i = 2:size
    RAANtemp = acosd(N(i,1)/abs_N(i,1));
    if N(i,2)<0
        RAANtemp  = 360 - RAANtemp;
    end
    RAAN = [RAAN;RAANtemp];
end

%ArgPerigee 
ArgPer1 = dot(N(1,:),e(1,:));
ArgPer2 = ArgPer1/abs_N(1,1)/abse(1,1);
ArgPer = acosd(ArgPer2);
if e(1,3) < 0
    ArgPer = 360 - ArgPer;
end
for i = 2:size
    ArgPer1temp = dot(N(i,:),e(i,:));
    ArgPer2temp = ArgPer1temp/abs_N(i,1)/abse(i,1);
    ArgPertemp = acosd(ArgPer2temp);
    if e(i,3) < 0
        ArgPertemp = 360 - ArgPertemp;
    end
    ArgPer = [ArgPer; ArgPertemp];
end

%TrueAnom
TrueAnom1 = dot(e(1,:),r(:,1)');
TrueAnom2 = TrueAnom1/absr(1,1)/abse(1,1);
TrueAnom = acosd(TrueAnom2);
if dot(r(:,1)',v(:,1)) <0 
    TrueAnom = 360 - TrueAnom;
end
for i=2:size
    TrueAnom1temp = dot(e(i,:),r(:,i)');
    TrueAnom2temp = TrueAnom1temp/absr(i,1)/abse(i,1);
    TrueAnomtemp = acosd(TrueAnom2temp);
    if dot(r(:,i)',v(:,i))<0
        TrueAnomtemp = 360 - TrueAnomtemp;
    end
    TrueAnom = [TrueAnom;TrueAnomtemp];
end

%time; a; abse; inc; RAAN; ArgPer; TrueAnom; 
Time = array2table(time);
Semi_Major_Axis = array2table(a);
Eccentricity = array2table(abse);
Inclination = array2table(inc);
RAAN = array2table(RAAN);
Arg_of_Perigee = array2table(ArgPer);
True_Anomaly = array2table(TrueAnom);
FinalTemp = [Time Semi_Major_Axis Eccentricity Inclination RAAN Arg_of_Perigee True_Anomaly];

clear A asbetemp absh absr absv ArgPer1 ArgPer1temp ArgPer2 ArgPer2temp ArgPertemp B e E 
clear e1temp absetemp e2 e2temp etemp hsubtemp htemp i itemp K mu e1
clear r RAANtemp size t TrueAnom1 TrueAnom1temp TrueAnom2 TrueAnom2temp TrueAnomtemp
clear v vx Vx vy Vy vz Vz x X y Y z Z trueanom argper
clear time a abse inc RAAN ArgPer TrueAnom abs_N h N Nabstemp Ntemp
clear Time Semi_Major_Axis Eccentricity Inclination RAAN Arg_of_Perigee True_Anomaly

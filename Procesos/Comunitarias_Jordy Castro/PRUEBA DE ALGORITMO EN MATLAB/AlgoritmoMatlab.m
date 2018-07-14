
clc;
clear all;
close all;
format short
syms q1 q2 beta1 beta2 l1 l2 alpha1 alpha2 real
 
H20=H_r2gdl(); 
 
disp('Transformación homogénea del robot 2gdl');
disp(H20);
 
[R20, cinemat_r2gdl,cero, c]=H_DH(H20);
disp('Matriz de rotación'); disp(R20);
disp('cinemática directa');
disp(cinemat_r2gdl);
[x0, y0, z0] = cinematica_r2gdl(beta1,l1,q1,beta2,l2,q2);
jac_r2gdl=jacobian([x0; y0], [q1;q2]);
det_r2gdl=simplify(det(jac_r2gdl));
% det[J]=l_1l_2 sin(q_2)
 
%ejemplo numérico 
%t=0:0.0001:4.7122;
t=0:0.0001:6.2830;
 
%parámetros del círculo:[x_c,y_c]'=[0.3,-0.3]' y radio r=0.2
xc=0.3; yc=-0.3; r=0.20;
l1=0.45; l2=0.45;
beta1=0.1; beta2=0.1;
q1=[]; q2=[];
 
% ecuación del círculo
x=xc+r*sin(t);
%x=xc+r*sin(2*t);
 
%y=yc+r*cos(t);
y=yc+r*cos(t)+r*sin(t);
 
% cinemática inversa
[q1,q2]=cinv_r2gdl(l1,l2,x,y);
%coordenas cartesianas del extremo final del robot de 2 gdl
[x0, y0, z0] = cinematica_r2gdl(beta1,l1,q1,beta2,l2, q2);
figure 
plot(x0,y0,'LineWidth',2,'color',[0.8,0.1,0.1]) 
El código 1, hace llamadas a otras funciones, las cuales son:


function H=H_r2gdl()
syms q1 q2 beta1 beta2 l1 l2 alpha1 alpha2 real
disp('Parámetros Denavit-Hartenberg del robot planar vertical de 2 gdl')
disp('[ l alpha d q]')
dh=[l1, 0, beta1, q1; l2, 0, beta2, q2];
disp(dh)
%$H10=HDH{0}{q_1}{beta_1}{l_1}{0}
H10=HRz(q1)*HTz(beta1)*HTx(l1)*HRx(0);
%$H21=HDH{0}{q_2}{beta_2}{l_2}{0}
H21=HRz(q2)*HTz(beta2)*HTx(l2)*HRx(0);
H20=simplify(H10*H21); %H20=H10 H21
[R20, cinemat_r2gdl, cero, c]=H_DH(H20);
 H=[R20, cinemat_r2gdl; %R20(q_1,q_2),f_R(q_1,q_2)
 cero, c];
end

function [q1,q2]=cinv_r2gdl(l1,l2,x0,y0)
q2=acos((x0.*x0+y0.*y0-l1*l1-l2*l2)/(2*l1*l2));
q1=atan(y0./x0)-atan((l2*sin(q2))./(l1+l2*cos(q2)));
end


function [x0, y0, z0]=cinematica_r2gdl(beta1,l1,q1,beta2,l2,q2)
dato1=whos('beta1'); dato2=whos('l1');
dato3=whos('q1');
dato4=whos('beta2'); dato5=whos('l2');
dato6=whos('q2');
v1=strcmp(dato1.class, 'sym');
v2=strcmp(dato2.class, 'sym');
v3=strcmp(dato3.class, 'sym');
v4=strcmp(dato4.class, 'sym');
v5=strcmp(dato5.class, 'sym');
v6=strcmp(dato6.class, 'sym');
digits(3);
 if (v1 & v2 & v3 & v4 & v5 & v6) %caso simbólico
  x0=simplify(vpa(l1*cos(q1)+l2*cos(q1+q2),3));
  y0=simplify(vpa(l1*sin(q1)+l2*sin(q1+q2),3));
  z0=vpa(beta1+beta2,3);
 else %caso numérico
  x0=l1*cos(q1)+l2*cos(q1+q2);
  y0=l1*sin(q1)+l2*sin(q1+q2);
  z0=beta1+beta2;
 end
end


function RHz=HRz(theta)
dato=whos('theta');
if strcmp(dato.class, 'sym') %variables simbólicas
 RHz=[cos(theta), -sin(theta), 0, 0;
 sin(theta), cos(theta), 0, 0;
 0, 0, 1, 0;
 0, 0, 0, 1];
else
 digits(3); %cálculos numéricos
 RHz=round([ vpa(cos(theta),3), vpa(-sin(theta),3), 0, 0;
 vpa(sin(theta),3),vpa(cos(theta),3), 0, 0;
 0, 0, 1, 0;
 0, 0, 0, 1]);
end
end

%Función HTz para el cálculo de la variable simbólica d.
function Tz=HTz(d)
Tz=[ 1 0 0 0; 0 1 0 0;
0 0 1 d; 0 0 0 1];
end

%Función HTx para el cálculo de la variable simbólica d.
function Tx=HTx(d)
Tx=[ 1 0 0 d; 0 1 0 0;
0 0 1 0; 0 0 0 1];
end

%Función HTy para el cálculo de la variable simbólica d.
function Tz=HTy(d)
Tz=[ 1 0 0 0; 0 1 0 d;
0 0 1 0; 0 0 0 1];
end

%Función HRx para el cálculo de la variable simbólica theta. 
function RHx=HRx(theta)
dato=whos('theta');
if strcmp(dato.class, 'sym') %variables simbólicas
RHx=[1, 0, 0, 0;
 0, cos(theta), -sin(theta), 0;
 0, sin(theta), cos(theta), 0;
 0, 0, 0, 1];
else
digits(3); %cálculos numéricos
 RHx=round([1, 0, 0, 0;
 0, vpa(cos(theta),3), vpa(-sin(theta),3), 0;
 0, vpa(sin(theta),3),vpa(cos(theta),3), 0;
 0, 0, 0, 1]);
end
end

%Función HRy para el cálculo de la variable simbólica theta.
function RHy=HRy(theta)
dato=whos('theta');
if strcmp(dato.class, 'sym') %variables simbólicas
RHy=[cos(theta), 0, sin(theta), 0;
 0, 1, 0, 0;
 -sin(theta), 0, cos(theta), 0;
 0, 0, 0, 1];
else
 digits(3); %cálculos numéricos
 RHy=round([ vpa(cos(theta),3), 0, vpa(sin(theta),3), 0;
 0, 1, 0, 0;
 vpa(-sin(theta),3), 0, vpa(cos(theta),3), 0;
 0, 0, 0, 1]);
end
end

%Función H_DH para el cálculo de la matriz de transformación homogénea.
function [R vect_d vect_cero c]=H_DH(H)
for i=1:3
   for j=1:3 R(i,j)=H(i,j);
   end
end
%estructura de la matriz de transformación homogénea
vect_d=[H(1,4); H(2,4); H(3,4)];
vect_cero=[0;0;0]';
c=1;
end

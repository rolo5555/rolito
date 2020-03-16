%% definicion de los parametros de mi mecanismo 
clear all,clc,close all;
a= 350 %[mm]
b= 2000 %[mm]
c= 1280 %[mm]
d= 2000 %[mm]
%% resolucion 
n_steps=500;
teta2=linspace(0,2*pi,n_steps);
K1=(d/a),K2=(d/c),K3=(a^2-b^2+c^2+d^2)/(2*a*c);
K4=(d/b),K5=(c^2-d^2-a^2-b^2)/(2*a*b);
for i=1:n_steps
A=cos(teta2(i))-K1-K2*cos(teta2(i))+K3;
B=-2*sin(teta2(i));
C=K1-(K2+1)*cos(teta2(i))+K3;
teta4(i)=2*atan((-B+sqrt(B^2-4*A*C))/(2*A));
D=cos(teta2(i))-K1+K4*cos(teta2(i))+K5;
E=-2*sin(teta2(i)); 
F=K1+(K4-1)*cos(teta2(i))+K5;
teta3(i)=2*atan((-E+sqrt(E^2-4*D*F))/(2*D));
end
%% chorresion del angulo
ang_add=pi-atan(1600/1200);
for i=1:n_steps
teta2(i)=teta2(i)+2*pi-ang_add;
teta3(i)=teta3(i)+2*pi-ang_add;
teta4(i)=teta4(i)+2*pi-ang_add;
end
%% convierto a grados
teta2=rad2deg(teta2);
teta3=rad2deg(teta3);
teta4=rad2deg(teta4);
%% Angulos de extremos de 
mu_1=acos((b^2+c^2-(d+a)^2)/(2*b*c));
mu_2=acos((b^2+c^2-(d-a)^2)/(2*b*c));
mu_1=rad2deg(mu_1);
mu_2=rad2deg(mu_2);
%% gráficos
figure(1)
hold on
title('Angulo theta 3')
grid on
xlabel('theta2')
ylabel('theta3')
plot(teta2,teta3)  
hold off

figure(2)
hold on
title('Angulo theta 4')
grid on
xlabel('theta2')
ylabel('theta4')
plot(teta2,teta4)  
hold off
%% analisis de movimiento 
O2= [0 0]
axis(gca,'equal');
axis([-3000 3000 -3000 3000]);
O4=[-1200 1600]
for i=1:n_steps
    punto_A= a*[cosd(teta2(i)) sind(teta2(i))];
    crank = line([O2(1),punto_A(1)],[O2(2),punto_A(2)]);grid on;
    %intervalo de actualizacion
    %punto_B= punto_A+b*[cosd(teta4(i)) sind(teta4(i))]
    %biela= line([punto_A(1) punto_B(1)],[punto_A(2) punto_B(2)])
    %balancin=line([punto_B(1) O4(1)],[punto_B(2) O4(2)])
    %mod_1(i)=sqrt((punto_A(1)-O2(1))^2+(punto_A(2)-O2(2))^2);
    %mod_2(i)=sqrt((punto_A(1)-punto_B(1))^2+(punto_A(2)-punto_B(2))^2);
    %mod_3(i)=sqrt((punto_B(1)-O4(1))^2+(punto_B(2)-O4(2))^2);
    pause(0.01);       
    %BORRO LO PREVIO
    delete(crank);
    %delete(biela);
    %delete(balancin);
end



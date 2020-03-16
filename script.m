%% definicion de los parametros de mi mecanismo 
clear all,clc,close all;
a= 350 %[mm]
b= 2000 %[mm]
c= 1280 %[mm]
d= 2000 %[mm] 
%% otra manera de resolverlo
%r=[b, a, c, d];
n_steps=500;
teta1=pi-atan(1600/1200);
teta2=linspace(0,2*pi,n_steps);
teta0=[0 0];
for i=1:n_steps
xxx=fsolve(@(teta) myfun(teta,teta2(i),teta1),teta0);
x(i,1)=xxx(1);
x(i,2)=xxx(2);
end
%% Angulos de extremos de mi mecanismo
mu_1=acos((b^2+c^2-(d+a)^2)/(2*b*c));
mu_2=acos((b^2+c^2-(d-a)^2)/(2*b*c));
mu_1=rad2deg(mu_1);
mu_2=rad2deg(mu_2);
%% gráficos de los angulos 
figure(1)
hold on
title('Angulo theta 3')
grid on
xlabel('theta2 [°]')
ylabel('theta3 [°]')
plot(rad2deg(teta2),rad2deg(x(:,1)))  
legend('theta3')
hold off
figure(2)
hold on
title('Angulo theta 4')
grid on
xlabel('theta2 [°]')
ylabel('theta4 [°]')
plot(rad2deg(teta2),rad2deg(x(:,2)))
legend('theta4')
hold off
%% analisis de movimiento animacion
O2= [0 0]
axis(gca,'equal');
axis([-3000 3000 -3000 3000]);
O4=[-1200 1600]
for i=1:n_steps
    punto_A= a*[cos(teta2(i)) sin(teta2(i))];
    crank = line([O2(1),punto_A(1)],[O2(2),punto_A(2)]);grid on;
    %intervalo de actualizacion
    punto_B= punto_A+b*[cos(x(i,1)) sin(x(i,1))]
    biela= line([punto_A(1) punto_B(1)],[punto_A(2) punto_B(2)])
    balancin=line([punto_B(1) O4(1)],[punto_B(2) O4(2)])
    bancada= line([O4(1) O2(1)],[O4(2) O2(2)])
    %mod_1(i)=sqrt((punto_A(1)-O2(1))^2+(punto_A(2)-O2(2))^2);
    %mod_2(i)=sqrt((punto_A(1)-punto_B(1))^2+(punto_A(2)-punto_B(2))^2);
    %mod_3(i)=sqrt((-punto_B(1)+O4(1))^2+(-punto_B(2)+O4(2))^2);
    pause(0.001);       
    %BORRO LO PREVIO
    delete(crank);
    delete(biela);
    delete(balancin);
end
%% posiciones de agarrotamiento 
clear dif_2_3;clear dif_3_4;
for i=1:n_steps   
   dif_2_3(i)=rem(abs(teta2(i)-x(i,1)),pi);
   dif_3_4(i)=rem(abs(x(i,2)-x(i,1)-pi),pi);
end
%% busco los minimos de los vectores
[value_agarrot(1),pos_agarrot_i(1)]=min(dif_2_3(1:300))
[value_agarrot(2),pos_agarrot_i(2)]=min(dif_2_3(200:500))
teta_agarrot_1=rad2deg(teta2(pos_agarrot_i(1)))
teta_agarrot_2=rad2deg(teta2(pos_agarrot_i(2)))
%%
figure(1)
hold on
title('Posiciones de agarrotamiento')
plot(rad2deg(teta2),rad2deg(dif_2_3))
%plot(rad2deg(teta2),dif_3_4)
xlabel('Theta 2 [°]')
ylabel('diferencia Theta2-Theta3 [°]')
hold off
%% coordenadas de los puntos A y B
for i=1:n_steps
    punto_A(i,1)= 350*cos(teta2(i));
    punto_A(i,2)= 350*sin(teta2(i));
    punto_B(i,1)= punto_A(i,1)+2000*cos(x(i,1)); 
    punto_B(i,2)= punto_A(i,2)+2000*sin(x(i,1)); 
    %aux_angle(i)=x(i,2)-pi-rad2deg(14)
end
hold on 
title('Coordenada Y del punto B')
plot(rad2deg(teta2),punto_B(:,2))
grid on
xlabel('theta2 [°]')
ylabel('Coordena Y del punto B [mm]')
hold off
[M,I] = min(punto_B(:,2))
[M,I]= max(punto_B(:,2))
%% grafico de las coordenadas de las juntas A y B
hold on
grid on
title('Coordenadas de las juntas A y B')
plot(punto_A(:,1),punto_A(:,2))
plot(punto_B(:,1),punto_B(:,2))
xlabel('Coordenada X [mm]')
ylabel('Coordenada Y [mm]')
legend('Punto A','Punto B')
hold off
%% posiciones de los centros de gravedad
long_O4_CG4=1980;
O4=[-1200 1600];
for i=1:n_steps
    CG_2(i,1) =(350/2)*cos(teta2(i));
    CG_2(i,2) =(350/2)*sin(teta2(i));
    punto_A= 350*[cos(teta2(i)) sin(teta2(i))];
    CG_3(i,1)= punto_A(1) + (2000/2)*cos(x(i,1));
    CG_3(i,2)= punto_A(2) + (2000/2)*sin(x(i,1));
    punto_B= punto_A+2000*[cos(x(i,1)) sin(x(i,1))]; 
    aux_angle(i)=x(i,2)-deg2rad(14);
    CG_4(i,1)= O4(1)+ 1980*cos(aux_angle(i));
    CG_4(i,2)= O4(2)+ 1980*sin(aux_angle(i));
end
%%
hold on
title('Coordenadas de los centros de gravedad')
plot(CG_2(:,1),CG_2(:,2))
plot(CG_3(:,1),CG_3(:,2))
plot(CG_4(:,1),CG_4(:,2))
xlabel('Coordenada X [mm]')
ylabel('Coordenada Y [mm]')
legend('CG_2','CG_3','CG_4')
grid on
hold off
%% angulo de avance para P
plot(teta2,punto_B(:,1))
%% verificacion con metodo grafico
 %angulo de 30   
    punto_A= a*[cos(teta2(43)) sin(teta2(43))];
    crank = line([O2(1),punto_A(1)],[O2(2),punto_A(2)]);grid on;
    punto_B= punto_A+b*[cos(x(43,1)) sin(x(43,1))]
    biela= line([punto_A(1) punto_B(1)],[punto_A(2) punto_B(2)])
    balancin=line([punto_B(1) O4(1)],[punto_B(2) O4(2)])
    bancada= line([O4(1) O2(1)],[O4(2) O2(2)])
    
    show1=rad2deg(x(43,1))
    show2=rad2deg(x(43,2))
%% verificacion con metodo grafico
 %angulo de 120   
    punto_A= a*[cos(teta2(167)) sin(teta2(167))];
    crank = line([O2(1),punto_A(1)],[O2(2),punto_A(2)]);grid on;
    punto_B= punto_A+b*[cos(x(167,1)) sin(x(167,1))]
    biela= line([punto_A(1) punto_B(1)],[punto_A(2) punto_B(2)])
    balancin=line([punto_B(1) O4(1)],[punto_B(2) O4(2)])
    bancada= line([O4(1) O2(1)],[O4(2) O2(2)])
 
    show1=rad2deg(x(167,1))
    show2=rad2deg(x(167,2))
%% PARTE DE LA VELOCIDAD
%%velocidad de omega 3 y omega 4 en funcion de omega 2
clc;
omega2=(4*2*pi/60); 
omega0=[0 0];
for i=1:n_steps
vector_aux = fsolve(@(omega) velocity(omega,x(i,:),teta2(i),omega2),omega0); %teta,teta1,teta2,omega2,omega
omega_ans(i,1)=vector_aux(1);
omega_ans(i,2)=vector_aux(2);
end 
%% grafico de las velocidades angulares de los eslabones
hold on
plot(teta2,omega_ans(:,1))
plot(teta2,omega_ans(:,2))
title('Velocidades angulares en funcion de \theta_{2}','FontSize',11)
xlabel('Angulo \theta_{2} [rad]','FontSize',11)
ylabel('Velocidad angular [rad/s]','FontSize',11)
legend('\omega_{3}','\omega_{4}','FontSize',11)
grid on
hold off
%% velocidades juntas
a= 350 %[mm]
b= 2000 %[mm]
c= 1280 %[mm]
d= 2000 %[mm] 
for i=1:n_steps
V_Ax(i)= -omega2*a*sin(teta2(i));
V_Ay(i)= omega2*a*cos(teta2(i));
V_A_mod(i)=sqrt(V_Ax(i)^2+V_Ay(i)^2);
%V_Bx(i)= V_Ax(i)+ b*omega_ans(i,1)*cos(x(i,1)) 
%V_By(i)= V_Ay(i)+ -b*omega_ans(i,1)*sin(x(i,1))
V_Bx1(i)= -c*sin(x(i,2))*omega_ans(i,2);
V_By1(i)= c*cos(x(i,2))*omega_ans(i,2);
V_B_mod(i)=sqrt(V_Bx1(i)^2+V_By1(i)^2);
end
%% velocidades de las juntas A y B
% de la junta A
hold on
plot(teta2,V_Ax)
plot(teta2,V_Ay)
%plot(teta2,omega_ans(:,2))
title('Velocidad de la junta A','FontSize',11)
xlabel('Angulo \theta_{2} [rad]','FontSize',11)
ylabel('Velocidad [mm/s]','FontSize',11)
legend('V_Ax','V_Ay','FontSize',11)
grid on
hold off
%% comprobacion de los resultados de las cvelocidades
%para la junta B
hold on 
plot(teta2,V_Bx1)
plot(teta2,V_By1)
title('Velocidad de la junta B','FontSize',11)
xlabel('Angulo \theta_{2} [rad]','FontSize',11)
ylabel('Velocidad [mm/s]','FontSize',11)
legend('V_{Bx}','V_{By}','FontSize',11)
grid on
hold off
%% comprobacion de los resultados de las cvelocidades
%para los modulos
hold on 
plot(teta2,V_A_mod)
plot(teta2,V_B_mod,'LineWidth',1.5)
title('Velocidad de la juntas','FontSize',11)
xlabel('Angulo \theta_{2} [rad]','FontSize',11)
ylabel('Velocidad [mm/s]','FontSize',11)
legend('V_{B}','V_{A}','FontSize',11)
hold off
%% velocidades de los centros de gravedad 
for i=1:n_steps
V_CG2_x(i)= omega2*-330*sin(teta2(i));
V_CG2_y(i)= omega2*330*cos(teta2(i));
V_CG3_x(i)= V_Ax(i)-(b/2)*sin(x(i,1))*omega_ans(i,1);
V_CG3_y(i)= V_Ay(i)+(b/2)*cos(x(i,1))*omega_ans(i,1);
V_CG4_x(i)= -1980*sin(x(i,2)-deg2rad(14))*omega_ans(i,2);
V_CG4_y(i)= 1980*cos(x(i,2)-deg2rad(14))*omega_ans(i,2);
end
%% graficos de las velocidades de los centros de gravedad
%para todos los centros de gravedad
%PARA CG2
hold on 
plot(teta2,V_CG2_x)
plot(teta2,V_CG2_y)
title('Velocidad del CG2','FontSize',11)
xlabel('Angulo \theta_{2} [rad]','FontSize',11)
ylabel('Velocidad [mm/s]','FontSize',11)
legend('VCG2_{x}','VCG2_{y}','FontSize',11)
grid on
hold off
%% PARA CG3 
hold on 
plot(teta2,V_CG3_x)
plot(teta2,V_CG3_y)
title('Velocidad del CG3','FontSize',11)
xlabel('Angulo \theta_{2} [rad]','FontSize',11)
ylabel('Velocidad [mm/s]','FontSize',11)
legend('VCG3_{x}','VCG3_{y}','FontSize',11)
grid on
hold off
%% PARA CG4 
hold on 
plot(teta2,V_CG4_x)
plot(teta2,V_CG4_y)
title('Velocidad del CG4','FontSize',11)
xlabel('Angulo \theta_{2} [rad]','FontSize',11)
ylabel('Velocidad [mm/s]','FontSize',11)
legend('VCG4_{x}','VCG4_{y}','FontSize',11)
grid on
hold off
%% PARTE DE LA ACELERACION
clc;
alfa2=0; 
alfa0=[0 0];
for i=1:n_steps
A=c*sin(x(i,2))
B=b*sin(x(i,1))
C=a*alfa2*sin(teta2(i))+a*omega2^2*cos(teta2(i))+b*(omega_ans(i,1))^2*cos(x(i,1))-c*(omega_ans(i,2))^2*cos(x(i,2))
D=c*cos(x(i,2))
E=b*cos(x(i,1))
F=a*alfa2*cos(teta2(i))-a*omega2^2*sin(teta2(i))-b*(omega_ans(i,1))^2*sin(x(i,1))+c*(omega_ans(i,2))^2*sin(x(i,2))
%vector_aux = fsolve(@(alfa) aceleration(alfa,alfa2,omega_ans,x(i,:),teta2(i),omega2),alfa0,options); %teta,teta1,teta2,omega2,omega
alfa_ans(i,1)=(C*D-A*F)/(A*E-B*D);
alfa_ans(i,2)=(C*E-B*F)/(A*E-B*D);
end 
%% ploteo de las aceleraciones angulares de los eslabones
hold on
plot(teta2,alfa_ans(:,1))
plot(teta2,alfa_ans(:,2))
title('Aceleraciones angulares en funcion de \theta_{2}','FontSize',11)
xlabel('Angulo \theta_{2} [rad]','FontSize',11)
ylabel('Aceleración angular [rad/s^{2}]','FontSize',11)
legend('\alpha_{3}','\alpha_{4}','FontSize',11)
grid on
hold off
%% aceleracion de las juntas A y B 
for i=1:n_steps
A_Ax(i)= -omega2^2*a*cos(teta2(i));
A_Ay(i)= -omega2^2*a*sin(teta2(i));
A_A_mod(i)=sqrt(A_Ax(i)^2+A_Ay(i)^2);
A_Bx(i)= cos(x(i,2))*(-c)*omega_ans(i,2)^2+cos((pi/2)-x(i,2))*-c*alfa_ans(i,2);
A_By(i)= sin(x(i,2))*-c*omega_ans(i,2)^2+sin((pi/2)-x(i,2))*c*alfa_ans(i,2);
A_B_mod(i)=sqrt(A_By(i)^2+A_Bx(i)^2);
end
%% ploteo de las aceleracion de las junta a
hold on
plot(teta2,A_Ax)
plot(teta2,A_Ay)
title('Aceleraciones de vinculo A en funcion de \theta_{2}','FontSize',11)
xlabel('Angulo \theta_{2} [rad]','FontSize',11)
ylabel('Aceleración [mm/s^{2}]','FontSize',11)
legend('Aa_{x}','Aa_{y}','FontSize',11)
grid on
hold off
%% ploteo de las aceleracion de la junta b
hold on
plot(teta2,A_Bx)
plot(teta2,A_By)
title('Aceleraciones de vinculo B en funcion de \theta_{2}','FontSize',11)
xlabel('Angulo \theta_{2} [rad]','FontSize',11)
ylabel('Aceleración [mm/s^{2}]','FontSize',11)
legend('Ab_{x}','Ab_{y}','FontSize',11)
grid on
hold off
%% Aceleracion de los centros de gravedad 
for i=1:n_steps
A_CG2_x(i)= omega2^2*-330*cos(teta2(i));
A_CG2_y(i)= -omega2^2*330*sin(teta2(i));
A_CG3_x(i)= A_Ax(i)-(2000/2)*cos(x(i,1))*omega_ans(i,1)^2-(2000/2)*alfa_ans(i,1)*sin(x(i,1));
A_CG3_y(i)= A_Ay(i)-(2000/2)*sin(x(i,1))*omega_ans(i,1)^2+(2000/2)*alfa_ans(i,1)*cos(x(i,1));
A_CG3_mod(i)= sqrt(A_CG3_x(i)^2+A_CG3_y(i)^2);
A_CG4_x(i)= -1980*omega_ans(i,2)^2*cos(x(i,2)-deg2rad(14))+-alfa_ans(i,2)*1980*sin(x(i,2)-deg2rad(14));
A_CG4_y(i)= -1980*omega_ans(i,2)^2*sin(x(i,2)-deg2rad(14))+alfa_ans(i,2)*1980*cos(x(i,2)-deg2rad(14));
A_CG4_mod(i)= sqrt(A_CG4_x(i)^2+A_CG4_y(i)^2);
end
%% graficos de las aceleraciones de los centros de gravedad
%para todos los centros de gravedad
%PARA CG2
hold on 
plot(teta2,A_CG2_x)
plot(teta2,A_CG2_y)
title('Aceleracion del CG2','FontSize',11)
xlabel('Angulo \theta_{2} [rad]','FontSize',11)
ylabel('Aceleracion [mm^2/s]','FontSize',11)
legend('ACG2_{x}','ACG2_{y}','FontSize',11)
grid on
hold off
%% PARA CG3
hold on 
plot(teta2,A_CG3_x)
plot(teta2,A_CG3_y)
plot(teta2,A_CG3_mod)
title('Aceleracion del CG3','FontSize',11)
xlabel('Angulo \theta_{2} [rad]','FontSize',11)
ylabel('Aceleracion [mm^2/s]','FontSize',11)
legend('ACG3_{x}','ACG3_{y}','ACG3_{mod}','FontSize',11)
grid on
hold off
%% PARA CG4
hold on 
plot(teta2,A_CG4_x)
plot(teta2,A_CG4_y)
plot(teta2,A_CG4_mod)
title('Aceleracion del CG4','FontSize',11)
xlabel('Angulo \theta_{2} [rad]','FontSize',11)
ylabel('Aceleracion [mm^2/s]','FontSize',11)
legend('ACG4_{x}','ACG4_{y}','ACG4_{mod}','FontSize',11)
grid on
hold off
%% PARTE DE FUERZAS 
%% definicion de los parametros de mi sistema 
g=9.8; %[m/s^2]
m2=2660/g; %[N]
m3=480/g; %[N]
m4=12000/g; %[N]
icg_2=1.33; %[N.m/s^2] 
icg_3=17; %[N.m/s^2] 
icg_4=1210; %[N.m/s^2] 
r12=(330/1000);r32=(20/1000); %[m]
r23=(1000/1000);r43=(1000/1000); %[m] 
O4_metro=[-1200 1600]/1000;punto_B_m=punto_B/1000;CG_4_m=CG_4/1000; % todo en metro como buen cristiano
A_CG2_x_m=A_CG2_x/1000;A_CG2_y_m=A_CG2_y/1000;A_CG3_x_m=A_CG3_x/1000;A_CG3_y_m=A_CG3_y/1000;A_CG4_x_m=A_CG4_x/1000;A_CG4_y_m=A_CG4_y/1000; % todo en m/s^2 como buen cristiano
r_torque_4=-2.4-CG_4_m(:,1);
%% definicion de la fuerzas en carrera ascendete y descendente
for i=1:n_steps
if(131<i && i<371)
    f_as_des(i)=-13200; %[N]
else f_as_des(i)=-10200; %[N]
end 
end 
%% ahora defino a mi sitema en forma matricial
for i=1:n_steps 
 fila1=[1 0 1 0 0 0 0 0 0];
 fila2=[0 1 0 1 0 0 0 0 0];
 fila3=[r12*sin(teta2(i)) -r12*cos(teta2(i)) -r32*sin(teta2(i)) r32*cos(teta2(i)) 0 0 0 0 1];
 fila4=[0 0 -1 0 1 0 0 0 0];
 fila5=[0 0 0 -1 0 1 0 0 0];
 fila6=[0 0 -r23*sin(x(i,1)) r23*cos(x(i,1)) -r43*sin(x(i,1)) r43*cos(x(i,1)) 0 0 0];
 fila7=[0 0 0 0 -1 0 1 0 0];
 fila8=[0 0 0 0 0 -1 0 1 0];
 fila9=[0 0 0 0 (punto_B_m(i,2)-CG_4_m(i,2)) -(punto_B_m(i,1)-CG_4_m(i,1)) -(-CG_4_m(i,2)+O4_metro(2)) (-CG_4_m(i,1)+O4_metro(1)) 0]
 A=[fila1;fila2;fila3;fila4;fila5;fila6;fila7;fila8;fila9];
 b1=m2*A_CG2_x_m(i);
 b2=m2*A_CG2_y_m(i)+m2*g;
 b3=icg_2*alfa2;
 b4=m3*A_CG3_x_m(i);
 b5=m3*A_CG3_y_m(i)+m3*g;
 b6=icg_3*alfa_ans(i,1);
 b7=m4*A_CG4_x_m(i);
 b8=m4*A_CG4_y_m(i)+m4*g-f_as_des(i);
 b9=icg_4*alfa_ans(i,2)-f_as_des(i)*r_torque_4(i);
 b=[b1;b2;b3;b4;b5;b6;b7;b8;b9];
 x_fuerzas_aux=inv(A)*b;
 x_fuerzas(1,i)=x_fuerzas_aux(1);
 x_fuerzas(2,i)=x_fuerzas_aux(2);
 x_fuerzas(3,i)=x_fuerzas_aux(3);
 x_fuerzas(4,i)=x_fuerzas_aux(4);
 x_fuerzas(5,i)=x_fuerzas_aux(5);
 x_fuerzas(6,i)=x_fuerzas_aux(6);
 x_fuerzas(7,i)=x_fuerzas_aux(7);
 x_fuerzas(8,i)=x_fuerzas_aux(8);
 x_fuerzas(9,i)=x_fuerzas_aux(9);
end
%% calculo de las normas de las fuerzas en la juntas 
for i=1:n_steps
f32_A_norma(i)=norm([x_fuerzas(3,i),x_fuerzas(4,i)]);
f43_A_norma(i)=norm([x_fuerzas(5,i),x_fuerzas(6,i)]);
f_O2_norma(i)=norm([x_fuerzas(1,i),x_fuerzas(2,i)]);
f_O4_norma(i)=norm([x_fuerzas(7,i),x_fuerzas(8,i)]);
end 
%% graficos de las fuerzas en las juntas PARA A
hold on 
title('Fuerza en la junta A','FontSize',11)
subplot(311)
plot(teta2,x_fuerzas(3,:))
xlabel('Angulo \theta_{2} [rad]','FontSize',11)
ylabel('FA_{x} [N]','FontSize',11)
subplot(312)
plot(teta2,x_fuerzas(4,:))
xlabel('Angulo \theta_{2} [rad]','FontSize',11)
ylabel('FA_{y} [N]','FontSize',11)
subplot(313)
plot(teta2,f32_A_norma)
xlabel('Angulo \theta_{2} [rad]','FontSize',11)
ylabel('FA [N]','FontSize',11)
grid on
hold off
%% graficos de las fuerzas en las juntas PARA B
hold on 
title('Fuerza en la junta B','FontSize',11)
subplot(311)
plot(teta2,x_fuerzas(5,:))
xlabel('Angulo \theta_{2} [rad]','FontSize',11)
ylabel('FB_{x} [N]','FontSize',11)
subplot(312)
plot(teta2,x_fuerzas(6,:))
xlabel('Angulo \theta_{2} [rad]','FontSize',11)
ylabel('FB_{y} [N]','FontSize',11)
subplot(313)
plot(teta2,f43_A_norma)
xlabel('Angulo \theta_{2} [rad]','FontSize',11)
ylabel('FB [N]','FontSize',11)
grid on
hold off
%% graficos de las fuerzas en las juntas PARA 02
hold on 
title('Fuerza en la junta 02','FontSize',11)
subplot(311)
plot(teta2,x_fuerzas(1,:))
xlabel('Angulo \theta_{2} [rad]','FontSize',11)
ylabel('F02_{x} [N]','FontSize',11)
subplot(312)
plot(teta2,x_fuerzas(2,:))
xlabel('Angulo \theta_{2} [rad]','FontSize',11)
ylabel('F02_{y} [N]','FontSize',11)
subplot(313)
plot(teta2,f_O2_norma)
xlabel('Angulo \theta_{2} [rad]','FontSize',11)
ylabel('F02 [N]','FontSize',11)
grid on
hold off
%% graficos de las fuerzas en las juntas PARA 04
hold on 
title('Fuerza en la junta 04','FontSize',11)
subplot(311)
plot(teta2,x_fuerzas(7,:))
xlabel('Angulo \theta_{2} [rad]','FontSize',11)
ylabel('F04_{x} [N]','FontSize',11)
subplot(312)
plot(teta2,x_fuerzas(8,:))
xlabel('Angulo \theta_{2} [rad]','FontSize',11)
ylabel('F04_{y} [N]','FontSize',11)
subplot(313)
plot(teta2,f_O4_norma)
xlabel('Angulo \theta_{2} [rad]','FontSize',11)
ylabel('F04 [N]','FontSize',11)
grid on
hold off
%% Grafico de torque requerido en la manivela de entrada 
hold on 
plot(teta2,x_fuerzas(9,:))
title('Torque requerido en la manivela de entrada','FontSize',11)
xlabel('Angulo \theta_{2} [rad]','FontSize',11)
ylabel('Torque [N.m]','FontSize',11)
legend('Torque requerido','FontSize',11)
grid on
hold off
%% Potencia requerida en la junta O2
for i=1:n_steps
   Potencia_2(i)=x_fuerzas(9,i)*omega2; %potencia en W 
end
%% graficos de la potencia requerida en el eslabon para una vuelta de la manevla
hold on 
plot(teta2,Potencia_2)
title('Potencia requerida en el eslabón 2 en funcion de \theta_{2}','FontSize',11)
xlabel('Angulo \theta_{2} [rad]','FontSize',11)
ylabel('Potencia [W]','FontSize',11)
legend('Pontencia','FontSize',11)
grid on
hold off
%% Fuerza que el eslabon 3 ejerce sobre el eslabón 4 
hold on 
plot(teta2,Potencia_2)
title('Potencia requerida en el eslabón 2 en funcion de \theta_{2}','FontSize',11)
xlabel('Angulo \theta_{2} [rad]','FontSize',11)
ylabel('Potencia [W]','FontSize',11)
legend('Pontencia_{x}','FontSize',11)
grid on
hold off
%% fuerza de contacto normal a la superficie de la viga en B
hold on
%title('fuerza de contacto normal a la superficie de la viga en B','Fontsize',11)
plot(teta2,f43_A_norma)
xlabel('Angulo \theta_{2} [rad]','FontSize',11)
ylabel('F [N]','FontSize',11)
hold off
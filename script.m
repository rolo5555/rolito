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
%% Angulos de extremos de 
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

%%
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
    aux_angle(i)=x(i,2)-deg2rad(14)
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
 %angulo de 20   
    punto_A= a*[cos(teta2(28)) sin(teta2(28))];
    crank = line([O2(1),punto_A(1)],[O2(2),punto_A(2)]);grid on;
    punto_B= punto_A+b*[cos(x(28,1)) sin(x(28,1))]
    biela= line([punto_A(1) punto_B(1)],[punto_A(2) punto_B(2)])
    balancin=line([punto_B(1) O4(1)],[punto_B(2) O4(2)])
    bancada= line([O4(1) O2(1)],[O4(2) O2(2)])
%%
 %angulo de 120   
    punto_A= a*[cos(teta2(167)) sin(teta2(167))];
    crank = line([O2(1),punto_A(1)],[O2(2),punto_A(2)]);grid on;
    punto_B= punto_A+b*[cos(x(167,1)) sin(x(167,1))]
    biela= line([punto_A(1) punto_B(1)],[punto_A(2) punto_B(2)])
    balancin=line([punto_B(1) O4(1)],[punto_B(2) O4(2)])
    bancada= line([O4(1) O2(1)],[O4(2) O2(2)])
 
    
    
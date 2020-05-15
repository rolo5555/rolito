%% eje de entrada 
% Ax=b asi es el sistema a resolver  x= rax ray rbx rby
wt1=3528.8 % [N]
wr1=1645.5 % [N]
fp= 2669.13; %[N]
n_steps= 500;
teta= linspace(0,2*pi,n_steps);
for i=1:n_steps
b=[-wt1+fp*cos(teta(i));-wr1-fp*sin(teta(i));-wr1*(248.7/1000);-wt1*(248.7/1000)];
A=[1 0 1 0;0 1 0 1;0 (95/1000) 0 (365/1000);(95/1000) 0 (365/1000) 0];
end
%% calculo los modulos 
for i=1:n_steps
ra(i)= sqrt(xxx(1,i)^2+xxx(2,i)^2);
rb(i)= sqrt(xxx(3,i)^2+xxx(4,i)^2);
end
%% Ploteo de las reacciones en funcion del angulo de paso   
figure(1)
hold on
plot(teta,ra)
plot(teta,rb)
title('Modulo de las reacciones en funcion del angulo de la polea','FontSize',14)
xlabel('Anuglo de la polea \theta [rad]','FontSize',14)
ylabel('Fuerza [N]','FontSize',14)
legend('Ra','Rb','FontSize',14)
grid on
hold off
%% cuando son iguales
k=1;
for i=1:n_steps
    if abs(ra(1,i)-rb(1,i)) < 0.1
    algo(k)=i
    k=k+1;
    end
end
%% Calculo del eje de salida
%defino mi sistema 
f21x=-1900; %[N] %es la fuerza que ejerce cuando el 2 al 1 cuando el torque es max
f21y=12000; %[N]
wt2=11216.4; %[N]
wr2=5230.3; %[N]
b2=[(wt2-f21x);(wr2-f21y);(-wr2*0.2572-f21y*0.5*0.3);(wt2*(257.2/1000)-f21x*0.5*0.6)];
A2=[1 0 1 0;0 1 0 1; 0 (-85/1000) 0 (-515/1000); (85/1000) 0 (515/1000) 0]; 
x=inv(A2)*b2
rc=sqrt(x(1,1)^2+x(2,1)^2)
rd=sqrt(x(3,1)^2+x(4,1)^2)



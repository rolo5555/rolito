function A =aceleration(alfa,alfa2,omega,teta,teta2,omega2)
r=[2000 350 2000 1280]; 
A(1)= -r(2)*alfa2*sin(teta2)-r(2)*(omega2^2)*cos(teta2)-r(3)*alfa(1)*sin(teta(1))-r(3)*(omega(1)^2)*cos(teta(1))+r(4)*alfa(2)*sin(teta(2))+r(4)*(omega(2)^2)*cos(teta(2));
A(2)= r(1)*alfa(1)*cos(teta(1))-r(4)*alfa(2)*cos(teta(2))-r(2)*(omega2^2)*sin(teta2)-r(3)*(omega(1)^2)*sin(teta(1))+r(4)*(omega(2)^2)*sin(teta(2));
end 

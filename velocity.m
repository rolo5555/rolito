function V =velocity(omega,teta,teta2,omega2)
r=[2000 350 2000 1280]; 
V(1)= -r(2)*omega2*sin(teta2)-r(3)*omega(1)*sin(teta(1))+r(4)*omega(2)*sin(teta(2));
V(2)= r(2)*omega2*cos(teta2)+r(3)*omega(1)*cos(teta(1))-r(4)*omega(2)*cos(teta(2));
end 

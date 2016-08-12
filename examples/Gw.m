function out= Gw(I0,N,m,ASPprms)

b=ASPprms(1);
psi=ASPprms(2);
alpha=ASPprms(3);
t=linspace(-pi/15,pi/15,N);
p=linspace(-pi/15,pi/15,N);
[T,P]=meshgrid(t,p);
Trot=T*cosd(psi)+P*sind(psi);

out= I0.*(1 + m*cos(b*Trot+alpha*pi/180));





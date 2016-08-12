function out=Amp(N,sc)

sigx=sc*pi/15;
sigy=sc*pi/15;

t=linspace(-pi/15,pi/15,N);
p=linspace(-pi/15,pi/15,N);
[T,P]=meshgrid(t,p);


 out=1/(2*pi*sigx*sigy)*exp(-0.5*((T/sigx).^2+(P/sigy).^2));

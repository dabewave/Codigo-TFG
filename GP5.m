% SOLVES THE 1D GPE VIA THE SPLIT-STEP FOURIER METHOD
clear all;clf; %Clear workspace and figure

hbar=1.054e-34; %Physical constants

c=1;
t=0; c1=-c; c2=c; x0=5;
M=1500; Nx=2*M+1;
dx=double(1e-2); x=(-M:1:M)*dx; %Define spatial grid
dk=pi/(M*dx); k=(-M:1:M)*dk; %Define k-space grid; co
dt=double(1e-3); Nt=2*x0/(dt*c2); %Define time step and number
V=-1; %Define potential

Nframe=100; %Data saved every Nframe steps
spacetime1=[]; spacetime2=[]; psi1_col=[]; psi2_col=[]; %Initialization

%Initial wavefunction
phi1=0; phi2=pi/2;
mu1_1=(3*c1^2+2*sqrt(3*c1^2+4)-4)/sqrt((3*c1^2+4)^(3/2)+18*c1^2-8);
mu2_1=3*c1*sqrt(4-c1^2)/sqrt((3*c1^2+4)^(3/2)+18*c1^2-8);
mu_1=3*(c1^2-4)*(sqrt(3*c1^2+4)-2)/((sqrt(3*c1^2+4)+4)*(3*c1^2+2*sqrt(3*c1^2+4)-4));
psi_1 = ((1i*mu1_1+mu2_1.*tanh(0.5.*sqrt(4-c1^2).*(x-x0-c1*t)))./sqrt(mu_1*tanh(0.5.*sqrt(4-c1^2).*(x-x0-c1*t)).^2+1))*exp(1i*phi1);

mu1_2=(3*c2^2+2*sqrt(3*c2^2+4)-4)/sqrt((3*c2^2+4)^(3/2)+18*c2^2-8);
mu2_2=3*c2*sqrt(4-c2^2)/sqrt((3*c2^2+4)^(3/2)+18*c2^2-8);
mu_2=3*(c2^2-4)*(sqrt(3*c2^2+4)-2)/((sqrt(3*c2^2+4)+4)*(3*c2^2+2*sqrt(3*c2^2+4)-4));
psi_2 = ((1i*mu1_2+mu2_2.*tanh(0.5.*sqrt(4-c2^2).*(x+x0-c2*t)))./sqrt(mu_2*tanh(0.5.*sqrt(4-c2^2).*(x+x0-c2*t)).^2+1))*exp(1i*phi2);

psi_0=abs(psi_1).^2+abs(psi_2).^2+2.*abs(psi_1).*abs(psi_2).*cos(phi2-phi1);

for itime=1:Nt %Time-stepping with split-step Fourier method
psi1 = ((1i*mu1_1+mu2_1.*tanh(0.5.*sqrt(4-c1^2).*(x-x0-c1*t)))./sqrt(mu_1*tanh(0.5.*sqrt(4-c1^2).*(x-x0-c1*t)).^2+1))*exp(1i*phi1);
psi1=psi1.*exp(-0.5*1i*dt*(V-abs(psi1).^4));
psi1_k=fftshift(fft(psi1)/Nx);
psi1_k=psi1_k.*exp(-0.5*dt*1i*k.^2);
psi1=ifft(ifftshift(psi1_k))*Nx;
psi1=psi1.*exp(-0.5*1i*dt*(V-abs(psi1).^4));
if mod(itime,Nt/Nframe) == 0 %Save wavefunction every Nframe steps
    spacetime1=vertcat(spacetime1,psi1); t
end
if (mod(itime,Nt/2) == 0) && (itime < Nt) 
    psi1_col=psi1;
end
t=t+dt;
end
t=0;
for itime=1:Nt %Time-stepping with split-step Fourier method
psi2=((1i*mu1_2+mu2_2.*tanh(0.5.*sqrt(4-c2^2).*(x+x0-c2*t)))./sqrt(mu_2*tanh(0.5.*sqrt(4-c2^2).*(x+x0-c2*t)).^2+1))*exp(1i*phi2);
psi2=psi2.*exp(-0.5*1i*dt*(V-abs(psi2).^4));
psi2_k=fftshift(fft(psi2)/Nx);
psi2_k=psi2_k.*exp(-0.5*dt*1i*k.^2);
psi2=ifft(ifftshift(psi2_k))*Nx;
psi2=psi2.*exp(-0.5*1i*dt*(V-abs(psi2).^4));
if mod(itime,Nt/Nframe) == 0 %Save wavefunction every Nframe steps
    spacetime2=vertcat(spacetime2,psi2); t
end
if (mod(itime,Nt/2) == 0) && (itime < Nt) 
    psi2_col=psi2;
end
t=t+dt;
end
% Gráficas
psi=abs(psi1).^2+abs(psi2).^2+2.*abs(psi1).*abs(psi2).*cos(phi2-phi1);
psi_col=abs(psi1_col).^2+abs(psi2_col).^2+2.*abs(psi1_col).*abs(psi2_col).*cos(phi2-phi1);
spacetime=abs(spacetime1).^2+abs(spacetime2).^2+2.*abs(spacetime1).*abs(spacetime2).*cos(phi2-phi1);
figure(1)
subplot(3,1,3)
plot(x,psi_0,'b',x,psi_col,'r',x,psi,'g--','LineWidth',1.5);
legend('\psi(x,0)','\psi(x,t_c)','\psi(x,t_f)');xlabel('x');ylabel('|\psi|^2');
ylim([0 max([max(psi_0),0.1])])
xlim([-10 10])
subplot(3,1,2)
dt_large=dt*double(Nt/Nframe);
surf(x,dt_large*(1:1:Nframe),spacetime); shading interp; colorbar;
xlabel('x'); ylabel('t'); zlabel('|\psi|^2')
xlim([-10 10])
clim([0 max([max(psi_0),0.1])])
subplot(3,1,1)
pcolor(x,dt_large*(1:1:Nframe),spacetime); shading interp; colorbar;
xlabel('x'); ylabel('t');
xlim([-10,10])
clim([0 max([max(psi_0),0.1])])

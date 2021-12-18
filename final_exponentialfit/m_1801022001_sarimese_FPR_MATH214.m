clear all;
close all;
format long;
Vs=2;
R=14.2;
L=0.98;
f=importdata('fprdata.dat');
x=f(:,1); % Voltage
y=f(:,2); % Current
%------------least square approximation (Exponential)-------------
xi2yi=0;yilnyi=0;xiyi=0;xiyilnyi=0;yi=0;
for i=1:1:length(f)
    xi2yi=xi2yi+((x(i)^2)*y(i));
    yilnyi=yilnyi+(y(i)*log(y(i)));  
    xiyi=xiyi+(x(i)*y(i));
    xiyilnyi=xiyilnyi+(x(i)*y(i)*log(y(i)));
    yi=yi+y(i);
end
a=exp(((xi2yi*yilnyi)-(xiyi*xiyilnyi))/((yi*xi2yi)-(xiyi^2)));
b=((yi*xiyilnyi)-(xiyi*yilnyi))/((yi*xi2yi)-(xiyi^2));
%------------plot-fit---------------
tx=0:0.025:1.3;
p=a*exp(b*tx);
figure
grid on;
hold on;
plot(x,y,'*');
plot(tx,p);
xlabel('Voltage (V)');
ylabel('Current (A)');
g1=legend('given values','exponential fit');
title(g1,'Diode Voltage-Current Exponential Fit');
%-----------runge kutta-Order4--timestep=0.025---------
T1=0:0.025:0.6;
idt=zeros(0,length(T1));
idt(1)=0;
Vd(1)=0;
for i=2:1:length(T1)
    k1=(0.025)*((1/L)*(Vs-Vd(i-1)-(idt(i-1)*R)));
    k2=(0.025)*((1/L)*(Vs-Vd(i-1)-((idt(i-1)+k1/2)*R)));
    k3=(0.025)*((1/L)*(Vs-Vd(i-1)-((idt(i-1)+k2/2)*R)));
    k4=(0.025)*((1/L)*(Vs-Vd(i-1)-((idt(i-1)+k3)*R)));
    idt(i)=idt(i-1)+(1/6)*(k1+2*k2+2*k3+k4);
    Vd(i)=(log(idt(i)/a))/b;
    VR(i)=R*idt(i);
    VL(i)=Vs-Vd(i)-VR(i);
end
%-----------plot-results--timestep=0.025------------
tx2=0:0.025:0.6;
figure
grid on;
hold on;
plot(tx2,idt,'r');
ylabel('Current(A)');
xlabel('Time (s)');
g2=legend('i(t) timestep size=25ms');
title(g2,'Current by time');
figure
grid on;
hold on;
plot(tx2,VR,'g');
plot(tx2,VL,'r');
plot(tx2,Vd,'b');
ylabel('Voltage(V)');
xlabel('Time (s)');
g3=legend('V1(t)','V2(t)','VD(t)');
title(g3,'Voltages by time(timestep size=25ms)');
%-------------runge kutta-Order4--timestep=0.0025------
T1=0:0.0025:0.6;
idt=zeros(0,length(T1));
idt(1)=0;
Vd(1)=0;
for i=2:1:length(T1)
    k1=(0.025)*((1/L)*(Vs-Vd(i-1)-(idt(i-1)*R)));
    k2=(0.025)*((1/L)*(Vs-Vd(i-1)-((idt(i-1)+k1/2)*R)));
    k3=(0.025)*((1/L)*(Vs-Vd(i-1)-((idt(i-1)+k2/2)*R)));
    k4=(0.025)*((1/L)*(Vs-Vd(i-1)-((idt(i-1)+k3)*R)));
    idt(i)=idt(i-1)+(1/6)*(k1+2*k2+2*k3+k4);
    Vd(i)=(log(idt(i)/a))/b;
    VR(i)=R*idt(i);
    VL(i)=Vs-Vd(i)-VR(i);
end
%-----------plot-results--timestep=0.0025-----------
tx2=0:0.0025:0.6;
figure
grid on;
hold on;
plot(tx2,idt,'r');
ylabel('Current(A)');
xlabel('Time (s)');
g2=legend('i(t) timestep size=2.5ms');
title(g2,'Current by time');
figure
grid on;
hold on;
plot(tx2,VR,'g');
plot(tx2,VL,'r');
plot(tx2,Vd,'b');
ylabel('Voltage(V)');
xlabel('Time (s)');
g3=legend('V1(t)','V2(t)','VD(t)');
title(g3,'Voltages by time(timestep size=2.5ms)');
%proof ODE numerical solution - given data
%test=0.057;
%test1=(log(test/a))/b
%test2=Vs-test*R-test*L
%error of fitted function
%for i=1:1:5
   % display(abs(y(i)-a*exp(b*x(i))));
%end
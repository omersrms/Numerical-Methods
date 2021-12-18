clear all;
close all;
format long;
V=2;
R=14.2;
L=0.98;
%------------analytical solution---------------
T1=0:0.05:0.6;
T2=0:0.025:0.6;
sonuc1=zeros(0,length(T1)+1);
sonuc2=zeros(0,length(T2)+1);
sonuc1(1)=0.1;
sonuc2(1)=0.1;
for i=2:1:length(T1)
    sonuc1(i)=(V/R)*(1-(exp(-1*(R/L)*T1(i))));
end
for i=2:1:length(T2)
    sonuc2(i)=(V/R)*(1-(exp(-1*(R/L)*T2(i))));
end
figure
grid on;
hold on;
plot(T1,sonuc1);
plot(T2,sonuc2);
ylabel('Current (Ampere)');
xlabel('time (second)');
a=legend('Current Value(h=0.05 second)','Current Value(h=0.025 second)');
title(a,'Current by Time on RL circuit');
%--------------Euler-h=0.05--------------------------
Euler1=zeros(0,length(T1));
Euler1(1)=0.1;
for i=2:1:length(T1)
    Euler1(i)=Euler1(i-1)+(0.05)*(((-14.4898)*Euler1(i-1))+12.2449);
end
%--------------Euler-h=0.025--------------------------
Euler2=zeros(0,length(T2));
Euler2(1)=0.1;
for i=2:1:length(T2)
    Euler2(i)=Euler2(i-1)+(0.025)*(((-14.4898)*Euler2(i-1))+12.2449);
end
%-----------Modified-Euler-h=0.05---------------------
M_Euler1=zeros(0,length(T1));
M_Euler1(1)=0.1;
for i=2:1:length(T1)
    M_Euler1(i)=M_Euler1(i-1)+(0.05/2)*((((-14.4898)*M_Euler1(i-1))+12.2449)+(((-14.4898)*(M_Euler1(i-1)+((0.05)*(((-14.4898)*M_Euler1(i-1))+12.2449))))+12.2449));
end
%-----------Modified-Euler-h=0.025---------------------
M_Euler2=zeros(0,length(T2));
M_Euler2(1)=0.1;
for i=2:1:length(T2)
    M_Euler2(i)=M_Euler2(i-1)+(0.05/2)*((((-14.4898)*M_Euler2(i-1))+12.2449)+(((-14.4898)*(M_Euler2(i-1)+((0.05)*(((-14.4898)*M_Euler2(i-1))+12.2449))))+12.2449));
end
%----------Midpoint-h=0.05-----------------------------
Midpoint1=zeros(0,length(T1));
Midpoint1(1)=0.1;
for i=2:1:length(T1)
    Midpoint1(i)=Midpoint1(i-1)+(0.05)*(((-14.4898)*(Midpoint1(i-1)+((0.05/2)*(((-14.4898)*Midpoint1(i-1))+12.2449))))+12.2449);
end
%----------Midpoint-h=0.025-----------------------------
Midpoint2=zeros(0,length(T2));
Midpoint2(1)=0.1;
for i=2:1:length(T2)
    Midpoint2(i)=Midpoint2(i-1)+(0.05)*(((-14.4898)*(Midpoint2(i-1)+((0.05/2)*(((-14.4898)*Midpoint2(i-1))+12.2449))))+12.2449);
end
%--------Runge-Kutte-Order4-h=0.05----------------------
RKO1=zeros(0,length(T1));
RKO1(1)=0.1;
for i=2:1:length(T1)
    k1=(0.05)*(((-14.4898)*RKO1(i-1))+12.2449);
    
    k2=(0.05)*(((-14.4898)*(RKO1(i-1)+(k1/2)))+12.2449);
    
    k3=(0.05)*(((-14.4898)*(RKO1(i-1)+(k2/2)))+12.2449);
    
    k4=(0.05)*(((-14.4898)*(RKO1(i-1)+(k3)))+12.2449);
    RKO1(i)=RKO1(i-1)+(1/6)*(k1+2*k2+2*k3+k4);
end
%--------Runge-Kutte-Order4-h=0.025----------------------
RKO2=zeros(0,length(T2));
RKO2(1)=0.1;
for i=2:1:length(T2)
    k1=(0.05)*(((-14.4898)*RKO2(i-1))+12.2449);
    
    k2=(0.05)*(((-14.4898)*(RKO2(i-1)+(k1/2)))+12.2449);
    
    k3=(0.05)*(((-14.4898)*(RKO2(i-1)+(k2/2)))+12.2449);
    
    k4=(0.05)*(((-14.4898)*(RKO2(i-1)+(k3)))+12.2449);
    RKO2(i)=RKO2(i-1)+(1/6)*(k1+2*k2+2*k3+k4);
end
%---------------methods-graphs--------------------
figure
hold on;
grid on;
plot(T1,sonuc1,'--');
plot(T1,Euler1);
plot(T1,M_Euler1);
plot(T1,Midpoint1);
plot(T1,RKO1);
ylabel('Current (Ampere)');
xlabel('time (second)');
a=legend('Exact Value','Euler','Modified Euler','Midpoint','Runge-Kutta order4');
title(a,'METHODS on RL circuit(h=0.05 second)');

figure
hold on;
grid on;
plot(T2,sonuc2,'--');
plot(T2,Euler2);
plot(T2,M_Euler2);
plot(T2,Midpoint2);
plot(T2,RKO2);
ylabel('Current (Ampere)');
xlabel('time (second)');
a=legend('Exact Value','Euler','Modified Euler','Midpoint','Runge-Kutta order4');
title(a,'METHODS on RL circuit(h=0.025 second)');
%---------bound error-------------------------------------------------
B_Error=zeros(0,length(T2));
for i=1:1:length(T2)
    h2=0.025;
    Lipt=14.49;
    M=(-156.361)*exp(-8.694);
    B_Error(i)=((h2*M)/(2*Lipt))*((exp(T2(i)*Lipt))-1);
end
figure
hold on;
grid on;
display(B_Error);
plot(T2,B_Error);
ylabel('Error Bound');
xlabel('time (second)');
a=legend('Error Bound');
title(a,'Error Bound of Methods on RL circuit(h=0.025 second)');

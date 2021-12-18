clear all;
close all;
format long;
L=0.1;
f=importdata('pr3data.dat');
%-------current-voltage graph-------------------------
figure 
grid on;
hold on;
plot(f(:,1),f(:,2));
plot(f(:,1),f(:,3));
ylabel('Value');
xlabel('time(second)');
a=legend('Current Value(Ampere)','Voltage Value(Volt)');
title(a,'Current-Voltage by Time');
%-------Power graph---------------------------------
power=zeros(1,length(f));
total=0;
for i=1:1:length(f)
    power(i)=f(i,2)*f(i,3);
    total=power(i)+power(i);
end
display("power of inductance:");
display(total);
figure 
grid on;
hold on;
plot(f(:,1),power(:));
ylabel('Power Value');
xlabel('time(second)');
a=legend('Power Value(Watt)');
title(a,'Power by Time');
%----------Stored Energy by data----------------------
Energy_D=zeros(1,length(f));
for i=1:1:length(f)
    Energy_D(i)=(1/2)*L*((f(i,2))^2);
end
figure 
grid on;
hold on;
plot(f(:,1),Energy_D(:));
ylabel('Energy(Joule)');
xlabel('time(second)');
a=legend('Energy(Joule)');
title(a,'Energy by Time');
EXACT_VALUE=(0.5)*(0.1)*((f(length(f),2)-f(1,2))^2); % stored energy
%----------composite trapezoidal----------------------
s3=0;
n=length(f);
h=((f(length(f),1)-f(1,1))/n);
for i=1:1:n-1
    s3=s3+(f(i,2)*f(i,3));
end
result_T=(h/2)*((f(1,2)*f(1,3))+(f(length(f),2)*f(length(f),3))+(2*s3));
display(result_T);
%----------composite midpoint--------------------------
s4=(f(1,2)*f(1,3)); % i=0 için sonuç
n=length(f)-1;
h=((f(length(f),1)-f(1,1))/n+2);
for i=1:1:n/2
    s4=s4+(f((2*i),2)*f((2*i),3));
end
result_M=2*h*s4;
display(result_M);
%--------composite simpson----------------------------
s1=0;
s2=0;
n=length(f)+1;
h=((f(length(f),1)-f(1,1))/n);
for i=1:1:(n/2)-1 % 2-40 arasi çiftler
    s1=s1+(f(2*i,2)*f(2*i,3));
end
for i=1:1:(n/2) % 1-41 arasi tekler
    s2=s2+(f((2*i)-1,2)*f((2*i)-1,3));
end
result_S=(h/3)*((f(1,2)*f(1,3))+(f(length(f),2)*f(length(f),3))+(2*s1)+(4*s2));
display(result_S);
%-------------------------------------------------------

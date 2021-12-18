clear all;
close all;
format long;
f=importdata('pr5data.dat');
x=f(:,1)/1000;
y=f(:,2)*1000;
%----------linear-----------
m=length(f);
a0=0;
a1=0;
yi=y(1);
xi=x(1);
xi2=x(1)^2;
xiyi=x(1)*y(1);
for j=2:1:m
    xi2=(x(j)^2)+xi2;
end
for j=2:1:m
    xiyi=(x(j)*y(j))+xiyi;
end
for j=2:1:m
    yi=y(j)+yi;
end
for j=2:1:m
    xi=x(j)+xi;
end
a0=((xi2*yi)-(xiyi*xi))/((m*xi2)-xi^2);
a1=((m*xiyi)-(xi*yi))/((m*xi2)-xi^2);
yl=zeros(1,m);
for i=1:1:m
    yl(i)=(a1*x(i))+a0;
end
%--------degree n=3------
n=3;
A=zeros(n+1,n+1);

for a=1:1:n+1                 % For the First Equation 
    xi=x(1)^(a-1);            % because of matlab cannot start loop
    for t=2:1:m               % with initial value equals zero.
        xi=xi+(x(t)^(a-1));
    end
    A(1,a)=xi;
end

for j=1:1:n                   % For the (2...n+1)th Equations
    count=1;
    for i=j:1:j+n 
        xi=x(1)^i;           
        for k=2:1:m  
            xi=xi+(x(k)^(i));
        end
        A(j+1,count)=xi;
        count=count+1;
    end
end
B=zeros(n+1,1);
for c=1:1:n+1
    yixi=y(1)*(x(1)^(c-1));
    for t=2:1:m
        yixi=yixi+(y(t)*(x(t)^(c-1)));
    end
    B(c,1)=yixi;
end
inv_A=inv(A);
Coef=(inv_A)*B;
p=[Coef(4) Coef(3) Coef(2) Coef(1)];
%------------------------n=2---------------------------------
n=2;
A2=zeros(n+1,n+1);

for a=1:1:n+1                 % For the First Equation 
    xi=x(1)^(a-1);            % because of matlab cannot start loop
    for t=2:1:m               % with initial value equals zero.
        xi=xi+(x(t)^(a-1));
    end
    A2(1,a)=xi;
end

for j=1:1:n                   % For the (2...n+1)th Equations
    count=1;
    for i=j:1:j+n 
        xi=x(1)^i;           
        for k=2:1:m  
            xi=xi+(x(k)^(i));
        end
        A2(j+1,count)=xi;
        count=count+1;
    end
end
B2=zeros(n+1,1);
for c=1:1:n+1
    yixi=y(1)*(x(1)^(c-1));
    for t=2:1:m
        yixi=yixi+(y(t)*(x(t)^(c-1)));
    end
    B2(c,1)=yixi;
end
inv_A2=inv(A2);
Coef2=(inv_A2)*B2;
p2=[Coef2(3) Coef2(2) Coef2(1)];
%------------plot-----------------
figure
grid on;
hold on;
plot(x,y,'.');
plot(x,yl,'g');
plot(x,polyval(p,x));
plot(x,polyval(p2,x));
xlabel('Distance (Km)');
ylabel('Signal Level (mV)');
a=legend('Exact Values','Linear Polynomial','Second Degree Polynomial','Third Degree Polynomial');
title(a,'Comparison');
%------------Error------------
E_linear=0;
for i=1:1:m
    E_linear=E_linear+(y(i)-yl(i))^2;
end
display(E_linear);

E_poly2=0;
for i=1:1:m
    pxi=Coef2(1)+Coef2(2)*x(i)+Coef2(3)*x(i)^2;
    E_poly2=E_poly2+((y(i)-pxi)^2);
end
display(E_poly2);

E_poly3=0;
for i=1:1:m
    pxi=Coef(1)+Coef(2)*x(i)+Coef(3)*x(i)^2+Coef(4)*x(i)^3;
    E_poly3=E_poly3+((y(i)-pxi)^2);
end
display(E_poly3);
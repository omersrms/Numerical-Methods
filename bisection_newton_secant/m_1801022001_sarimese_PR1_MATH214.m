clear all;
close all;

format long; % reason why is be able to display more decimal digits 
tol=10^(-10); % tolerance value

%-----------------------Bisection---------------------

UB=10; % Upper Bound Value
LB=-3; % Lower Bound Value
n=1; % iteration variable
temp=1; % temporary value
a_1=zeros(20,5); % space for saving values

while(temp>tol && temp~=0)
    % saving LB,UB values before creating table...
    a_1(n,1)=n;
    a_1(n,2)=LB;
    a_1(n,3)=UB;
 
    mid=LB+(UB-LB)/2;
    temp=abs(UB-LB)/(2^n);
    if(f(UB)*f(mid)>0)
        UB=mid;
    else
        LB=mid;
    n=n+1;
    % saving res,temp values before creating table...    
    a_1(n-1,4)=mid;
    a_1(n-1,5)=temp;
    end  
end

disp("__Bisection Method__");
% creating & arranging & displaying table 
T=array2table(a_1);
T.Properties.VariableNames={'iteration' 'Lowerbound' 'Upperbound' 'midpoint' 'error'}; 
disp(T);

%---------------------------newton--------------------------

prev=-3; 
next=1;  
i=1; % iteration variable
temp=1; % temporary value
a_2=zeros(20,3); % space for saving values
a_2(1,2)=prev;
a_2(1,3)=abs(next-prev);

while(tol<temp)
    a_2(i,1)=i;
    next=(prev)-(f(prev)/fd(prev));
    temp=abs((next)-(prev));
    prev=next;
    i=i+1;
    a_2(i,2)=next;
    a_2(i,3)=temp;
end
a_2(i,1)=i;

disp("__Newton Method__");
% creating & arranging & displaying table
T=array2table(a_2);
T.Properties.VariableNames={'iteration' 'root' 'error'};
disp(T);
%---------------------------secant--------------------------

prev=-3;
mid=10;
next=0;
temp=1; % temporary value
i=1; % iteration variable
a_3=zeros(20,3); % space for saving values

while(temp>tol)
    next=mid-((f(mid)*(prev-mid))/(f(prev)-f(mid)));
    temp=abs(next-prev);
    prev=mid;
    mid=next;
    a_3(i,1)=i;
    a_3(i,2)=next;
    a_3(i,3)=temp;
    i=i+1;
end

disp("__Secant Method__");
% creating & arranging & displaying table
T=array2table(a_3);
T.Properties.VariableNames={'iteration' 'root' 'error'}; 
disp(T);

%--------------------------graphs---------------------------
% i didnt find any other way to plot double array values.
% i could not use loglog() & semilog().
hold on;
grid on;
p=plot(a_1(:,5));
p=plot(a_2(:,3));
p=plot(a_3(:,3));
legend('Bisection method','Newton method','Secant method');
xlabel('iteration number');
ylabel('error value');
axis([0 10 -1 10]);

%-------------------------functions-------------------------
function val=f(x)
val=(1/(4*pi*(1/36*pi)*10^(-9)))*((13*(x+7)/(abs((x+7)^3)))+(9*(x+4)/(abs((x+4)^3)))+(5*(x-11)/(abs((x-11)^3)))+(3*(x-15)/(abs((x-15)^3))));
end

function deg=fd(x)
deg=(1/(4*pi*(1/36*pi)*10^(-9)))*((3/((abs(x-15))^3))+(5/((abs(x-11))^3))+(9/((abs(x+4))^3))+(13/((abs(x+7))^3))-((9*(x-15)^2)/((abs(x-15))^5))-((15*(x-11)^2)/((abs(x-11))^5))-((27*(x+4)^2)/((abs(x+4))^5))-((39*(x+7)^2)/((abs(x+7))^5)));
end
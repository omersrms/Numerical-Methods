clear all;
close all;
format long;

C1=importdata('current1.dat');
C2=importdata('current2.dat');
C3=importdata('current3.dat');
C4=importdata('current4.dat');
L=0.98;
R=14.2;

%---------------Derivatives of current--------------------

%---------------Forward difference-C1---------------------
D_F1=zeros(1,length(C1)-1);
D_F1(length(C1))=nan;
for i=1:1:length(C1)-1 
   D_F1(i)=((C1(i+1,2))-(C1(i,2)))/((C1(i+1,1))-(C1(i,1))); 
end
%---------------Forward difference-C2---------------------
D_F2=zeros(1,length(C2)-1);
D_F2(length(C2))=nan;
for i=1:1:length(C2)-1
    D_F2(i)=((C2(i+1,2))-(C2(i,2)))/((C2(i+1,1))-(C2(i,1))); 
end
%---------------Forward difference-C3---------------------
D_F3=zeros(1,length(C3)-1);
D_F3(length(C3))=nan;
for i=1:1:length(C3)-1
    D_F3(i)=((C3(i+1,2))-(C3(i,2)))/((C3(i+1,1))-(C3(i,1))); 
end
%---------------Forward difference-C4---------------------
D_F4=zeros(1,length(C4)-1);
D_F4(length(C4))=nan;
for i=1:1:length(C4)-1
    D_F4(i)=((C4(i+1,2))-(C4(i,2)))/((C4(i+1,1))-(C4(i,1))); 
end

%--------------Backward difference-C1---------------------
D_B1=zeros(1,length(C1));
D_B1(1)=nan;
for i=1:1:length(C1)-1
    D_B1(i+1)=((C1(i+1,2))-(C1(i,2)))/((C1(i,1))-(C1(i+1,1)));
end
%--------------Backward difference-C2---------------------
D_B2=zeros(1,length(C2));
D_B2(1)=nan;
for i=1:1:length(C2)-1
    D_B2(i+1)=((C2(i+1,2))-(C2(i,2)))/((C2(i,1))-(C2(i+1,1)));
end
%--------------Backward difference-C3---------------------
D_B3=zeros(1,length(C3));
D_B3(1)=nan;
for i=1:1:length(C3)-1
    D_B3(i+1)=((C3(i+1,2))-(C3(i,2)))/((C3(i,1))-(C3(i+1,1)));
end
%--------------Backward difference-C4---------------------
D_B4=zeros(1,length(C4));
D_B4(1)=nan;
for i=1:1:length(C4)-1
    D_B4(i+1)=((C4(i+1,2))-(C4(i,2)))/((C4(i,1))-(C4(i+1,1)));
end
%-------------Three point-C1-------------------------------
D_T1=zeros(1,length(C1));
for i=1:1:length(C1)
h=C1(2,1)-C1(1,1);
    if i==1
       D_T1(i)=(1/(2*h))*(((-3)*C1(i,2))+(4*C1(i+1,2))-((-1)*C1(i+2,2)));
    elseif i==length(C1)
       D_T1(i)=(1/(-2*h))*(((-3)*C1(i,2))+(4*C1(i-1,2))-((-1)*C1(i-2,2)));
    else
        D_T1(i)=(1/(2*h))*(C1(i+1,2)-C1(i-1,2));
    end
end
%-------------Three point-C2-------------------------------
D_T2=zeros(1,length(C2));
for i=1:1:length(C2)
h=C2(2,1)-C2(1,1);
    if i==1
       D_T2(i)=(1/(2*h))*(((-3)*C2(i,2))+(4*C2(i+1,2))-((-1)*C2(i+2,2)));
    elseif i==length(C2)
       D_T2(i)=(1/(-2*h))*(((-3)*C2(i,2))+(4*C2(i-1,2))-((-1)*C2(i-2,2)));
    else
        D_T2(i)=(1/(2*h))*(C2(i+1,2)-C2(i-1,2));
    end
end
%-------------Three point-C3-------------------------------
D_T3=zeros(1,length(C3));
for i=1:1:length(C3)
h=C3(2,1)-C3(1,1);
    if i==1
       D_T3(i)=(1/(2*h))*(((-3)*C3(i,2))+(4*C3(i+1,2))-((-1)*C3(i+2,2)));
    elseif i==length(C3)
       D_T3(i)=(1/(-2*h))*(((-3)*C3(i,2))+(4*C3(i-1,2))-((-1)*C3(i-2,2)));
    else
        D_T3(i)=(1/(2*h))*(C3(i+1,2)-C3(i-1,2));
    end
end
%-------------Three point-C4-------------------------------
D_T4=zeros(1,length(C4));
for i=1:1:length(C4)
h=C4(2,1)-C4(1,1);
    if i==1
       D_T4(i)=(1/(2*h))*(((-3)*C4(i,2))+(4*C4(i+1,2))-((-1)*C4(i+2,2)));
    elseif i==length(C4)
       D_T4(i)=(1/(-2*h))*(((-3)*C4(i,2))+(4*C4(i-1,2))-((-1)*C4(i-2,2)));
    else
        D_T4(i)=(1/(2*h))*(C4(i+1,2)-C4(i-1,2));
    end
end

%------------Voltages-------------------------------

%-----------Voltage-C1------------------------------
V_F1=zeros(1,length(C1));
for i=1:1:length(C1)
    V_F1(i)=(C1(i,2)*R)+(D_F1(i)*L);
end
V_B1=zeros(1,length(C1));
for i=1:1:length(C1)-1
    V_B1(i)=(C1(i,2)*R)+(D_B1(i)*L);
end
V_T1=zeros(1,length(C1));
for i=1:1:length(C1)
    V_T1(i)=(C1(i,2)*R)+(D_T1(i)*L);
end
%-----------Voltage-C2------------------------------
V_F2=zeros(1,length(C2));
for i=1:1:length(C2)
    V_F2(i)=(C2(i,2)*R)+(D_F2(i)*L);
end
V_B2=zeros(1,length(C2));
for i=1:1:length(C2)-1
    V_B2(i)=(C2(i,2)*R)+(D_B2(i)*L);
end
V_T2=zeros(1,length(C2));
for i=1:1:length(C2)
    V_T2(i)=(C2(i,2)*R)+(D_T2(i)*L);
end
%-----------Voltage-C3------------------------------
V_F3=zeros(1,length(C3));
for i=1:1:length(C3)
    V_F3(i)=(C3(i,2)*R)+(D_F3(i)*L);
end
V_B3=zeros(1,length(C3));
for i=1:1:length(C3)-1
    V_B3(i)=(C3(i,2)*R)+(D_B3(i)*L);
end
V_T3=zeros(1,length(C3));
for i=1:1:length(C3)
    V_T3(i)=(C3(i,2)*R)+(D_T3(i)*L);
end
%-----------Voltage-C4------------------------------
V_F4=zeros(1,length(C4));
for i=1:1:length(C4)
    V_F4(i)=(C4(i,2)*R)+(D_F4(i)*L);
end
V_B4=zeros(1,length(C4));
for i=1:1:length(C4)-1
    V_B4(i)=(C4(i,2)*R)+(D_B4(i)*L);
end
V_T4=zeros(1,length(C4));
for i=1:1:length(C4)
    V_T4(i)=(C4(i,2)*R)+(D_T4(i)*L);
end
%-----------------Graphs--------------------------
figure
grid on;
hold on;
plot(C1(:,1),D_F1);
plot(C1(:,1),D_B1);
plot(C1(:,1),D_T1);
ylabel('Derivative of current');
xlabel('time(second)');
xlim auto
ylim auto
a=legend('Forward diff.','Backward diff.','Three point diff.');
title(a,'Derivatives of Current');

figure
grid on;
hold on;
plot(C1(:,1),V_F1);
plot(C1(:,1),V_B1);
plot(C1(:,1),V_T1);
ylabel('Voltage (V)'); 
xlabel('time(second)'); 
xlim auto
ylim auto
a=legend('Forward diff.','Backward diff.','Three point diff.');
title(a,'impressed voltages');
%---------------Tables-------------------------------
%--------------Forward Difference Method-------------
T=zeros(length(C4),3);
for i=1:1:length(C4)
    T(i,1)=C4(i,1);    
    T(i,2)=D_F4(i);
    if i>1
        T(i,3)=D_F4(i-1)-D_F4(i);  
    end
end
Forward_method=array2table(T);
Forward_method.Properties.VariableNames={'Time' 'Derivative_of_current' 'error_convergence'};
display('Forward Difference Method Converge Table');
display(Forward_method);
%--------------Backward Difference Method-------------
T2=zeros(length(C4),3);
for i=1:1:length(C4)
    T2(i,1)=C4(i,1);    
    T2(i,2)=D_B4(i);
    if i>1
        T2(i,3)=D_B4(i-1)-D_B4(i);  
    end
end
Backward_method=array2table(T2);
Backward_method.Properties.VariableNames={'Time' 'Derivative_of_current' 'error_convergence'};
display('Backward Difference Method Converge Table');
display(Backward_method);
%-------------Three point mid/end point Method----------
T3=zeros(length(C4),3);
for i=1:1:length(C4)
    T3(i,1)=C4(i,1);    
    T3(i,2)=D_T4(i);
    if i>1
        T3(i,3)=D_T4(i-1)-D_T4(i);  
    end
end
Three_point_method=array2table(T3);
Three_point_method.Properties.VariableNames={'Time' 'Derivative_of_current' 'error_convergence'};
display('Three point mid/end point Method Converge Table');
display(Three_point_method);
%-------------------------------------------------------





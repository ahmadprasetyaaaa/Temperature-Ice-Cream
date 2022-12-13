clc
clear
load DataEskrimFix_2.mat

T19E_2 = D2_H1D9_T_eskrim;
t19_2 = D2_H1D9_t_menit;
T19L_2 = D2_H1D9_T_lingkungan; 

T57E_2 = D2_H5D7_T_eskrim;
t57_2 = D2_H5D7_t_menit;
T57L_2 = D2_H5D7_T_lingkungan; 

T115E_2 = D2_H11D5_T_eskrim;
t115_2 = D2_H11D5_t_menit;
T115L_2 = D2_H11D5_t_lingkungan; 
n2 = length(t115_2);
    
%% Pengolahan Data 
%1 Curve Fitting + diferentiation
% H1D9
x = t19_2;
y = T19E_2 - T19L_2;
l =  T19E_2;
n = n2;
disp('Data 2 : H1D9 ');
[a, r2] = linearRegression_khusus((x),log(y));
ynn_0 = @(tn1) -45.916870*exp(-0.023442*tn1) + mean(T19L_2);
yd = ynn_0(x);
%Diferentiation
dTdt1 = diff(yd)./diff(x);
n = length(x);
xm1 = (x(1:n-1)+x(2:n))./2;
% figure(3)
% plot(xm,dTdt )
% Bisection 
ynn_akhir = @(tn1) -34.128903*exp(-0.020657*tn1) + mean(T19L_2) - mean(T19L_2) + 0.1;
[root1,fx,ea,iter] = biSection(ynn_0,-5,20,0.001);
[root2,fx,ea,iter] = biSection(ynn_akhir,-5,1000000,0.001);
%Lagrange 
T0_1 = lagrange(l,x,0);

fprintf('T saat lapisan eskrim 1 cm cair (newton law cooling) :%f menit \n',root1);
fprintf('T saat lapisan eskrim 1 cm cair (Interpolasi lagrange) :%f menit \n',T0_1 );
fprintf('T saat lapisan eskrim 1 cm mencapai saturasi :%f menit \n',root2);
disp(' ');

% H11D5
x = t115_2;
y = T115E_2 - T115L_2;
l = T115E_2;
n = n2;
disp('Data 2 : H11D5 ');
[a, r2] = linearRegression_khusus((x),log(y));
% tn1 = linspace(0,2000,200);
% Tn1 = -43.075267*exp(-0.010126*tn1) + mean(T57L_2);
ynn_0 = @(tn1) -43.075267*exp(-0.010126*tn1) + mean(T57L_2);
yd = ynn_0(x);
%Diferentiation
dTdt3 = diff(yd)./diff(x);
n = length(x);
xm3 = (x(1:n-1)+x(2:n))./2;
%Lagrange 
T0_2 = lagrange(l,x,0);
%Bisection
ynn_akhir = @(tn1) -43.075267*exp(-0.010126*tn1) + mean(T57L_2) - mean(T57L_2) + 0.1;
[root1,fx,ea,iter] = biSection(ynn_0,-5,60,0.001);
[root2,fx,ea,iter] = biSection(ynn_akhir,-5,1000000,0.001);
fprintf('T saat lapisan eskrim 11 cm cair :%f menit \n',root1);
fprintf('T saat lapisan eskrim 11 cm cair (Interpolasi lagrange) :%f menit \n',T0_2 );
fprintf('T saat lapisan eskrim 11 cm mencapai saturasi :%f menit \n',root2);
disp(' ');
%----------------------------------------------------------------------------------------------
%2 Perbandingan Kurva estimasi Data 2
tn1_2 = linspace(0,2000,200);
tn2_2 = linspace(0,60,200);%Interval hasil curve fitting
Tn1_2 = -45.916870*exp(-0.023442*tn1_2) + mean(T19L_2);
Tn1_2_1 = -45.916870*exp(-0.023442*tn2_2) + mean(T19L_2);%Kurve fitting
Tn3_2 = -43.075267*exp(-0.010126*tn1_2) + mean(T115L_2);
Tn3_2_3 = -43.075267*exp(-0.010126*tn2_2) + mean(T115L_2);%Kurve fitting
figure(1);
title('Ploting Data Temperatur Es krim Variasi Kedalaman dan Diameter Cup');
plot(t19_2,T19E_2,'ro');
hold on
% plot(t19_2,T57E_2,'r*');
plot(t19_2,T115E_2,'kd');
legend('h = 1cm, d = 9 cm','h = 11cm , d = 9 cm')
grid on
title('Kurva temperatur mencapai posisi setimbang')
xlabel('Waktu(menit)')
ylabel('Temperatur es krim(C)')
hold off

figure(2);
plot(t19_2,T19E_2,'ro');
hold on
plot(t19_2,T115E_2,'kd');
plot(tn2_2,Tn1_2_1,'r-');
plot(tn2_2,Tn3_2_3,'k-');
legend('h = 1cm, d = 9 cm','h = 11cm , d = 9 cm')
grid on
title('Curve Fitting Temperatur Es krim Variasi Kedalaman dan Diameter Cup');
xlabel('Waktu(menit)')
ylabel('Temperatur es krim(C)')
hold off

figure(3)
plot(tn1_2,Tn1_2,'r-')
hold on
% plot(tn1_2,Tn2_2,'r-')
plot(tn1_2,Tn3_2,'k-')
legend('h = 1cm, d = 9 cm','h = 11cm , d = 9 cm')
grid on
title('Kurva temperatur mencapai posisi setimbang')
xlabel('Waktu(menit)')
ylabel('Temperatur es krim(C)')
hold off

figure(4)
plot(xm1,dTdt1,'ro')
hold on
plot(xm3,dTdt3,'ko')
title('Laju Kenaikan Temperatur es Krim')
xlabel('Waktu(menit)')
ylabel('DT/dt (C / menit)')
grid on
legend('h = 1cm, d = 9 cm','h = 11cm , d = 9 cm')
hold off

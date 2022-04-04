%%1-D Self Potential linear Inversion (SP forward modelling : El-Kaliboy dan Al-Gami (2009))
%Inclined sheet anomaly using Simulated Annealing 
% Algorithm
%Mohammad Rheza Zamani
clear all;
clc;
%Parameter
k  = 100; %Amplitudo polarisasi
z = 15; %Kedalaman dari permukaan ke titik tengah sheet
x0 = 5;  % Jarak horizontal dari sheet
alpha = 30; %Sudut inklinasi dari sheet
a = 10; %1/2 jarak lebar dari sheet

%Jarak pengukuran
x = -100:1:100;

%Data Sintetik
[V_obs] = fwd_SP(k,x0,z,alpha,a,x);
%Definisi Ruang Model
nitr = 300; 
k_min = 1;
k_max = 200;
x0_min = 1;
x0_max = 10;
z_min = 1;
z_max = 30;
alpha_min = 1;
alpha_max = 80;
a_min = 1;
a_max = 20;
T = 5;
dec = 0.05;

%Membuat Video
v = VideoWriter('Simulated Annealing Anomaly SP.avi');
open(v);
%Model acak
%Membuat model acak
rand1 = rand;
rand2 = rand;
rand3 = rand;
rand4 = rand;
rand5 = rand;
model1(1,1) = k_min + rand1*(k_max-k_min);
model1(1,2) = x0_min + rand2*(x0_max-x0_min);
model1(1,3) = z_min + rand3*(z_max - z_min);
model1(1,4) = alpha_min + rand4*(alpha_max - alpha_min);
model1(1,5) = a_min + rand5*(a_max - a_min);
%model1(1,1) = 75;
%model1(1,2) = 8;
%model1(1,3) = 12;
%model1(1,4) = 50;
%model1(1,5) = 15;
V_cal1(1,:) = fwd_SP(model1(1),model1(2),model1(3),model1(4),model1(5),x);
misfit = misfit_SP(V_obs(1,:),V_cal1(1,:));
E1 = misfit;

for itr = 1 : nitr
    model2(1,1) = k_min + rand*(k_max-k_min);
    model2(1,2) = x0_min + rand*(x0_max-x0_min);
    model2(1,3) = z_min + rand*(z_max - z_min);
    model2(1,4) = alpha_min + rand*(alpha_max - alpha_min);
    model2(1,5) = a_min + rand*(a_max - a_min);
    V_cal2(1,:) = fwd_SP(model2(1),model2(2),model2(3),model2(4),model2(5),x);
    E2 = misfit_SP(V_obs(1,:),V_cal2(1,:));
    delta_E = E2 - E1;
    if delta_E<0
        model1 = model2;
        E1 = E2;
        V_cal1 = V_cal2;
            if model1(1)<k_min 
                model1(1) = k_min;
            end
            if model1(2)<x0_min
                model1(2) = x0_min;
            end
            if model1(3)<z_min
                model1(3) = z_min;
            end
            if model1(4)<alpha_min
                model1(4) = alpha_min;
            end
            if model1(5)<a_min
                model1(5) = a_min;
            end
            if model1(1)>k_max 
                model1(1) = k_max;
            end
            if model1(2)>x0_max
                model1(2) = x0_max;
            end
            if model1(3)>z_max
                model1(3) = z_max;
            end
            if model1(4)>alpha_max
                model1(4) = alpha_max;
            end
            if model1(5)>a_max
                model1(5) = a_max;
            end
    else
        P = exp(-delta_E/T); 
        if P>= rand
           model1 = model2;
           E1 = E2;
           V_cal1 = V_cal2;
                       if model1(1)<k_min 
                model1(1) = k_min;
            end
            if model1(2)<x0_min
                model1(2) = x0_min;
            end
            if model1(3)<z_min
                model1(3) = z_min;
            end
            if model1(4)<alpha_min
                model1(4) = alpha_min;
            end
            if model1(5)<a_min
                model1(5) = a_min;
            end
            if model1(1)>k_max 
                model1(1) = k_max;
            end
            if model1(2)>x0_max
                model1(2) = x0_max;
            end
            if model1(3)>z_max
                model1(3) = z_max;
            end
            if model1(4)>alpha_max
                model1(4) = alpha_max;
            end
            if model1(5)>a_max
                model1(5) = a_max;
            end
    end
    end
   Egen(itr) = E1;
   T = T*(1-dec);
   Temperature(itr) = T;

figure(1)
plot(x,V_obs,'r.',x,V_cal1,'b-','MarkerSize',20,'MarkerFaceColor','r','LineWidth',2.5);
axis([-100 100 -150 50])
legend({'Forward Model','Inversion Model'},'Location','Southeast')
title(['\bf \fontsize{10}\fontname{Times}SP Anomaly  || Misfit : ', num2str(Egen(itr)),' || iteration : ', num2str(itr)],'FontWeight','bold');
subtitle(['k = ',num2str(model1(1)),' || x0 = ',num2str(model1(2)),' || z = ',num2str(model1(3)),' || alpha = ',num2str(model1(4)),' || a = ',num2str(model1(5))],'FontWeight','bold')
xlabel('Distance (m)','FontWeight','bold');
ylabel('SP Anomaly (mV)','FontWeight','bold');
grid on
set(gcf, 'Position', get(0, 'Screensize'));
frame = getframe(gcf);
writeVideo(v,frame);
end
close(v);
%Plot grafik misfit
figure(2)
plot(1:nitr,Egen,'r','Linewidth',1.5)
xlabel('Iteration Number','FontSize',10,'FontWeight','Bold');
ylabel('RSME','FontSize',10,'FontWeight','Bold');
title('\bf \fontsize{12} Grafik Misfit ');
grid on

%Plot Temperature
figure(3)
plot(1:nitr,Temperature,'b','Linewidth',1.5)
xlabel('Iteration Number','FontSize',10,'FontWeight','Bold');
ylabel('Temperature','FontSize',10,'FontWeight','Bold');
title('\bf \fontsize{12} Grafik Penurunan Temperature ');
grid on

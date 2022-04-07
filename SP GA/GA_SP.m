%%1-D Self Potential linear Inversion (SP forward modelling : El-Kaliboy dan Al-Gami (2009))
%Inclined sheet anomaly using Genetic Algorithm
% Algorithm
%Inspired by M. Randy Caesario H. GA MATLAB code (https://github.com/mrch-hub/geoph-inversion)
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
npop = 100;
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
prob = 0.9; %Probabilitas crossover

%Membuat Video
v = VideoWriter('GA Anomaly SP.avi');
open(v);

for ipop = 1 : npop
    model(ipop,1) = k_min + rand*(k_max-k_min);
    model(ipop,2) = x0_min + rand*(x0_max-x0_min);
    model(ipop,3) = z_min + rand*(z_max - z_min);
    model(ipop,4) = alpha_min + rand*(alpha_max - alpha_min);
    model(ipop,5) = a_min + rand*(a_max - a_min);
    V_cal1(ipop,:) = fwd_SP(model(ipop,1),model(ipop,2),model(ipop,3),model(ipop,4),model(ipop,5),x);
    misfit(ipop) = misfit_SP(V_obs(1,:),V_cal1(ipop,:))
end


for itr = 1 : nitr
    %fitness
    fitness = 1./misfit;
    %fitness ternomalisasi
    fitness_norm = fitness./sum(fitness);
    sc = 0;
    for i = 1 : npop
        sc = sc + fitness_norm(i);
        cumm(i) = sc;
    end
    indx = 1;
    %Roulette wheel pemilihan induk
    for i=1:npop/2
      R1 = rand;
     for j=1:npop
        if R1 < cumm(j)
            ipar1 = j;
            break
        end
    end
    
    R2 = rand;
    for j=1:npop
        if R2 < cumm(j)
            ipar2 = j;
            break
        end
    end
    
    %Crossover dan offspring
    R3 = rand;
    if R3 < prob
        i1 = rand;
        i2 = rand;
        i3 = rand;
        i4 = rand;
        i5 =  rand;
        model_new(indx,1) = i1*model(ipar1,1)+(1-i1)*model(ipar2,1);
        model_new(indx+1,1) = i1*model(ipar2,1)+(1-i1)*model(ipar1,1);
        model_new(indx,2) =  i2*model(ipar1,2)+(1-i2)*model(ipar2,2);
        model_new(indx+1,2) = i2*model(ipar2,2)+(1-i2)*model(ipar1,2);
        model_new(indx,3) = i3*model(ipar1,3)+(1-i3)*model(ipar2,3);
        model_new(indx+1,3) = i3*model(ipar2,3)+(1-i3)*model(ipar1,3);
        model_new(indx,4) = i4*model(ipar1,4)+(1-i4)*model(ipar2,4);
        model_new(indx+1,4) = i4*model(ipar2,4)+(1-i4)*model(ipar1,4);
        model_new(indx,5) = i5*model(ipar1,5)+(1-i5)*model(ipar2,5);
        model_new(indx+1,5) = i5*model(ipar2,5)+(1-i5)*model(ipar1,5);
    else
        model_new(indx,:) = model(ipar1,:);
        model_new(indx+1,:) = model(ipar2,:);
    end
    %Menambah index untuk setiap iterasi
     indx = indx + 2; 
    end
    model = model_new;
    for i=1:npop
        V_cal1(i,:) = fwd_SP(model(i,1),model(i,2),model(i,3),model(i,4),model(i,5),x);
        misfit(i) = misfit_SP(V_obs(1,:),V_cal1(i,:));
    end
    Egen(itr) = misfit(i);
    fitgen(itr) = fitness(i);

figure(1)
plot(x,V_obs,'r.',x,V_cal1,'b-','MarkerSize',20,'MarkerFaceColor','r','LineWidth',2.5);
axis([-100 100 -150 50])
legend({'Forward Model','Inversion Model'},'Location','Southeast')
title(['\bf \fontsize{10}\fontname{Times}SP Anomaly  || Misfit : ', num2str(Egen(itr)),' || iteration : ', num2str(itr)],'FontWeight','bold');
subtitle(['k = ',num2str(model(i,1)),' || x0 = ',num2str(model(i,2)),' || z = ',num2str(model(i,3)),' || alpha = ',num2str(model(i,4)),' || a = ',num2str(model(i,5))],'FontWeight','bold')
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

%Plot grafik fitness
figure(3)
plot(1:nitr,fitgen,'b','Linewidth',1.5)
xlabel('Iteration Number','FontSize',10,'FontWeight','Bold');
ylabel('Fitness','FontSize',10,'FontWeight','Bold');
title('\bf \fontsize{12} Grafik Fitness ');
grid on

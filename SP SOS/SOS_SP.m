%%1-D Self Potential linear Inversion (Rumus SP forward modelling : El-Kaliboy dan Al-Gami (2009))
%Inclined sheet anomaly using mSOS algorithm
%Mohammad Rheza Zamani
clear all;
clc;
%Parameter
k  = 100; %Amplitudo polarisasi
z = 15; %Kedalaman dari permukaan ke titik tengah sheet
x0 = 5;  % Jarak horizontal dari sheet
alpha = 40; %Sudut inklinasi dari sheet
a = 10; %1/2 jarak lebar dari sheet

%Jarak pengukuran
x = -100:1:100;

%Data Sintetik
[V_obs] = fwd_SP(k,x0,z,alpha,a,x);
%Definisi Ruang Model
npop = 100; 
nitr = 200; 
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

%Membuat Video
v = VideoWriter('Data sintetik Anomaly SP.avi');
open(v);
%Model acak
%Membuat model acak
for ipop = 1 : npop
    model(ipop,1) = k_min + rand*(k_max-k_min);
    model(ipop,2) = x0_min + rand*(x0_max-x0_min);
    model(ipop,3) = z_min + rand*(z_max-z_min);
    model(ipop,4) = alpha_min + rand*(alpha_max-alpha_min);
    model(ipop,5) = a_min + rand*(a_max-a_min);
    V_cal(ipop,:) = fwd_SP(model(ipop,1),model(ipop,2),model(ipop,3),model(ipop,4),model(ipop,5),x);
    misfit = misfit_SP(V_obs(1,:),V_cal(ipop,:));
    err(ipop) = misfit;
end

%Proses Inversi
for itr = 1 : nitr
    for i = 1 : npop
        idx = find(err ==min(err));
        model_best = model(idx(1),:);
        %Mutualisme
        j = randi(npop,1);
        k = randi(npop,1);
        if j==i || k==i
            j = randi(npop,1);
            k = randi(npop,1);
        end
        model_mut = [model(i,:);model(j,:)];
        mv_m =(model(i,:)+model(j,:))/2;
        bf = 1;
        for l = 1 : 2
            mod_mut(l,:) = model_mut(l,:) + rand*(model(k)-mv_m*bf);
            if mod_mut(l,1)<k_min 
                mod_mut(l,1) = k_min;
            end
            if mod_mut(l,2)<x0_min
                mod_mut(l,2) = x0_min;
            end
            if mod_mut(l,3)<z_min
                mod_mut(l,3) = z_min;
            end
            if mod_mut(l,4)<alpha_min
                mod_mut(l,4) = alpha_min;
            end
            if mod_mut(l,5)<a_min
                mod_mut(l,5) = a_min;
            end
            if mod_mut(l,1)>k_max 
                mod_mut(l,1) = k_max;
            end
            if mod_mut(l,2)>x0_max
                mod_mut(l,2) = x0_max;
            end
            if mod_mut(l,3)>z_max
                mod_mut(l,3) = z_max;
            end
            if mod_mut(l,4)>alpha_max
                mod_mut(l,4) = alpha_max;
            end
            if mod_mut(l,5)>a_max
                mod_mut(l,5) = a_max;
            end
        end
        %Hitung model untuk prosedur mutualisme
        for l = 1 : 2
            [Vcal_mut] = fwd_SP(mod_mut(l,1),mod_mut(l,2),mod_mut(l,3),mod_mut(l,4),mod_mut(l,5),x);
            err_mut = misfit_SP(V_obs,Vcal_mut);
            Em(l) = err_mut;
            %Update model jika  nilai misfit lebih baik proses mutualisme
            if l == 1
                if Em(l)<err(i)
                    model(i,:) = mod_mut(l,:);
                    err(i) = Em(l);
                    V_cal(i,:) = Vcal_mut;
                end
            else
                if Em(l)<err(j)
                    model(j,:) = mod_mut(l,:);
                    err(j) = Em(l);
                    V_cal(j,:) = Vcal_mut;
                end
            end
        end
        %Komensalisme
        j = randi(npop,1);
        if j == i
            j = randi(npop,1);
        end
        mod_com = model(i) +(0.4+0.9*rand)*(model_best-model(j));
            if mod_com(1)<k_min 
                mod_com(1) = k_min;
            end
            if mod_com(2)<x0_min
                mod_com(2) = x0_min;
            end
            if mod_com(3)<z_min
                mod_com(3) = z_min;
            end
            if mod_com(4)<alpha_min
                mod_com(4) = alpha_min;
            end
            if mod_com(5)<a_min
                mod_com(5) = a_min;
            end
            if mod_com(1)>k_max 
                mod_com(1) = k_max;
            end
            if mod_com(2)>x0_max
                mod_com(2) = x0_max;
            end
            if mod_com(3)>z_max
                mod_com(3) = z_max;
            end
            if mod_com(4)>alpha_max
                mod_com(4) = alpha_max;
            end
            if mod_com(5)>a_max
                mod_com(5) = a_max;
            end
        %Perhitungan misfit untuk prosedur komensalisme
        [Vcal_com] = fwd_SP(mod_com(1),mod_com(2),mod_com(3),mod_com(4),mod_com(5),x);
         Ec = misfit_SP(V_obs,Vcal_com);
         %Update model jika  nilai misfit lebih baik proses komensalisme
         if Ec < err(i)
             model(i,:) = mod_com(1,:);
             err(i) = Ec;
             V_cal(i,:) = Vcal_com(1,:);
         end
         %Parasitisme
         j = randi(npop,1);
         if j == i 
             j = randi(npop,1);
        end
         mod_par = model(i,:);
         p1 = randi(5,1);
         if p1 == 1
            mod_par(1) = k_min + rand*(k_max-k_min);
         elseif p1 == 2
            mod_par(2) = x0_min + rand*(x0_max-x0_min);
         elseif p1 == 3
            mod_par(3) = z_min + rand*(z_max-z_min);
         elseif p1 == 4
            mod_par(4) = alpha_min + rand*(alpha_max-alpha_min);
         else
            mod_par(5) = a_min + rand*(a_max-a_min);
         end
         %Perhitungan misfit untuk tahap parasitisme
         [Vcal_par] = fwd_SP(mod_par(1),mod_par(2),mod_par(3),mod_par(4),mod_par(5),x);
         Ep = misfit_SP(V_obs,Vcal_par);
         %Update model jika  nilai misfit lebih baik proses parasitisme
         if Ep < err(j)
             model(j,:) = mod_par(1,:);
             err(j) = Ep;
             V_cal(j,:) = Vcal_par(1,:);
         end
    end
    %Update model terbaik untuk setiap iterasi
    Emin = 100;
    for ipop = 1 : npop
        Emin = err(ipop);
        model_baru = model(ipop,:);
        V_model = V_cal(ipop,:);
    end
    %Nilai misfit terbaik
    Egen(itr)=Emin;

figure(1)
plot(x,V_obs,'r.',x,V_model,'b-','MarkerSize',20,'MarkerFaceColor','r','LineWidth',2.5);
axis([-100 100 -150 50])
legend({'Forward Model','Inversion Model'},'Location','Southeast')
title(['\bf \fontsize{10}\fontname{Times}SP Anomaly  || Misfit : ', num2str(Egen(itr)),' || iteration : ', num2str(itr)],'FontWeight','bold');
subtitle(['k = ',num2str(model_baru(1)),' || x0 = ',num2str(model_baru(2)),' || z = ',num2str(model_baru(3)),' || alpha = ',num2str(model_baru(4)),' || a = ',num2str(model_baru(5))],'FontWeight','bold')
xlabel('Distance (m)','FontWeight','bold');
ylabel('SP Anomaly (mV)','FontWeight','bold');
grid on
set(gcf, 'Position', get(0, 'Screensize'));
frame = getframe(gcf);
writeVideo(v,frame);
end
close(v);

%plot misfit
figure(2)
plot(1:nitr,Egen,'r','Linewidth',1.5)
xlabel('Iteration Number','FontSize',10,'FontWeight','Bold');
ylabel('misfit','FontSize',10,'FontWeight','Bold');
title('\bf \fontsize{12} Grafik Misfit ');
grid on
clc; clear;
set(0,'defaultTextInterpreter','latex');
set(gca,'nextplot','replacechildren');  set(gcf,'Renderer','zbuffer');  set(gca, 'FontName', 'Helvetica'); set(gca, 'ColorOrder', colormap(gray(12)));
global betai betaj gamma sigma alpha12 alpha21 delta nu a mu w beta alpha_All alpha_P12 alpha_P21 alpha_P
run importtraveldata.m
a=0;
w=0;
betai = 0.74288;
betaj = 0.53171;
sigmai = 0.21524;
sigmaj = 0.09358;
deltai = 0.08708;
deltaj = 0.10717;
gammaij = 0.51289;
gammaji = 0.51707;
nui = 0.29835;
nuj = 0.22016;

etai = 0.29014;
mui = 14.39439;
xii =0.59203;
etaj = 0.26917;
muj =17.98416;
xij=0.85078;
alpha_P12 = alpha_P(1,2);
alpha_P21 = alpha_P(2,1);
t1 = 7;

ft = 1:1:123; % Generate t for f 
beta_effL = betai*(etai+(1-etai)./(1+exp(xii*(ft-t1-mui))));
gt = 1:1:123; % Generate t for f 
beta_effY = betaj*(etaj+(1-etaj)./(1+exp(xij*(gt-t1-muj))));
day = ([1:123;]');
N1 = 795430; %total leeds population estimate 2020
N2= 211116; %total york population estimate 2020

i01 = 3; %no of total cases on 17/03
e01 = 1;
s01 = N1 - i01;
r01 = 0;
d01 = 1 ;% no of total deaths on 17/03

i02 = 2; %no of total cases on 17/03
e02 = 1;
s02 = N2 - i02;
r02 = 0;
d02 = 1; % no of total deaths on 17/03



alpha_12 = [alpha_P(1,2)*.75,alpha_P(1,2)*0.5,alpha_P(1,2)*0.25,alpha_P(1,2)*0.1,alpha_P(1,2)*0.05];
alpha_21 = [alpha_P(2,1)*0.75,alpha_P(2,1)*0.5,alpha_P(2,1)*0.25,alpha_P(2,1)*0.1,alpha_P(2,1)*0.05];
alpha_12A = [alpha_All(1,2)*0.75,alpha_All(1,2)*0.5,alpha_All(1,2)*0.25,alpha_All(1,2)*0.1,alpha_All(1,2)*0.05];
alpha_21A = [alpha_All(2,1)*0.75,alpha_All(2,1)*0.5,alpha_All(2,1)*0.25,alpha_All(2,1)*0.1,alpha_All(2,1)*0.05];
%% 
 kk = 1; 
  
 for ii = 1:5
   for jj = 1:5 
 for hh = 1:5
     for mm = 1:5
     alpha_P12 = alpha_12(ii);
     alpha_P21 = alpha_21(jj);
     alpha12 = alpha_12A(hh);
     alpha21 = alpha_21A(mm);
     kk    = kk + 1;

       tableA(kk-1,:)=[alpha_P12 alpha_P21 alpha12 alpha21];
  
   end 
 end
   end
 end
 %%
 k = 1;

  for i = [1,157,313,469,625]
     alpha_P12 = tableA(i,1);
     alpha_P21 = tableA(i,2);
     alpha12 = tableA(i,3);
     alpha21 = tableA(i,4);

TSPAN = [1:0.5:123]; % Solve from t=1 to t=5
IC = [s01,e01,i01,r01,d01,s02,e02,i02,r02,d02];

[T Y] = ode45(@(t,y) myODEsensitivity(t, y, ft, beta_effL, gt,beta_effY), TSPAN, IC); % Solve ODE

Si = Y(:,1);
Ei = Y(:,2);
Ii = Y(:,3);
Ri = Y(:,4);
Di = Y(:,5);
Sj = Y(:,6);
Ej = Y(:,7);
Ij = Y(:,8);
Rj = Y(:,9);
Dj = Y(:,10);



[Peak21, PeakIdx] = findpeaks(Ii); [Peak22, PeakIdx2] = findpeaks(Ij);
Peak21 = Peak21;   Peak22 = Peak22; 
PeakIdx = PeakIdx; PeakIdx2 = PeakIdx2;

peaks21 = [' Peak At t = ', num2str(T(PeakIdx), '%.2g'), ' Days with ', num2str(Peak21,'%.2g'), ' infections'];
peaks22 = [' Peak At t = ', num2str(T(PeakIdx2), '%.2g'), ' Days with ', num2str(Peak22, '%.2g'),' infections'];

if i == 1
    percenttxt = 25;
    txt = ['I_1, ' peaks21 ', ' num2str(percenttxt) '% Reduction in Travel'];
    txt2 = ['I_2, '  peaks22 ', ' num2str(percenttxt) '% Reduction in Travel'];
elseif i == 157
    percenttxt = 50;
    txt = ['I_1, ' peaks21 ', ' num2str(percenttxt) '% Reduction in Travel'];
    txt2 = ['I_2, '  peaks22 ', ' num2str(percenttxt) '% Reduction in Travel'];
elseif i == 313 
    percenttxt = 75;
    txt = ['I_1, ' peaks21 ', ' num2str(percenttxt) '% Reduction in Travel'];
    txt2 = ['I_2, '  peaks22 ', ' num2str(percenttxt) '% Reduction in Travel'];
elseif i == 469
    percenttxt = 90;
    txt = ['I_1, ' peaks21 ', ' num2str(percenttxt) '% Reduction in Travel'];
    txt2 = ['I_2, '  peaks22 ', ' num2str(percenttxt) '% Reduction in Travel'];
else
    percenttxt = 95;
    txt = ['I_1, ' peaks21 ', ' num2str(percenttxt) '% Reduction in Travel'];
    txt2 = ['I_2, '  peaks22 ', ' num2str(percenttxt) '% Reduction in Travel'];
end 

colororder({'#7a0177','#ae017e','#dd3497','#f768a1' ,'#fa9fb5','#fcc5c0','#feebe2'})
figure(1)

hold on
plot(T,Ii,'LineWidth',1.25,'DisplayName',txt)
xlabel('Time, Days','FontSize',30)
 ylabel('Number of Individuals in Infected Compartment in Leeds','interpreter','latex',FontSize=30)
xlim([0 123])
legend show
legend(Location="northeast", FontSize=20)

g = gcf;
g.WindowState = 'maximized';
ax = gca;
ax.XAxis.FontSize = 15;
ax.XLabel.FontSize = 30;
ax.YAxis.FontSize = 15;
ax.YLabel.FontSize = 30;

figure(2)

hold on
plot(T,Ij,'LineWidth',1.25,'DisplayName',txt2) 
xlim([0 123])
 
ylabel('Number of Individuals in Infected Compartment in York','interpreter','latex', FontSize=30)
xlabel('Time, Days','FontSize',30)
legend show
legend(Location="northeast", FontSize=20)

ax = gca;
ax.XAxis.FontSize = 15;
ax.XLabel.FontSize = 30;
ax.YAxis.FontSize = 15;
ax.YLabel.FontSize = 30;
g = gcf;
g.WindowState = 'maximized';
  end 

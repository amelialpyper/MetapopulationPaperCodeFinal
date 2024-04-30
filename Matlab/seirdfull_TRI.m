clc; clear;
set(0,'defaultTextInterpreter','latex');
set(gca,'nextplot','replacechildren');  set(gcf,'Renderer','zbuffer');  set(gca, 'FontName', 'Helvetica');
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
alpha12 = alpha_All(1,2);
alpha21 = alpha_All(2,1);
t1 = 7;

run r03_fullORIGINAL.m %make sure file path reflects where this file is saved

R0 = double(maxeig);
%%
ft = 1:1:123; % Generate t for f 
beta_effL = betai*(etai+(1-etai)./(1+exp(xii*(ft-t1-mui))));
gt = 1:1:123; % Generate t for f 
beta_effY = betaj*(etaj+(1-etaj)./(1+exp(xij*(gt-t1-muj))));
day = ([1:106;]');
%coviddataleedsyork.Day = day;
%day1 = coviddataleedsyork.Day;
N1 = 795425; %total leeds population estimate 2020
N2= 211113; %total york population estimate 2020

i01 = 3; %no of total cases on 17/03
e01 = 1;
s01 = N1 - i01-e01;
r01 = 0;
d01 = 1 ;% no of total deaths on 17/03

i02 = 2; %no of total cases on 17/03
e02 = 1;
s02 = N2 - i02-e01;
r02 = 0;
d02 = 1; % no of total deaths on 17/03



TSPAN = [1 123]; % Solve from t=1 to t=123
IC = [s01,e01,i01,r01,d01,s02,e02,i02,r02,d02];

[T Y] = ode45(@(t,y) myODE(t, y, ft, beta_effL, gt,beta_effY), TSPAN, IC); % Solve ODE

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
 peaks21 = ['Peak At t = ', num2str(T(PeakIdx), '%.2g'), ' Days with ', num2str(Peak21,'%.2g'), ' infections'];
 peaks22 = ['Peak At t = ', num2str(T(PeakIdx2), '%.2g'), ' Days with ', num2str(Peak22, '%.2g'),' infections'];
colororder({'#7a0177','#ae017e','#dd3497','#f768a1' ,'#fa9fb5','#fcc5c0','#feebe2'})
hold on

figure(1)

hold on
plot(T(PeakIdx),Peak21,'Marker','*','MarkerSize',15,'Color','black','linewidth',1.5,LineStyle='none')
plot(T,Ei,'-',T,Ii,'--',T,Ri,':',T,Di,'-.','LineWidth',3,'MarkerSize',.5)
ylabel('Number of Individuals in each \\ Compartment in Leeds','Color', 'black',FontSize=30)
xlabel('Time, Days', FontSize=30); xlim([0 123])
lg = legend({peaks21,'$E_1$','$I_1$','$R_1$','$D_1$'},'FontSize',25, 'NumColumns',1,'Location','northwest','Interpreter','latex','Box','on');
legend show
ax = gca;
ax.XAxis.FontSize = 15;
ax.XLabel.FontSize = 30;
ax.YAxis.FontSize = 15;
ax.YLabel.FontSize = 30;
g = gcf;
g.WindowState = 'maximized';

figure(2)
colororder({'#7a0177','#ae017e','#dd3497','#f768a1' ,'#fa9fb5','#fcc5c0','#feebe2'})
hold on 
plot(T(PeakIdx2),Peak22,'Marker','.','MarkerSize',40,'Color','black',LineStyle='none')
plot(T,Ej,'-',T,Ij,'--',T,Rj,':',T,Dj,'-.','LineWidth',3,'MarkerSize',.5)
xlabel('Time, Days', FontSize=30); ylabel('Number of Individuals in each Compartment in York', FontSize=30)
xlim([0 123])
legend({peaks22,'$E_2$','$I_2$','$R_2$','$D_2$'},'FontSize',25, 'NumColumns',1,'Location','northwest','Interpreter','latex','Box','on'); 
legend show
ax = gca;
ax.XAxis.FontSize = 15;
ax.XLabel.FontSize = 30;
ax.YAxis.FontSize = 15;
ax.YLabel.FontSize = 30;
g = gcf;
g.WindowState = 'maximized';

%Save images to current wd
%exportgraphics(figure(1), seirdfullmeta_leeds.pdf','Resolution',300)
%exportgraphics(figure(2), 'seirdfullmeta_york.pdf','Resolution',300)

function dydt = myODE(t, x, ft, beta_effL,gt, beta_effY)
run importdata.m
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
alpha12 = alpha_All(1,2);
alpha21 = alpha_All(2,1);

beta_effL = interp1(ft, beta_effL, t); % Interpolate the data set (ft, f) at times t %leeds
beta_effY = interp1(gt, beta_effY, t); % Interpolate the data set (ft, f) at times t %york
Si = x(1);
Ei = x(2);
Ii = x(3);
Ri = x(4);
Di = x(5);
Sj = x(6);
Ej= x(7);
Ij = x(8);
Rj = x(9);
Dj = x(10);

Ni = Si + Ei + Ii + Ri;
Nj = Sj + Ej + Ij + Rj;

dsi = a*Ni-((beta_effL)*Si*Ii)/Ni - gammaij*alpha_P(1,2)*(Si)*((alpha_P(1,2)*Ii/Ni)+(alpha_P(2,1)*Ij/Nj)) - (beta_effL)*alpha21*(Si)*(Ij/Nj) - (beta_effY)*alpha12*(Si)*(Ij/Nj) - beta_effY*(alpha12)^2*Si*(Ii/Ni)+ w*Ri;
dei = (beta_effL)*Si*(Ii/Ni) + gammaij*alpha_P(1,2)*(Si)*((alpha_P(1,2)*Ii/Ni)+(alpha_P(2,1)*Ij/Nj)) + (beta_effL)*alpha21*(Si)*(Ij/Nj) + (beta_effY)*alpha12*(Si)*(Ij/Nj) + beta_effY*(alpha12)^2*Si*(Ii/Ni) - nui*Ei;
dii = nui*Ei - deltai*Ii - (sigmai)*Ii;
dri = deltai*Ii- w*Ri;
dni = dsi + dei + dii + dri;
ddi = -dni;

dsj = a*Nj -((beta_effY)*Sj*Ij)/Nj - gammaji*alpha_P(2,1)*(Sj)*((alpha_P(1,2)*Ii/Ni)+(alpha_P(2,1)*Ij/Nj)) - (beta_effY)*alpha12*(Sj)*(Ii/Ni) - (beta_effL)*alpha21*(Sj)*(Ii/Ni) - beta_effL*(alpha21)^2*Sj*(Ij/Nj) + w*Rj;
dej = ((beta_effY)*Sj*Ij)/Nj + gammaji*alpha_P(2,1)*(Sj)*((alpha_P(1,2)*Ii/Ni)+(alpha_P(2,1)*Ij/Nj)) + (beta_effY)*alpha12*(Sj)*(Ii/Ni) + (beta_effL)*alpha21*(Sj)*(Ii/Ni)  + beta_effL*(alpha21)^2*Sj*(Ij/Nj) - nuj*Ej; 
dij = nuj*Ej - deltaj*Ij  - (sigmaj)*Ij;
drj = deltaj*Ij - w*Rj;
dnj = dsj + dej + dij + drj;
ddj = -dnj;

dydt = [dsi;dei;dii;dri;ddi;dsj;dej;dij;drj;ddj];
end
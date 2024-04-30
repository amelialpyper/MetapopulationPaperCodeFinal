clc; clear;
set(0,'defaultTextInterpreter','latex');
set(gca,'nextplot','replacechildren');  set(gcf,'Renderer','zbuffer');  set(gca, 'FontName', 'Helvetica');
run importtraveldata.m
a=0;
w=0;
betai = 0.25; 
betaj = 0.25;
sigmai =  0.025;
sigmaj = 0.025;
deltai = 1/10;
deltaj = 1/10;
nui = 1/5.2;
nuj = 1/5.2;
alpha_P12 = 0;
alpha_P21 = 0;
alpha12 = 0;
alpha21 = 0;
gammaij=0;
gammaji=0;

run r03_fullORIGINAL.m

R0 = double(maxeig);
%%
N1 = 795425; %total leeds population estimate 2020
N2= 211116; %total york population estimate 2020

i01 = 3; %no of total cases on 17/03
e01 = 1;
s01 = N1 - i01-e01;
r01 = 0;
d01 = 0 ;% no of total deaths on 17/03

i02 = 2; %no of total cases on 17/03
e02 = 1;
s02 = N2 - i02-e02;
r02 = 0;
d02 = 0; % no of total deaths on 17/03



TSPAN = [1 365]; % Solve from t=1 to t=123
IC = [s01,e01,i01,r01,d01,s02,e02,i02,r02,d02];

[T Y] = ode45(@(t,y) myODE(t, y), TSPAN, IC); % Solve ODE

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
 peaks21 = ['Peak At t = ', num2str(T(PeakIdx), '%.3g'), ' Days with ', num2str(Peak21,'%.5g'), ' infections'];
 peaks22 = ['Peak At t = ', num2str(T(PeakIdx2), '%.3g'), ' Days with ', num2str(Peak22, '%.5g'),' infections'];
 %%
colororder({'#7a0177','#ae017e','#dd3497','#f768a1' ,'#fa9fb5','#fcc5c0','#feebe2'})
hold on

hold on
figure(1)
plot(T(PeakIdx),Peak21,'Marker','*','MarkerSize',15,'Color','black','linewidth',1.5,LineStyle='none')
plot(T,Ei,'-',T,Ii,'--',T,Ri,':',T,Di,'-.','LineWidth',3,'MarkerSize',.5)
ylabel({'Number of Individuals in each Compartment in City 1'},'Color', 'black',FontSize =30, Interpreter='latex')
xlabel('Time, Days', FontSize=30); xlim([0 365])
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

xlabel('Time, Days', FontSize=30); ylabel({'Number of Individuals in each Compartment in City 2'}, FontSize=30, Interpreter='latex')

xlim([0 365])
legend({peaks22,'$E_2$','$I_2$','$R_2$','$D_2$'},'FontSize',25, 'NumColumns',1,'Location','northwest','Interpreter','latex','Box','on'); 
legend show
ax = gca;
ax.XAxis.FontSize = 15;
ax.XLabel.FontSize = 30;
ax.YAxis.FontSize = 15;
ax.YLabel.FontSize = 30;
 g = gcf;
 g.WindowState = 'maximized';

%exportgraphics(figure(1), 'seirdtoy_leeds_notrav.pdf','Resolution',300)
%exportgraphics(figure(2), 'seirdtoy_york_notrav.pdf','Resolution',300)

function dydt = myODE(t, x)
run importtraveldata.m
a=0;
w=0;
betai = 0.25; 
betaj = 0.25;
sigmai =  0.025;
sigmaj = 0.025;
deltai = 1/10;
deltaj = 1/10;
nui = 1/5.2;
nuj = 1/5.2;
alpha_P12 = 0;
alpha_P21 = 0;
alpha12 = 0;
alpha21 = 0;
gammaij=0;
gammaji=0;

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

dsi = a*Ni-((betai)*Si*Ii)/Ni - gammaij*alpha_P(1,2)*(Si)*((alpha_P(1,2)*Ii/Ni)+(alpha_P(2,1)*Ij/Nj)) - (betai)*alpha21*(Si)*(Ij/Nj) - (betaj)*alpha12*(Si)*(Ij/Nj) - betaj*(alpha12)^2*Si*(Ii/Ni)+ w*Ri;
dei = (betai)*Si*(Ii/Ni) + gammaij*alpha_P(1,2)*(Si)*((alpha_P(1,2)*Ii/Ni)+(alpha_P(2,1)*Ij/Nj)) + (betai)*alpha21*(Si)*(Ij/Nj) + (betaj)*alpha12*(Si)*(Ij/Nj) + betaj*(alpha12)^2*Si*(Ii/Ni) - nui*Ei;
dii = nui*Ei - deltai*Ii - (sigmai)*Ii;
dri = deltai*Ii- w*Ri;
dni = dsi + dei + dii + dri;
ddi = -dni;

dsj = a*Nj -((betaj)*Sj*Ij)/Nj - gammaji*alpha_P(2,1)*(Sj)*((alpha_P(1,2)*Ii/Ni)+(alpha_P(2,1)*Ij/Nj)) - (betaj)*alpha12*(Sj)*(Ii/Ni) - (betai)*alpha21*(Sj)*(Ii/Ni) - betai*(alpha21)^2*Sj*(Ij/Nj) + w*Rj;
dej = ((betaj)*Sj*Ij)/Nj + gammaji*alpha_P(2,1)*(Sj)*((alpha_P(1,2)*Ii/Ni)+(alpha_P(2,1)*Ij/Nj)) + (betaj)*alpha12*(Sj)*(Ii/Ni) + (betai)*alpha21*(Sj)*(Ii/Ni)  + betai*(alpha21)^2*Sj*(Ij/Nj) - nuj*Ej; 
dij = nuj*Ej - deltaj*Ij  - (sigmaj)*Ij;
drj = deltaj*Ij - w*Rj;
dnj = dsj + dej + dij + drj;
ddj = -dnj;

dydt = [dsi;dei;dii;dri;ddi;dsj;dej;dij;drj;ddj];
end
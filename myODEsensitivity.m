function dydt = myODEsensitivity(t, x, ft, beta_effL,gt, beta_effY)
global betai betaj gammaij gammaji sigmai alpha12 alpha21 deltai nui alpha_P12 alpha_P21 
run importtraveldata.m
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
a=0;
w=0;


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

dsi = a*Ni-((beta_effL)*Si*Ii)/Ni - gammaij*alpha_P12*(Si)*((alpha_P12*Ii/Ni)+(alpha_P21*Ij/Nj)) - (beta_effL)*alpha21*(Si)*(Ij/Nj) - (beta_effY)*alpha12*(Si)*(Ij/Nj) - beta_effY*(alpha12)^2*Si*(Ii/Ni)+ w*Ri;
dei = (beta_effL)*Si*(Ii/Ni) + gammaij*alpha_P12*(Si)*((alpha_P12*Ii/Ni)+(alpha_P21*Ij/Nj)) + (beta_effL)*alpha21*(Si)*(Ij/Nj) + (beta_effY)*alpha12*(Si)*(Ij/Nj) + beta_effY*(alpha12)^2*Si*(Ii/Ni) - nui*Ei;
dii = nui*Ei - deltai*Ii - (sigmai)*Ii;
dri = deltai*Ii- w*Ri;
dni = dsi + dei + dii + dri;
ddi = -dni;

dsj = a*Nj -((beta_effY)*Sj*Ij)/Nj - gammaji*alpha_P21*(Sj)*((alpha_P12*Ii/Ni)+(alpha_P12*Ij/Nj)) - (beta_effY)*alpha12*(Sj)*(Ii/Ni) - (beta_effL)*alpha21*(Sj)*(Ii/Ni) - beta_effL*(alpha21)^2*Sj*(Ij/Nj) + w*Rj;
dej = ((beta_effY)*Sj*Ij)/Nj + gammaji*alpha_P21*(Sj)*((alpha_P(1,2)*Ii/Ni)+(alpha_P12*Ij/Nj)) + (beta_effY)*alpha12*(Sj)*(Ii/Ni) + (beta_effL)*alpha21*(Sj)*(Ii/Ni)  + beta_effL*(alpha21)^2*Sj*(Ij/Nj) - nuj*Ej; 
dij = nuj*Ej - deltaj*Ij  - (sigmaj)*Ij;
drj = deltaj*Ij - w*Rj;
dnj = dsj + dej + dij + drj;
ddj = -dnj;

dydt = [dsi;dei;dii;dri;ddi;dsj;dej;dij;drj;ddj];
end

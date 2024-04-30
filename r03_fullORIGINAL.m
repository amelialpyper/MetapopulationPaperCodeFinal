%; clear;
% global alpha12 alpha21 a
syms Si Ei Ii Ri Sj Ej Ij Rj Ni Nj
alphaA12 = alpha12;
alphaA21=alpha21;
alphaP12=alpha_P12;
alphaP21=alpha_P21;
F3 = [(betai)*Si*(Ii/Ni) + (betai)*alphaA21*(Si)*(Ij/Nj) + (betaj)*alphaA12*(Si)*(Ij/Nj) + betaj*(alphaA12)^2*Si*(Ii/Ni)+ gammaij*alphaP12*(Si)*((alphaP12*Ii/Ni)+(alphaP21*Ij/Nj)), 0,((betaj)*Sj*Ij)/Nj + gammaji*alphaP21*(Sj)*((alphaP12*Ii/Ni)+(alphaP21*Ij/Nj)) + (betaj)*alphaA12*(Sj)*(Ii/Ni) + (betai)*alphaA21*(Sj)*(Ii/Ni)  + betai*(alphaA21)^2*Sj*(Ij/Nj),0];
V3 = [nui*Ei, -nui*Ei + deltai*Ii + sigmai*Ii, nuj*Ej, -nuj*Ej + deltaj*Ij + sigmaj*Ij];
Fbold3 = jacobian(F3,[Ei,Ii,Ej,Ij]);
FBold3 = subs(Fbold3,[Si,Sj],[Ni,Nj]);
%%
Vbold3 = jacobian(V3, [Ei,Ii,Ej,Ij]);
Vinv3 = inv(Vbold3);

G3 = FBold3*Vinv3;
E3 = eig(G3);
maxeig = E3(4,:);

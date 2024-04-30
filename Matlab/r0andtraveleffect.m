clc; clear;
set(0,'defaultTextInterpreter','latex');
set(gca,'nextplot','replacechildren');  set(gcf,'Renderer','zbuffer');  set(gca, 'FontName', 'Helvetica');
global betai betaj gammaij gammaji sigmai sigmaj alpha12 alpha21 deltai deltaj nui nuj a mu w beta alpha_All alpha_P12 alpha_P21 alpha_P alphaA21 alphaA12 alphaP12 alphaP21
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


alpha_12 = [alpha_P(1,2),alpha_P(1,2)*0.5,alpha_P(1,2)*0.25,alpha_P(1,2)*0.1,alpha_P(1,2)*0.05];
alpha_21 = [alpha_P(2,1),alpha_P(2,1)*0.5,alpha_P(2,1)*0.25,alpha_P(2,1)*0.1,alpha_P(2,1)*0.05];
alpha_12A = [alpha_All(1,2),alpha_All(1,2)*0.5,alpha_All(1,2)*0.25,alpha_All(1,2)*0.1,alpha_All(1,2)*0.05];
alpha_21A = [alpha_All(2,1),alpha_All(2,1)*0.5,alpha_All(2,1)*0.25,alpha_All(2,1)*0.1,alpha_All(2,1)*0.05];

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

x = []; 
y = [];
ll=1;

  for i = [1,152,313,469,625]
     alpha_P12 = tableA(i,1);
     alpha_P21 = tableA(i,2);
     alpha12 = tableA(i,3);
     alpha21 = tableA(i,4);
     ll = ll + 1;

tableB(ll-1,:) = [alpha_P12 alpha_P21 alpha12 alpha21];

for betai = [0.3]%0.295%0.3
    for betaj = [0.2]%0.195%0.2
       run r03_fullORIGINAL.m;
     r03 = subs(maxeig,{alphaP12,alphaP21,alphaA12,alphaA21,betai,betaj,deltai,deltaj,sigmai,sigmaj,nui,nuj,gammaij,gammaji},{alpha_P12,alpha_P21,alpha12,alpha21,betai,betaj,deltai,deltaj,sigmai,sigmaj,nui,nuj,gammaij,gammaji});
     r03 = double(r03);
     x = [x,r03];
    end  
end
  end
%%
%For Legend Purposes      
hold on
txt1 = ['R_0 = ',num2str(x(:,1), '%.4g')];
txt2 = ['R_0 = ',num2str(x(:,2), '%.4g')]; %For beta = 0.3, 0.2
%txt2 = ['R_0 = ',num2str(x(:,2), '%.3g')];% For beta = 0.295, 0.195
txt3 = ['R_0 = ',num2str(x(:,3), '%.4g')]; %For beta = 0.3, 0.2
%txt3 = ['R_0 = ',num2str(x(:,3), '%.3g')]; % For beta = 0.295, 0.195
txt4 = ['R_0 = ',num2str(x(:,4), '%.3g')];
txt5 = ['R_0 = ',num2str(x(:,5), '%.3g')];


  hold on 
  colororder({'#49006a','#7a0177','#ae017e','#dd3497' ,'#f768a1','#fa9fb5','#fcc5c0'})
  figure(1)
      plot(tableB(1,1),x(:,1),marker = "*",MarkerSize=20,DisplayName=txt1,LineStyle="none")
      plot(tableB(2,1),x(:,2),marker = ".", MarkerSize=30,DisplayName=txt2,LineStyle="none")
      plot(tableB(3,1),x(:,3),marker = "o", MarkerSize=15,DisplayName=txt3,LineStyle="none")
      plot(tableB(4,1),x(:,4),marker = "+", MarkerSize=20,DisplayName=txt4,LineStyle="none")
      plot(tableB(5,1),x(:,5),marker = "x", MarkerSize=20,DisplayName=txt5,LineStyle="none")
       legend show
  legend(FontSize=20, Location="best")
    ax = gca;
  ax.XAxis.FontSize = 12;
ax.XLabel.FontSize = 30;
ax.YAxis.FontSize = 12;
ax.YLabel.FontSize = 30;
  yline(1,'-','R_0 Threshold','HandleVisibility','off',FontSize=22,LabelHorizontalAlignment='left')
  xlabel('Percentage (\%) Reduction of Travel',FontSize=30)
ylim([0.97 1.03])
  ylabel('$R_0$ Value',FontSize=30)
   values = [3.09387412922414e-05 6.18774825844827e-05 0.000154693706461207 0.000309387412922414 0.000618774825844827];
   labels = ({'95%','90%','75%','50%','0%'});
   set(gca, 'XTick',sort(values), 'XDir','reverse','XTickLabel',labels)
   set(gcf, 'Position',  [100, 100, 800, 700])
   %exportgraphics(figure(1), 'r0changing.pdf','Resolution',300)

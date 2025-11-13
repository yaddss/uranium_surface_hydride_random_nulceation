clc;
clear;
n = 40; 
x = linspace(0,1,n);
y = linspace(0,1,n);
[c,eta]=meshgrid(x,y);
G1 = 1.0e+02*[2.350916205236954  -0.146697171206786   0.002288475870826];
G2 = 1.0e+02*[3.687711694489345  -5.568444658678911   2.074807827005719];
a1 = G1(1) ;%; 
b1 = G1(2) ;%; 
c1 = G1(3) ;%; 
a2 = G2(1) ;%; 
b2 = G2(2) ;%; 
c2 = G2(3) ;%;
w = 48;
p = 20;
h_eta = eta.^3.*(6*eta.^2-15.*eta+10);
g_eta = 2*w*eta.^2.*(1-eta).^2;
cm = (c-h_eta*(b1-b2)/2/a2)./(1-h_eta+a1/a2*h_eta);
ch = (c-(1-h_eta)*(b2-b1)/2/a1)./(h_eta+a2/a1*(1-h_eta));
F1 = (G1(1)*cm.^2+G1(2)*cm+G1(3))/w;  %metal matrix
F2 = (G2(1)*ch.^2+G2(2)*ch+G2(3))/w;  %hydride precipit  ate
fbulk = (1-h_eta).*F1+h_eta.*F2+g_eta;
ei00 = [0.1259 0 0; 0 0.1259 0; 0 0 0.3777];
C11m = 207.6; %J/mm^3
C12m = 50.6; %J/mm^3
C22m = 185.9; %J/mm^3
C44m = 108.4; %J/mm^3
C11h = 222; %J/mm^3
C22h = C11h; %J/mm^3
C12h = 70; %J/mm^3
C44h = 58; %J/mm^3
C11 = C11m*(1-h_eta)+C11h*h_eta;
C22 = C22m*(1-h_eta)+C22h*h_eta;
C12 = C12m*(1-h_eta)+C12h*h_eta;
C44 = C44m*(1-h_eta)+C44h*h_eta;
epsilon1 = [-0.3 0; 0 0];
epsilon2 = [0 0; 0 0];
epsilon3 = [0.3 0; 0 0];
sigma_1 = [C11*epsilon1(1,1)+C12*epsilon1(2,2) 2*C44*epsilon1(1,2); 2*C44*epsilon1(1,2) C22*epsilon1(2,2)+C12*epsilon1(1,1)];
sigma_2 = [C11*epsilon2(1,1)+C12*epsilon2(2,2) 2*C44*epsilon2(1,2); 2*C44*epsilon2(1,2) C22*epsilon2(2,2)+C12*epsilon2(1,1)];
sigma_3 = [C11*epsilon3(1,1)+C12*epsilon3(2,2) 2*C44*epsilon3(1,2); 2*C44*epsilon3(1,2) C22*epsilon3(2,2)+C12*epsilon3(1,1)];
fel_1 = -(sigma_1(1,1)*ei00(1,1)+sigma_1(1,2)*ei00(1,2)+sigma_1(2,1)*ei00(2,1)+sigma_1(2,2)*ei00(2,2))*h_eta/w;
fel_2 = -(sigma_2(1,1)*ei00(1,1)+sigma_2(1,2)*ei00(1,2)+sigma_2(2,1)*ei00(2,1)+sigma_2(2,2)*ei00(2,2))*h_eta/w;
fel_3 = -(sigma_3(1,1)*ei00(1,1)+sigma_3(1,2)*ei00(1,2)+sigma_3(2,1)*ei00(2,1)+sigma_3(2,2)*ei00(2,2))*h_eta/w;
fbulk_plus_el_1 = fbulk+fel_1;
fbulk_plus_el_2 = fbulk+fel_2;
fbulk_plus_el_3 = fbulk+fel_3;

%%%%%%%%%%%%%%%%%%constraint c=[1-h_eta]c_m+h_eta*c_h%%%%%%%%%%%%%%%%%%%%%
c0 = linspace(0.05,0.75,n);
eta0 = linspace(0,1,n);
h_eta0 = eta0.^3.*(6*eta0.^2-15.*eta0+10);
g_eta0 = 2*w*eta0.^2.*(1-eta0).^2;
cm0 = (c0-h_eta0*(b1-b2)/2/a2)./(1-h_eta0+a1/a2*h_eta0);
ch0 = (c0-(1-h_eta0)*(b2-b1)/2/a1)./(h_eta0+a2/a1*(1-h_eta0));
F1_line = (G1(1)*cm0.^2+G1(2)*cm0+G1(3))/w;  %metal matrix
F2_line = (G2(1)*ch0.^2+G2(2)*ch0+G2(3))/w;  %hydride precipit  ate
fbulk_line = (1-h_eta0).*F1_line+h_eta0.*F2_line+g_eta0;

fel_1_line = -(sigma_1(1,1)*ei00(1,1)+sigma_1(1,2)*ei00(1,2)+sigma_1(2,1)*ei00(2,1)+sigma_1(2,2)*ei00(2,2))*h_eta0/w;
fbulk_plus_el_1_line = fbulk_line+fel_1_line;

fel_2_line = -(sigma_2(1,1)*ei00(1,1)+sigma_2(1,2)*ei00(1,2)+sigma_2(2,1)*ei00(2,1)+sigma_2(2,2)*ei00(2,2))*h_eta0/w;
fbulk_plus_el_2_line = fbulk_line+fel_2_line;

fel_3_line = -(sigma_3(1,1)*ei00(1,1)+sigma_3(1,2)*ei00(1,2)+sigma_3(2,1)*ei00(2,1)+sigma_3(2,2)*ei00(2,2))*h_eta0/w;
fbulk_plus_el_3_line = fbulk_line+fel_3_line;


% c_eq_mh_m = 0;
% c_eq_mh_h = 0.7499;
% A = [c_eq_mh_m,0,(G1(1)*c_eq_mh_m.^2+G1(2)*c_eq_mh_m+G1(3))/w];
% B = [c_eq_mh_h,1,(G2(1)*c_eq_mh_h.^2+G2(2)*c_eq_mh_h+G2(3))/w];

%green
C1(:,:,1) = 0*ones(n); % red 128
C1(:,:,2) = 0.5*ones(n); % green 0
C1(:,:,3) = 0*ones(n); % blue 128
%blue
C2(:,:,1) = 0.5*ones(n); % red 128
C2(:,:,2) = 0*ones(n); % green 0
C2(:,:,3) = 0.5*ones(n); % blue 128
%red
C3(:,:,1) = 1*ones(n); % red 255
C3(:,:,2) = 0*ones(n); % green 0
C3(:,:,3) = 0*ones(n); % blue 0

s1=surf(c,eta,fbulk_plus_el_1,C1);
hold on;
s2=surf(c,eta,fbulk_plus_el_2,C2);
hold on;
s3=surf(c,eta,fbulk_plus_el_3,C3);
shading interp;
s1.FaceAlpha = 0.7; %transparency
s2.FaceAlpha = 0.7; %transparency
s3.FaceAlpha = 0.7; %transparency
s1.EdgeColor = 'none';
s2.EdgeColor = 'none';
s3.EdgeColor = 'none';
hold on;
plot3(c0,eta0,fbulk_plus_el_1_line,'--','Color','#000000','LineWidth',3);%green #008000
hold on;
plot3(c0,eta0,fbulk_plus_el_2_line,'--','Color','#000000','LineWidth',3); %blue #0000FF
hold on;
plot3(c0,eta0,fbulk_plus_el_3_line,'--','Color','#000000','LineWidth',3);%red #FF0000
%hold on;
%plot3(c0,eta0,fbulk_line,'--','Color','#000000','LineWidth',3);%purple #800080
%plot3([A(1),B(1)],[A(2),B(2)],[A(3),B(3)],'-','Color','#008000','LineWidth',2);
hold on;
x_constraint = [0.05 0.05   0.75 0.75];
y_constraint = [0     0      1     1 ];
z_constraint = [100 -100  -100  100];
f1=fill3(x_constraint,y_constraint,z_constraint,'k');
f1.FaceAlpha = 0.8;
f1.FaceColor = '#A5A5A5';
f1.EdgeColor = 'none';


xlabel('H Concentration $c$','Fontname', 'Times New Roman','FontSize',20,'Interpreter','latex');
ylabel('Phase $\eta$','Fontname', 'Times New Roman','FontSize',20,'Interpreter','latex');
zlabel('Free energy density ($\omega$)','Fontname', 'Times New Roman','FontSize',20,'Interpreter','latex');
h=legend('$f_{\mathrm{bulk}}+f_{\mathrm{el}}(\varepsilon^a_{11} = -0.3)$','$f_{\mathrm{bulk}}+f_{\mathrm{el}}(\varepsilon^a_{11} = 0)$','$f_{\mathrm{bulk}}+f_{\mathrm{el}}(\varepsilon^a_{11} = 0.3)$','Fontname', 'Times New Roman','FontSize',16,'Interpreter','latex','position',[0.478955718927485,0.679777453289303,0.420605493851518,0.139688604782947]);
set(h,'Box','off')
axis([0 0.75 0 1.0 -2 10]);
set(gca,'fontsize',16,'fontname','Times');
set(gca,'XTick',(0:0.25:0.75));
set(gca,'YTick',(0:0.5:1));
set(gca,'xminortick','on');
set(gca,'yminortick','on');
set(gca,'position',[0.220275530586886,0.11,0.714724469413114,0.805])
box on;
ax = gca;
ax.BoxStyle = 'full';
view(45,10);
set(gcf,'unit','centimeters','position',[10 1 15 18])
title('$KKS$ Constraint','Fontname','Times New Roman','FontSize',20,'Interpreter','latex','position',[0.331107978123239,0.562468934321458,9.77801461517494]);
%set(gca,'ZScale','log')



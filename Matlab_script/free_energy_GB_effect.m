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
p =1;
h_eta = eta.^3.*(6*eta.^2-15.*eta+10);
g_eta = 2*w*eta.^2.*(1-eta).^2;
cm = (c-h_eta*(b1-b2)/2/a2)./(1-h_eta+a1/a2*h_eta);
ch = (c-(1-h_eta)*(b2-b1)/2/a1)./(h_eta+a2/a1*(1-h_eta));
F1 = (G1(1)*cm.^2+G1(2)*cm+G1(3))/w;  %metal matrix
F2 = (G2(1)*ch.^2+G2(2)*ch+G2(3))/w;  %hydride precipit  ate
fbulk = (1-h_eta).*F1+h_eta.*F2+g_eta;
theta1 = 1;
theta2 = 5;
theta3 = 15;
theta0 = 15;
gamma0 = 4;
gr1 = 0.5;
zeta1 = (gr1^2*(1-gr1)^2*theta1)/(gr1^2*(1-gr1)^2);
gamma1 =gamma0*(abs(theta1)/theta0).*(1-log(abs(theta1)/theta0)).*(abs(theta1)<=theta0)+(abs(theta1)>theta0);
fsg1 = -p*c*((-0.5*gr1^2+0.25*gr1^4)+(-0.5*(1-gr1)^2+0.25*(1-gr1)^4)+gamma1*(gr1^2*(1-gr1)^2)+0.25);
fbulk_plus_sg1 = fbulk+fsg1;

gr2 = 0.5;
zeta2 = (gr2^2*(1-gr2)^2*theta2)/(gr2^2*(1-gr2)^2);
gamma2 =gamma0*(abs(theta2)/theta0).*(1-log(abs(theta2)/theta0)).*(abs(theta2)<=theta0)+(abs(theta2)>theta0);
fsg2 = -p*c*((-0.5*gr2^2+0.25*gr2^4)+(-0.5*(1-gr2)^2+0.25*(1-gr2)^4)+gamma2*(gr2^2*(1-gr2)^2)+0.25);
fbulk_plus_sg2 = fbulk+fsg2;

gr3 = 0.5;
zeta3 = (gr3^2*(1-gr3)^2*theta3)/(gr3^2*(1-gr3)^2);
gamma3 =gamma0*(abs(theta3)/theta0).*(1-log(abs(theta3)/theta0)).*(abs(theta3)<=theta0)+(abs(theta3)>theta0);
fsg3 = -p*c*((-0.5*gr3^2+0.25*gr3^4)+(-0.5*(1-gr3)^2+0.25*(1-gr3)^4)+gamma3*(gr3^2*(1-gr3)^2)+0.25);
fbulk_plus_sg3 = fbulk+fsg3;


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

fsg1_line = -p*c0*((-0.5*gr1^2+0.25*gr1^4)+(-0.5*(1-gr1)^2+0.25*(1-gr1)^4)+gamma1*(gr1^2*(1-gr1)^2)+0.25);
fbulk_plus_sg1_line = fbulk_line+fsg1_line;

fsg2_line = -p*c0*((-0.5*gr2^2+0.25*gr2^4)+(-0.5*(1-gr2)^2+0.25*(1-gr2)^4)+gamma2*(gr2^2*(1-gr2)^2)+0.25);
fbulk_plus_sg2_line = fbulk_line+fsg2_line;

fsg3_line = -p*c0*((-0.5*gr3^2+0.25*gr3^4)+(-0.5*(1-gr3)^2+0.25*(1-gr3)^4)+gamma3*(gr3^2*(1-gr3)^2)+0.25);
fbulk_plus_sg3_line = fbulk_line+fsg3_line;


% c_eq_mh_m = 0;
% c_eq_mh_h = 0.7499;
% A = [c_eq_mh_m,0,(G1(1)*c_eq_mh_m.^2+G1(2)*c_eq_mh_m+G1(3))/w];
% B = [c_eq_mh_h,1,(G2(1)*c_eq_mh_h.^2+G2(2)*c_eq_mh_h+G2(3))/w];

%green
C1(:,:,1) = 0*ones(n); % red 128
C1(:,:,2) = 0.5*ones(n); % green 0
C1(:,:,3) = 0*ones(n); % blue 128
%blue
C2(:,:,1) = 0*ones(n); % red 0
C2(:,:,2) = 0*ones(n); % green 0
C2(:,:,3) = 1*ones(n); % blue 255
%red
C3(:,:,1) = 1*ones(n); % red 255
C3(:,:,2) = 0*ones(n); % green 0
C3(:,:,3) = 0*ones(n); % blue 0
%purple
C4(:,:,1) = 0.5*ones(n); % red 128
C4(:,:,2) = 0*ones(n); % green 0
C4(:,:,3) = 0.5*ones(n); % blue 128

s1=surf(c,eta,fbulk_plus_sg1,C1);
hold on;
s2=surf(c,eta,fbulk_plus_sg2,C2);
hold on;
s3=surf(c,eta,fbulk_plus_sg3,C3);
hold on;
s4=surf(c,eta,fbulk,C4);
shading interp;
s1.FaceAlpha = 0.7; %transparency
s2.FaceAlpha = 0.7; %transparency
s3.FaceAlpha = 0.7; %transparency
s4.FaceAlpha = 0.7; %transparency
s1.EdgeColor = 'none';
s2.EdgeColor = 'none';
s3.EdgeColor = 'none';
s4.EdgeColor = 'none';
hold on;
plot3(c0,eta0,fbulk_plus_sg1_line,'--','Color','#000000','LineWidth',3);%green #008000
hold on;
plot3(c0,eta0,fbulk_plus_sg2_line,'--','Color','#000000','LineWidth',3); %blue #0000FF
hold on;
plot3(c0,eta0,fbulk_plus_sg3_line,'--','Color','#000000','LineWidth',3);%red #FF0000
hold on;
plot3(c0,eta0,fbulk_line,'--','Color','#000000','LineWidth',3);%purple #800080
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
h=legend('$f_{\mathrm{bulk}}+f_{\mathrm{sg}}(\Delta\theta_{gg^\prime}=\pm1^\circ)$','$f_{\mathrm{bulk}}+f_{\mathrm{sg}}(\Delta\theta_{gg^\prime}=\pm5^\circ)$','$f_{\mathrm{bulk}}+f_{\mathrm{sg}}(\Delta\theta_{gg^\prime}=\pm15^\circ)$','$f_{\mathrm{bulk}}$','Fontname', 'Times New Roman','FontSize',16,'Interpreter','latex','position',[0.478955718927485,0.679777453289303,0.420605493851518,0.139688604782947]);
set(h,'Box','off')
%axis([0.7 0.75 0.95 1.0 -0.3 0.3]);
axis([0 0.75 0 1.0 -5 10]);
set(gca,'fontsize',16,'fontname','Times');
%set(gca,'XTick',(0.7:0.025:0.75));
%set(gca,'YTick',(0.95:0.025:1));
%set(gca,'ZTick',(-0.3:0.1:0.3));
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



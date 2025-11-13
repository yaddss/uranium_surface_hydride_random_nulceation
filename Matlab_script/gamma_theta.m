clc;
clear;
theta = [linspace(-20,-2,90) linspace(-2,2,70) linspace(2,20,90)];
theta_sp1 = [-1 1];
theta_sp2 = [-5.0 5.0];
theta_sp3 = [-14.99 14.99];
theta0 = 15;
gamma = (abs(theta)/theta0).*(1-log(abs(theta)/theta0)).*(abs(theta)<theta0)+(abs(theta)>theta0);
gamma_sp1 = (abs(theta_sp1)/theta0).*(1-log(abs(theta_sp1)/theta0)).*(abs(theta_sp1)<theta0)+(abs(theta_sp1)>theta0);
gamma_sp2 = (abs(theta_sp2)/theta0).*(1-log(abs(theta_sp2)/theta0)).*(abs(theta_sp2)<theta0)+(abs(theta_sp2)>theta0);
gamma_sp3 = (abs(theta_sp3)/theta0).*(1-log(abs(theta_sp3)/theta0)).*(abs(theta_sp3)<theta0)+(abs(theta_sp3)>theta0);
plot(theta,gamma,'-.','Color','#7E2F8E','LineWidth',2,'HandleVisibility','off');
hold on;
plot(theta_sp1,gamma_sp1,'bs','Color','#008000','LineWidth',2),
plot(theta_sp2,gamma_sp2,'bd','Color','#0000FF','LineWidth',2),
plot(theta_sp3,gamma_sp3,'bo','Color','#FF0000','LineWidth',2),
set(gca,'linewidth',2,'fontsize',20,'fontname','Times');
set(gca,'XTick',(-20:10:20));
set(gca,'YTick',(0.0:0.5:1.1));
set(gca,'yTickLabel',num2str(get(gca,'yTick')','%.1f'))
set(gca,'xminortick','on');
set(gca,'yminortick','on');
xlabel('GB Misorientation Angle  $\Delta\theta_{gg^\prime}$ $(^\circ)$','Fontname', 'Times New Roman','FontSize',20,'Interpreter','latex');
ylabel('$\gamma/\gamma_0$','Fontname', 'Times New Roman','FontSize',20,'Interpreter','latex');
legend('$1^\circ$','$5^\circ$','$15^\circ$','Fontname', 'Times New Roman','FontSize',18,'Interpreter','latex');
axis([-20 20 0 1.1]);
set(gcf,'unit','centimeters','position',[10 10 15 6])

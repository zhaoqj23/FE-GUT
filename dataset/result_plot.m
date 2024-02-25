% clc;
clear;
close all;
tdr = 0.04;
N = 1;
[t,pos,vel] = fload("navdata.txt",'r',[7 Inf],repmat('%f',1,7));
[a,tdd] = tload("tddata.txt",'r',[2 Inf],repmat('%f',1,2));
load('ekf_td.mat');
load("real_traj.mat");
load("ekf_pos.mat");
tddr = tdr*1000 + zeros(length(t),1);



perr = zeros(length(t),2);
perr_ekf = zeros(length(t),2);
perr(:,1) = sqrt(abs(pos(:,2) - p_traj(1,:)').^2+abs(pos(:,1) - p_traj(2,:)').^2);
perr(:,2) = abs(pos(:,3) + p_traj(3,:)');
perr_ekf(:,1) = sqrt(abs(r(:,1) - p_traj(1,:)').^2+abs(r(:,2) - p_traj(2,:)').^2);
perr_ekf(:,2) = abs(r(:,3) - p_traj(3,:)');

hor_ekf = sqrt(sum(perr_ekf(:,1).^2)/length(t));
ver_ekf = sqrt(sum(perr_ekf(:,2).^2)/length(t));
hor = sqrt(sum(perr(:,1).^2)/length(t));
ver = sqrt(sum(perr(:,2).^2)/length(t));
tderror_ekf = sqrt(sum((td-tdr).^2/length(t)))*1000;
tderror = sqrt(sum((tdd-tdr*1000).^2/length(tdd)));
str=['hor=' num2str(hor) '   ver=' num2str(ver) '   tderror=' num2str(tderror)];
disp(str);
str=['hor_ekf=' num2str(hor_ekf) '   ev=' num2str(ver_ekf) '   tderror_ekf=' num2str(tderror_ekf)];
disp(str);
str=['enhance hor=' num2str((hor_ekf-hor)/hor_ekf) '    enhance v=' num2str((ver_ekf-ver)/ver_ekf) '    enhance t=' num2str((tderror_ekf-tderror)/tderror_ekf)];
disp(str);


figure;
set(gcf,'position',[250 300 800 500])
plot(t,td*1000,'LineWidth',1,'Color',[233/255,196/255,107/255]);
% plot(t,td*1000,'LineWidth',0.5,'Color',[130/255,178/255,154/255]);
hold on;
% plot(a,tdd,'linewidth',2,'Color',[33/255,158/255,188/255]);
plot(a,tdd,'linewidth',2,'Color',[75/255,116/255,178/255]);
plot(t,tddr,'LineWidth',2.5,'Color',[193/255,18/255,33/255]);

legend({"EKF","FE-GUT","Nominal Time-Offset"},'Orientation','horizontal');
set(gca,'FontSize',14);
set(gca,'Fontname','times new Roman');
xlabel("Time [s]","FontSize",16);
ylabel("Time-Offset [ms]","FontSize",16);
set(legend,'Location','NorthWest');
set(gca, 'YGrid', 'on');
xlim([0 1200]);
% ylim([-100 200]);


% figure;
% plot(pos(:,2),pos(:,1)); % pos: North-East-Ground
% hold on;
% plot(p_traj(1,:),p_traj(2,:)); % p_traj and r: East-North-Up
% plot(r(:,1),r(:,2));
% legend("FE-GUT","Ground Truth","EKF");
% xlabel("East(m)","FontSize",14);
% ylabel("North(m)","FontSize",14);
% set(gca,'FontSize',14);
% 
% 
% 
% 
% figure;
% subplot(2,1,1);
% plot(t,perr(:,1));
% hold on;
% plot(t,perr_ekf(:,1));
% xlim([1 1200]);
% xlabel("time(s)","FontSize",14);
% ylabel("Horizontal Error(m)","FontSize",14);
% legend('FE-GUT','EKF');
% set(gca,'FontSize',14);
% set(gca, 'YGrid', 'on');



% subplot(2,1,2);
% plot(t,perr(:,2));
% hold on;
% plot(t,perr_ekf(:,2));
% xlim([1 1200]);
% xlabel("time(s)","FontSize",14);
% ylabel("Vertical Error(m)","FontSize",14);
% set(gca,'FontSize',14);
% set(gca, 'YGrid', 'on');
% legend('FE-GUT','EKF');



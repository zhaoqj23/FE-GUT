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
figure;
plot(a(N:end),tdd(N:end),'linewidth',3);
hold on;
plot(t(:),td(:)*1000);
plot(t(:),tddr,'linewidth',3);
legend("This Work","EKF","Nominal t_d");
xlabel("time(s)");
ylabel("td(ms)");
xlim([0 1200])


figure;
plot(pos(:,2),pos(:,1)); % pos: North-East-Ground
hold on;
plot(p_traj(1,:),p_traj(2,:)); % p_traj and r: East-North-Up
plot(r(:,1),r(:,2));
legend("This Work","Ground Truth","EKF");
xlabel("East(m)");
ylabel("North(m)");

perr = zeros(length(t),2);
perr_ekf = zeros(length(t),2);
perr(:,1) = sqrt(abs(pos(:,2) - p_traj(1,:)').^2+abs(pos(:,1) - p_traj(2,:)').^2);
perr(:,2) = abs(pos(:,3) + p_traj(3,:)');
perr_ekf(:,1) = sqrt(abs(r(:,1) - p_traj(1,:)').^2+abs(r(:,2) - p_traj(2,:)').^2);
perr_ekf(:,2) = abs(r(:,3) - p_traj(3,:)');
figure;
subplot(2,1,1);
plot(t,perr(:,1));
hold on;
plot(t,perr_ekf(:,1));
xlim([1 1200]);
xlabel("time(s)");
ylabel("Horizontal Error(m)");
legend('This Work','EKF');
subplot(2,1,2);
plot(t,perr(:,2));
hold on;
plot(t,perr_ekf(:,2));
xlim([1 1200]);
xlabel("time(s)");
ylabel("Vertical Error(m)");
legend('This Work','EKF');

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
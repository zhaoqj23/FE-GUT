clc;
clear;
close all;
load('Data\ekf_P.mat');
load('Data\ekf_X.mat');
Fs = 10;
dt = 1/Fs;
epoch_num = 200;
shiftnum = 100;
x_num_each_epoch = 12;

%数据加载
load("Data\GNSS_Simulation.mat");
load("Data\traj.mat");
p_traj = p_traj';v_traj = v_traj';

Sj = [0.4 0.4 0.4];St = 36;Sf = 0.01;Sdt = 0.01;%加速度、钟差、钟漂、时偏的功率谱密度
std_p = 2; std_pr = 0.1; std_uwb = 0.1;

recPos0 = [39.904987 116.405289 60.0352]; %初始位置
UWB_base_enu = [150,0,5;100,50,5;50,0,5;100,-50,5];
UWB_base_ecef = zeros(length(UWB_base_enu(:,1)),3);
%计算UWB基站在ECEF坐标系下的坐标
wgs84 = wgs84Ellipsoid('meter');
for i = 1:length(UWB_base_enu(:,1))
    [UWB_base_ecef(i,1),UWB_base_ecef(i,2),UWB_base_ecef(i,3)] = enu2ecef(UWB_base_enu(i,1),UWB_base_enu(i,2),UWB_base_enu(i,3),recPos0(1),recPos0(2),recPos0(3),wgs84);
end

TXk = Xk_save(1+shiftnum:(epoch_num+shiftnum),2:13);
PXk = Pk_save(1+shiftnum:(epoch_num+shiftnum),2:145);
N = size(TXk,1);

%状态递推噪声协方差矩阵Q
Q = zeros(x_num_each_epoch,x_num_each_epoch);
Q(1:3,1:3) = 1/20*diag(Sj)*dt^5;
Q(4:6,4:6) = 1/3*diag(Sj)*dt^3;
Q(7:9,7:9) = diag(Sj)*dt;
Q(4:6,1:3) = 1/8*diag(Sj)*dt^4; Q(1:3,4:6) = 1/8*diag(Sj)*dt^4;
Q(7:9,4:6) = 1/2*diag(Sj)*dt^2; Q(4:6,7:9) = 1/2*diag(Sj)*dt^2;
Q(7:9,1:3) = 1/6*diag(Sj)*dt^3; Q(1:3,7:9) = 1/6*diag(Sj)*dt^3;
Q(10:11,10:11) = [St*dt+1/3*Sf*dt^3 1/2*Sf*dt^2;1/2*Sf*dt^2 Sf*dt];
Q(12,12) = Sdt*dt;
Qc = chol(Q);

n = 6; %卫星数量
m = 4; %UWB锚点数量
%量测噪声协方差矩阵R
R = eye(2*n+m);
R(1:n,1:n) = eye(n)*(std_p)^2;
R(n+1:2*n,n+1:2*n) = eye(n)*(std_pr)^2;
R(2*n+1:2*n+m,2*n+1:2*n+m) = eye(m)*(std_uwb)^2;
Rg = eye(2*n);
R(1:n,1:n) = eye(n)*(std_p)^2;
R(n+1:2*n,n+1:2*n) = eye(n)*(std_pr)^2;
Rc = chol(R);

%状态递推矩阵
Fk = eye(x_num_each_epoch);
Fk(1:3,4:6) = eye(3)*dt;
Fk(1:3,7:9) = 1/2*eye(3)*dt^2;
Fk(4:6,7:9) = eye(3)*dt;
Fk(10,11) = dt;

%% 生成随机数
N_LON = 1000; % 蒙特卡洛次数


EXk = zeros(N,x_num_each_epoch,N_LON);
EXkm = zeros(N,x_num_each_epoch,N_LON);
EXk1 = EXk;
z = zeros(N,2*n+m,N_LON);
Hm = zeros(N,2*n+m,x_num_each_epoch,N_LON);

for i = 1:N
    Xx = TXk(i,:);
    Px = reshape(PXk(i,:), 12, 12);
    Pxc = chol(Px);
    EXk(i,:,:) = repmat(Xx',1,N_LON) + Pxc*randn(x_num_each_epoch,N_LON);
    Xx(12) = 0.04;
    EXkm(i,:,:) = repmat(Xx',1,N_LON) + Pxc*randn(x_num_each_epoch,N_LON);
end

for i = 1:N
    for j = 1:N_LON
        % 计算卫星测距信息
        ri = EXk(i,1:3,j);
        vi = EXk(i,4:6,j);
        ai = EXk(i,7:9,j);
        satpos = Satposition{i+shiftnum,2};
        tuk = EXk(i,10,j);
        fuk = EXk(i,11,j);
        tdk = EXk(i,12,j);
        cgnss = repmat(ri,n,1);%位置向量纵向复制n次
        rgnss = sqrt(sum((satpos - cgnss).^2, 2));%与卫星的距离向量
        HGn = -(satpos-cgnss)./repmat(rgnss,1,3);%距离向量横向复制3次，计算得到距离相对载体位置的Jacobian矩阵
        pd = rgnss + tuk;
        % 计算多普勒频移测量信息
        cv = repmat(vi,length(satpos(:,1)),1);
        prd = HGn*vi' + fuk;
        % 计算UWB测距信息
        dr = ri - vi*tdk - 0.5*ai*tdk^2;%矫正时偏之后的
        cuwb = repmat(dr,length(UWB_base_ecef(:,1)),1);%矫正后的位置向量纵向复制n次
        ruwb = sqrt(sum((UWB_base_ecef - cuwb).^2, 2));%与4个UWB基站的距离向量   ?
        HUm = -(UWB_base_ecef - cuwb)./repmat(ruwb,1,3);
        yd = [pd;prd;ruwb];
        z(i,:,j) = yd + Rc*randn(2*n+m,1);
        % 计算Jacobian矩阵
        tem = vi+ai*tdk;
        Htdk = -HUm*tem';
        H = zeros(2*n+m,x_num_each_epoch);
        H(1:n,1:3) = HGn;
        H(1:n,10) = ones(n,1);
        H((n+1):2*n,4:6) = HGn;
        H((n+1):2*n,11) = ones(n,1);
        H((2*n+1):(2*n+m),1:3) = HUm;
        H((2*n+1):(2*n+m),4:6) = -HUm*tdk;
        H((2*n+1):(2*n+m),7:9) = -0.5*HUm*tdk^2;
        H((2*n+1):(2*n+m),12) = Htdk;
        Hm(i,:,:,j) = H;
        if i > 1
            xk1 = Fk*(EXk(i-1,:,j))';
            EXk1(i,:,j) = xk1' + (Qc*randn(x_num_each_epoch,1))';
        end
    end
end
%% 计算离散时间状态空间模型非线性度
M = zeros(N-1,1);
for i = 1:N-1
    Cov_yy = zeros(x_num_each_epoch+2*n+m,x_num_each_epoch+2*n+m);
    Cov_xx = zeros(x_num_each_epoch,x_num_each_epoch);
    Cov_yx = zeros(x_num_each_epoch+2*n+m,x_num_each_epoch);
    Cov_xy = zeros(x_num_each_epoch,x_num_each_epoch+2*n+m);

    mean_x = zeros(x_num_each_epoch,1);
    mean_y = zeros(x_num_each_epoch+2*n+m,1);
    for j = 1:x_num_each_epoch
        mean_x(j) = mean(EXk(i,j,:));
    end

    for k = 1 : (x_num_each_epoch+2*n+m)
        if k <= x_num_each_epoch
            mean_y(k) = mean(EXk1(i+1,k,:));
        else
            mean_y(k) = mean(z(i,k-x_num_each_epoch,:));
        end
    end

    for l = 1 : N_LON
        x = EXk(i,:,l)';
        y = [EXk1(i+1,:,l)';z(i,:,l)'];

        Cov_xx = Cov_xx + (x-mean_x)*(x-mean_x)';
        Cov_yy = Cov_yy + (y-mean_y)*(y-mean_y)';
        B = squeeze(Hm(i,:,:,l));
        A = [Fk;B];
        Cov_yx = Cov_yx + A*(x-mean_x)*(x-mean_x)';
        Cov_xy = Cov_xy + (x-mean_x)*(A*(x-mean_x))';
    end

    Cov_xx = Cov_xx/(N_LON-1);
    Cov_yy = Cov_yy/(N_LON-1);
    Cov_yx = Cov_yx/(N_LON-1);
    Cov_xy = Cov_xy/(N_LON-1);

    M(i) = 1-trace(inv(Cov_yy)*Cov_yx*inv(Cov_xx)*Cov_xy)/(2*n+m+x_num_each_epoch);
end

%% 计算图状态空间模型非线性度
winsize=[5,18,30];
M_opt = zeros(3,N-1,1);
for iwinsize=1:3
    win_size=winsize(iwinsize);
    mean_x = [];
    mean_y = [];
    x = [];
    y = [];
    ynum = 2*n+m+x_num_each_epoch-1;
    xnum = x_num_each_epoch - 1;
    H_win = zeros(win_size*ynum+1,win_size*xnum+1,N_LON); 
    z_opt = zeros(N,2*n+m,N_LON);

    EXk_opt = EXk;
    td0_randn = EXkm(1,12,:);
    EXk_opt(:,12,:) = repmat(td0_randn,[epoch_num,1,1]);

    H_opt = zeros(N,2*n+m,x_num_each_epoch,N_LON);
    EXk1_opt = zeros(size(EXk));
    for i = 1 : N
        for j = 1 : N_LON
        % 计算卫星测距信息
        ri = EXk_opt(i,1:3,j);
        vi = EXk_opt(i,4:6,j);
        ai = EXk_opt(i,7:9,j);
        satpos = Satposition{i+shiftnum,2};
        tuk = EXk_opt(i,10,j);
        fuk = EXk_opt(i,11,j);
        tdk = EXk_opt(i,12,j);
        cgnss = repmat(ri,n,1);%位置向量纵向复制n次
        rgnss = sqrt(sum((satpos - cgnss).^2, 2));%与卫星的距离向量
        HGn = -(satpos-cgnss)./repmat(rgnss,1,3);%距离向量横向复制3次，计算得到距离相对载体位置的Jacobian矩阵
        pd = rgnss + tuk;
        % 计算多普勒频移测量信息
        cv = repmat(vi,length(satpos(:,1)),1);
        prd = HGn*vi' + fuk;
        % 计算UWB测距信息
        dr = ri - vi*tdk - 0.5*ai*tdk^2;%矫正时偏之后的
        cuwb = repmat(dr,length(UWB_base_ecef(:,1)),1);%矫正后的位置向量纵向复制n次
        ruwb = sqrt(sum((UWB_base_ecef - cuwb).^2, 2));%与4个UWB基站的距离向量   ?
        HUm = -(UWB_base_ecef - cuwb)./repmat(ruwb,1,3);
        yd = [pd;prd;ruwb];
        z_opt(i,:,j) = yd + Rc*randn(2*n+m,1);
        % 计算Jacobian矩阵
        tem = vi+ai*tdk;
        Htdk = -HUm*tem';
        H = zeros(2*n+m,x_num_each_epoch);
        H(1:n,1:3) = HGn;
        H(1:n,10) = ones(n,1);
        H((n+1):2*n,4:6) = HGn;
        H((n+1):2*n,11) = ones(n,1);
        H((2*n+1):(2*n+m),1:3) = HUm;
        H((2*n+1):(2*n+m),4:6) = -HUm*tdk;
        H((2*n+1):(2*n+m),7:9) = -0.5*HUm*tdk^2;
        H((2*n+1):(2*n+m),12) = Htdk;
        H_opt(i,:,:,j) = H;
            if i > 1
                xk1_opt = Fk*(EXk_opt(i-1,:,j))';
                EXk1_opt(i,:,j) = xk1_opt' + (Qc*randn(x_num_each_epoch,1))'; 
                % EXk1_opt(i,:,j) = xk1_opt';
            end
        end
    end
    mean_x(1,1) = mean(td0_randn);
    x(1,1:N_LON) = td0_randn;
    H_win(1,1,1:N_LON) = ones(1,N_LON);

    for k = 1:N-1
        td0_randn_k1 = squeeze(td0_randn)' + sqrt(Q(end,end))*randn(1,N_LON);
        % td0_randn_k1 = squeeze(td0_randn)';
        mean_y(1,1) = mean(td0_randn_k1);
        y(1,1:N_LON) = td0_randn_k1;

        Cov_yy = zeros(win_size*ynum+1,win_size*ynum+1);
        Cov_xx = zeros(win_size*xnum+1,win_size*xnum+1);
        Cov_yx = zeros(win_size*ynum+1,win_size*xnum+1);
        Cov_xy = zeros(win_size*xnum+1,win_size*ynum+1);

        if k > win_size
            mean_x(2:(xnum+1),:) = [];
            mean_y(2:(ynum+1),:) = [];
            x(2:(xnum+1),:) = [];
            y(2:(ynum+1),:) = [];
            for i = 1:xnum
                mean_x((win_size-1)*xnum+1+i,:) = mean(EXk_opt(k,i,:));
                x((win_size-1)*xnum+1+i,1:N_LON) = EXk_opt(k,i,:);
            end

            for i = 1:ynum
                if i <= xnum
                    mean_y((win_size-1)*ynum+1+i,:) = mean(EXk1_opt(k+1,i,:));
                    y((win_size-1)*ynum+1+i,1:N_LON) = EXk1_opt(k+1,i,:);
                else
                    mean_y((win_size-1)*ynum+1+i,:) = mean(z_opt(k,i-xnum,:));
                    y((win_size-1)*ynum+1+i,1:N_LON) = z_opt(k,i-xnum,:);
                end
            end
        
            H_win(2:(ynum+1),:,:) = [];
            H_win(:,2:(xnum+1),:) = [];

            for i = 1:N_LON
                H_win((win_size-1)*ynum+2:(win_size-1)*ynum+1+xnum,(win_size-1)*xnum+2:win_size*xnum+1,i) = Fk(1:11,1:11);
                H_win((win_size-1)*ynum+2+xnum:win_size*ynum+1,1,i) = H_opt(k,1:16,12,i);
                H_win((win_size-1)*ynum+2+xnum:win_size*ynum+1,(win_size-1)*xnum+2:win_size*xnum+1,i) = H_opt(k,1:16,1:11,i);
            end
        else
            for i = 1:xnum
                mean_x((k-1)*xnum+1+i,:) = mean(EXk_opt(k,i,:));
                x((k-1)*xnum+1+i,1:N_LON) = EXk_opt(k,i,:);
            end
            for i = 1:ynum
                if i <= xnum
                    mean_y((k-1)*ynum+1+i,:) = mean(EXk1_opt(k+1,i,:));
                    y((k-1)*ynum+1+i,1:N_LON) = EXk1_opt(k+1,i,:);
                else
                    mean_y((k-1)*ynum+1+i,:) = mean(z_opt(k,i-xnum,:));
                    y((k-1)*ynum+1+i,1:N_LON) = z_opt(k,i-xnum,:);
                end
            end
            for i = 1:N_LON
                H_win((k-1)*ynum+2:(k-1)*ynum+1+xnum,(k-1)*xnum+2:k*xnum+1,i) = Fk(1:11,1:11);
                H_win((k-1)*ynum+2+xnum:k*ynum+1,1,i) = H_opt(k,1:16,12,i);
                H_win((k-1)*ynum+2+xnum:k*ynum+1,(k-1)*xnum+2:k*xnum+1,i) = H_opt(k,1:16,1:11,i);
            end
        end
        if k >= win_size
            for i = 1:N_LON
                Cov_xx=Cov_xx+(x(:,i)-mean_x)*(x(:,i)-mean_x)';       
                Cov_yy=Cov_yy+(y(:,i)-mean_y)*(y(:,i)-mean_y)';        
                Cov_yx=Cov_yx+H_win(:,:,i)*(x(:,i)-mean_x)*(x(:,i)-mean_x)';
                Cov_xy=Cov_xy+(x(:,i)-mean_x)*(H_win(:,:,i)*(x(:,i)-mean_x))';
            end
            Cov_xx=Cov_xx/(N_LON-1);
            Cov_yy=Cov_yy/(N_LON-1);
            Cov_yx=Cov_yx/(N_LON-1);
            Cov_xy = Cov_xy/(N_LON-1);
            M_opt(iwinsize,k) = 1-trace(inv(Cov_yy)*Cov_yx*inv(Cov_xx)*Cov_xy)/(win_size*ynum+1); %#ok<*MINV>
        end
    end
end

time = (1:(N-1))*dt;
figure;
plot(time(1:N-1),M,'r-','LineWidth',2);%,M_wekf(win_size:N-1),time(win_size:N-1),M_opt(win_size:N-1));
hold;
% 逐个添加图例项
labelVarwin5 = 'GSSM with window size ' + string(5); % 假设value是一个你想在图例中显示的变量
labelVarwin18 = 'GSSM with window size ' + string(18); % 假设value是一个你想在图例中显示的变量
%labelVarwin20 = 'GSSM with  window size' + string(20)+'  window size'; % 假设value是一个你想在图例中显示的变量
labelVarwin30 = 'GSSM with window size ' + string(30); % 假设value是一个你想在图例中显示的变量
plot(time(winsize(1):N-1),M_opt(1,winsize(1):N-1,1),'b-','LineWidth',2);
plot(time(winsize(2):N-1),M_opt(2,winsize(2):N-1,1),'g-','LineWidth',2);
plot(time(winsize(3):N-1),M_opt(3,winsize(3):N-1,1),'k-','LineWidth',2);
hold;
set(gca, 'FontSize', 24);
ylabel('Nonlinearity','fontsize',24);xlabel('Time(s)','fontsize',24);
%box off;
grid on;
legend('TDTSSM',labelVarwin5,labelVarwin18,labelVarwin30,'FontSize', 24);%
% 获取当前图例的句柄
lgd = legend;
% 获取图例的当前位置和大小
pos = get(lgd, 'Position');
% 设置图例框的大小（根据需要调整数值）
pos(3) = pos(3) * 1.5; % 增加宽度为原来的1.5倍
pos(4) = pos(4) * 1.5; % 增加高度为原来的1.5倍
% 应用新的位置和大小
set(lgd, 'Position', pos);





















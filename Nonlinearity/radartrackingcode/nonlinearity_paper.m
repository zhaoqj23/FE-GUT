clc;
clear;
close all;

dt=0.05; % Sample time
epoch_num=400; %Epoch number
x_num_each_epoch=3;


filepath='Data\radardata2.txt';% data Storage
data_kf=dlmread([ filepath ]);

TXk=data_kf(:,2:3);
TXk(:,3)=ones(size(TXk,1),1).*1000;

N=size(TXk,1);%n epoch
%% Random number generate

N_LON=1000;%Monte Carlo time

sigma_X=10;
sigma_Z=2;
sigma_predict=5;


EXk=zeros(N,3,N_LON);
EXk1=EXk;
z=zeros(N,N_LON);
H=zeros(N,3,N_LON);

for i=1:N

    x(1,1,1:N_LON)=TXk(i,1).*ones(1,N_LON)+sigma_X.*randn(1,N_LON);
    v(1,1,1:N_LON)=TXk(i,2).*ones(1,N_LON)+sigma_X.*randn(1,N_LON);
    h(1,1,1:N_LON)=TXk(i,3).*ones(1,N_LON)+sigma_X.*randn(1,N_LON);

    EXk(i,1,:)=x;
    EXk(i,2,:)=v;
    EXk(i,3,:)=h;

end

Fk=[1 dt 0;
    0 1 0;
    0 0 1];

for k=1:N
    temp(1,1:N_LON)=sqrt(EXk(k,1,:).^2+EXk(k,3,:).^2);

    z(k,:)=temp+sigma_Z.*randn(1,N_LON);%

    xk(1,1:N_LON)=EXk(k,1,:);
    hk(1,1:N_LON)=EXk(k,3,:);
    H(k,:,:)=[xk./temp ;zeros(1,N_LON);hk./temp ];
%     Hx(k,:)=
%     for i=1:N_LON
%         Hx(k,i)=H(1:3,i)'*[xk(i);0 ;hk(i)]-temp(i);
%     end

    if k>1
        for i=1:N_LON
    %         xk()=EXk(k-1,:,i)
            xk1=Fk*EXk(k-1,:,i)';

            EXk1(k,:,i)=xk1+sigma_predict.*randn(3,1);
        end
    end
end


%% Non-linearity test
%win_size=30;
winsize=[5,18,30];

%% Graphical State Space Model------------------------------
M_opt=zeros(3,N-1,1);
for iwinsize=1:3
    win_size=winsize(iwinsize);
    mean_x=[];
    mean_y=[];
    x=[];
    y=[];
    a=[];
    H_win=zeros(2*win_size+2,win_size+2,i);
    z_opt=zeros(N,N_LON);

    v0_randn=TXk(1,2).*ones(1,N_LON)+sigma_X.*randn(1,N_LON);
    h0_randn=TXk(1,3).*ones(1,N_LON)+sigma_X.*randn(1,N_LON);

    EXk_opt=EXk;
    EXk_opt(:,2,:)=repmat(v0_randn,[400,1,1]);
    EXk_opt(:,3,:)=repmat(h0_randn,[400,1,1]);

    H_opt=zeros(N,2,N_LON);
    EXk1_opt=zeros(N,3,N_LON);
    for k=1:N
        temp(1,1:N_LON)=sqrt(EXk_opt(k,1,:).^2+EXk_opt(k,3,:).^2);
        z_opt(k,:)=temp+sigma_Z.*randn(1,N_LON);%

        xk(1,1:N_LON)=EXk(k,1,:);
        hk(1,1:N_LON)=EXk(k,3,:);

        H_opt(k,:,:)=[xk./temp ;hk./temp ];

        if k>1
            for i=1:N_LON
                xk1_opt=Fk*EXk_opt(k-1,:,i)';
                EXk1_opt(k,:,i)=xk1_opt+sigma_predict.*randn(3,1);
            end
        end
    end

    mean_x(1,1)=mean(v0_randn);
    mean_x(2,1)=mean(h0_randn);

    x(1,1:N_LON)=v0_randn;
    x(2,1:N_LON)=h0_randn;

    H_win(1,1,:)=ones(1,N_LON);
    H_win(2,2,:)=ones(1,N_LON);
    for k=1:N-1
        v0_randn_k1=v0_randn+sigma_predict.*randn(1,N_LON);
        h0_randn_k1=h0_randn+sigma_predict.*randn(1,N_LON);
        mean_y(1,1)=mean(v0_randn_k1);
        mean_y(2,1)=mean(h0_randn_k1);
        y(1,1:N_LON)=v0_randn_k1;
        y(2,1:N_LON)=h0_randn_k1;

        Cov_yy=zeros(2*win_size+2,2*win_size+2);
        Cov_xx=zeros(win_size+2,win_size+2);
        Cov_yx=zeros(2*win_size+2,win_size+2);

        if k>win_size
            mean_x(3,:)=[];%x=[v h x1 x2 x3 …… x_winsize]
            mean_y(3:4,:)=[];%y=[v h x1 z1 x2 z2 x3 z3 …… x_winsize z_winsize]
            mean_x(win_size+2,1)=mean(EXk_opt(k,1,:),3)';% Evaluation for  the last epoch
            mean_y(win_size*2+1:win_size*2+2,1)=[mean(EXk1_opt(k+1,1,:),3)';mean(z_opt(k,:),2)];% Evaluation for x, z in the last epoch

            x(3,:)=[];
            y(3:4,:)=[];
            x(win_size+2,1:N_LON)=EXk_opt(k,1,:);
            a(1,1:N_LON)=EXk1_opt(k+1,1,:);
            y((win_size-1)*2+1+2:win_size*2+2,1:N_LON)=[a;z_opt(k,:)];


            H_win(3:4,:,:)=[];
            H_win(:,3,:)=[];
            for i=1:N_LON
                H_win(win_size*2+1,1,i)=dt;
                H_win(win_size*2+2,2,i)=H_opt(k,2,i);
                H_win(win_size*2+1:win_size*2+2,win_size+2,i)=[1;H_opt(k,1,i)];
            end

        else
            mean_x(k+2,1)=mean(EXk_opt(k,1,:),3)';
            mean_y(2*k+1:2*k+2,1)=[mean(EXk1_opt(k+1,1,:),3)';mean(z_opt(k,:),2)];

            x(k+2,1:N_LON)=EXk_opt(k,1,:);
            a(1,1:N_LON)=EXk1_opt(k+1,1,:);
            y(2*k+1:2*k+2,1:N_LON)=[a;z_opt(k,:)];

            for i=1:N_LON
                H_win(k*2+1,1,i)=dt;
                H_win(k*2+2,2,i)=H_opt(k,2,i);
                H_win(k*2+1:k*2+2,k+2,i)=[1;H_opt(k,1,i)];
            end
        end

        if k>=win_size
            for i=1:N_LON
                Cov_xx=Cov_xx+(x(:,i)-mean_x)*(x(:,i)-mean_x)';
                Cov_yy=Cov_yy+(y(:,i)-mean_y)*(y(:,i)-mean_y)';
%             Cov_yx=Cov_yx+(y(:,i)-mean_y)*(x(:,i)-mean_x)';
                Cov_yx=Cov_yx+H_win(:,:,i)*(x(:,i)-mean_x)*(x(:,i)-mean_x)';
            end
            Cov_xx=Cov_xx/N_LON;
            Cov_yy=Cov_yy/N_LON;
            Cov_yx=Cov_yx/N_LON;

            M_opt(iwinsize,k)=sqrt(trace(Cov_yy-Cov_yx*inv(Cov_xx)*Cov_yx'))/sqrt(trace(Cov_yy));
        end
    end
end

%% Graphical State Space Model------------------------------

%% Caluation for single epoch Extended Kalman Filter------------------------------
M=zeros(N-1,1);
for k=1:N-1
    Cov_yy=zeros(4,4);
    Cov_xx=zeros(3,3);

    Cov_yx=zeros(4,3);

    mean_x=mean(EXk(k,:,:),3)';
     mean_y=[mean(EXk1(k+1,:,:),3)';mean(z(k,:),2)];
%    mean_y=[mean(EXk(k+1,:,:),3)';mean(z(k,:),2)];
    for i=1:N_LON
        x=EXk(k,:,i)';
         y=[EXk1(k+1,:,i)';z(k,i)];
%        y=[EXk(k+1,:,i)';z(k,i)];
        Cov_xx=Cov_xx+(x-mean_x)*(x-mean_x)';

        Cov_yy=Cov_yy+(y-mean_y)*(y-mean_y)';

%         Cov_yx=Cov_yx+(y-mean_y)*(x-mean_x)';
        Cov_yx=Cov_yx+[Fk;H(k,:,i)]*(x-mean_x)*(x-mean_x)';

    end

    Cov_xx=Cov_xx/N_LON;
    Cov_yy=Cov_yy/N_LON;
    Cov_yx=Cov_yx/N_LON;

    M(k)=sqrt(trace(Cov_yy-Cov_yx/Cov_xx*Cov_yx'))/sqrt(trace(Cov_yy));
end


figure(1);
%plot(data_kf(1:N-1,1),M,data_kf(win_size:N-1,1),M_kf(win_size:N-1,1),...
  %  data_kf(win_size:N-1,1),M_opt(win_size:N-1,1));
plot(data_kf(1:N-1,1),M, 'r-','LineWidth',2); 
hold;
% 逐个添加图例项
labelVarwin5 = 'GSSM with window size ' + string(5); % 假设value是一个你想在图例中显示的变量
labelVarwin18 = 'GSSM with window size ' + string(18); % 假设value是一个你想在图例中显示的变量
%labelVarwin20 = 'GSSM with  window size' + string(20)+'  window size'; % 假设value是一个你想在图例中显示的变量
labelVarwin30 = 'GSSM with window size ' + string(30); % 假设value是一个你想在图例中显示的变量
plot(data_kf(winsize(1):N-1,1),M_opt(1,winsize(1):N-1,1), 'b-','LineWidth',2);
plot(data_kf(winsize(2):N-1,1),M_opt(2,winsize(2):N-1,1), 'g-','LineWidth',2);
plot(data_kf(winsize(3):N-1,1),M_opt(3,winsize(3):N-1,1), 'k-','LineWidth',2);
%plot(data_kf(winsize(4):N-1,1),M_opt(4,winsize(4):N-1,1), 'k-');
hold;
% 设置X轴和Y轴的刻度标签字体大小
set(gca, 'FontSize', 24);
ylabel('Nonlinearity','fontsize',24);xlabel('Time(s)','fontsize',24);
box off;
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




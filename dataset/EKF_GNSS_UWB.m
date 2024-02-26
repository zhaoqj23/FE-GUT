% clc;
clear;
close all;

load("GNSS_Simulation.mat");
load("UWB_Range.mat");

Fs_g = 10; Fs_u = 200; ddt = 1/Fs_u; T = 1200;
time = 0:1/Fs_g:T;
utime = 0:1/Fs_u:T;
vini = 0;
recPos0 = [39.904987 116.405289 60.0352]; % Initial position
UWB_base_enu = [150,0,5;100,50,5;50,0,5;100,-50,5]; % UWB anchors
UWB_base_ecef = zeros(length(UWB_base_enu(:,1)),3);

wgs84 = wgs84Ellipsoid('meter');
for i = 1:length(UWB_base_enu(:,1))
    [UWB_base_ecef(i,1),UWB_base_ecef(i,2),UWB_base_ecef(i,3)] = enu2ecef(UWB_base_enu(i,1),UWB_base_enu(i,2),UWB_base_enu(i,3),recPos0(1),recPos0(2),recPos0(3),wgs84);
end

r = zeros(length(time),3);
v = zeros(length(time),3);
a = zeros(length(time),3);
speed = zeros(length(time),1);
td = zeros(length(time),1);
tu = zeros(length(time),1);
fu = zeros(length(time),1);


n = 6; % The number of satellites
m = 4; % The number of UWB anchors
ru0 = zeros(3,1);

au0 = zeros(3,1);
tu0 = 0;
fu0 = 0;
td0 = 0;



[ru0(1),ru0(2),ru0(3)] = geodetic2ecef(wgs84,recPos0(1),recPos0(2),recPos0(3));

R_ecef2enu = [-sin(recPos0(2)*pi/180),              -sin(recPos0(1)*pi/180)*cos(recPos0(2)*pi/180),    cos(recPos0(1)*pi/180)*cos(recPos0(2)*pi/180);
              cos(recPos0(2)*pi/180),     -sin(recPos0(1)*pi/180)*sin(recPos0(2)*pi/180),              cos(recPos0(1)*pi/180)*sin(recPos0(2)*pi/180);
              0,      cos(recPos0(1)*pi/180),      sin(recPos0(1)*pi/180)];

vu0 = R_ecef2enu*[0;vini;0];
X0 = [ru0;vu0;au0;tu0;fu0;td0]; 
P0 = 0.1*eye(length(X0));
Sj = [0.4 0.4 0.4]; St = 36;Sf = 0.01;Sdt = 0.01;
std_p = 0.4; std_pr = 0.05; std_uwb = 0.5;


% The Q matrix
Q = zeros(length(X0),length(X0));
Q(1:3,1:3) = 1/20*diag(Sj)*ddt^5;
Q(4:6,4:6) = 1/3*diag(Sj)*ddt^3;
Q(7:9,7:9) = diag(Sj)*ddt;
Q(4:6,1:3) = 1/8*diag(Sj)*ddt^4; Q(1:3,4:6) = 1/8*diag(Sj)*ddt^4;
Q(7:9,4:6) = 1/2*diag(Sj)*ddt^2; Q(4:6,7:9) = 1/2*diag(Sj)*ddt^2;
Q(7:9,1:3) = 1/6*diag(Sj)*ddt^3; Q(1:3,7:9) = 1/6*diag(Sj)*ddt^3;
Q(10:11,10:11) = [St*ddt+1/3*Sf*ddt^3 1/2*Sf*ddt^2;1/2*Sf*ddt^2 Sf*ddt];
Q(12,12) = Sdt*ddt;
% The R matrix
R = eye(2*n+m);
R(1:n,1:n) = eye(n)*(std_p)^2;
R(n+1:2*n,n+1:2*n) = eye(n)*(std_pr)^2;
R(2*n+1:2*n+m,2*n+1:2*n+m) = eye(m)*(std_uwb)^2;
Rg = eye(2*n);
R(1:n,1:n) = eye(n)*(std_p)^2;
R(n+1:2*n,n+1:2*n) = eye(n)*(std_pr)^2;
Ruwb = eye(m)*(std_uwb)^2;

Fk = eye(length(X0));
Fk(1:3,4:6) = eye(3)*ddt;
Fk(1:3,7:9) = 1/2*eye(3)*ddt^2;
Fk(4:6,7:9) = eye(3)*ddt;
Fk(10,11) = ddt;
X_old = X0;
P_old = P0;
index_g = 2;


for i = 2:length(utime)
    if utime(i) < time(index_g)
        % State prediction
        Xk_pri = Fk*X_old;
        Pk_pri = Fk*P_old*Fk' + Q;
        uwb = UWB_Range(i,2:5)';
        yk = uwb;
        % Calculate the Jacobian
        ri = Xk_pri(1:3)';
        vi = Xk_pri(4:6)';
        ai = Xk_pri(7:9)';
        tdk = Xk_pri(end);
        dr = ri - vi*tdk - 0.5*ai*tdk^2;
        cuwb = repmat(dr,length(UWB_base_ecef(:,1)),1);
        ruwb = sqrt(sum((UWB_base_ecef - cuwb).^2, 2));
        HUm = -(UWB_base_ecef - cuwb)./repmat(ruwb,1,3);
        tem = vi+ai*tdk;
        Htdk = -HUm*tem';
        H = zeros(m,length(X0));
        H(1:m,1:3) = HUm;
        H(1:m,4:6) = -HUm*tdk;
        H(1:m,7:9) = -0.5*HUm*tdk^2;
        H(1:m,12) = Htdk;
        % Calculate the Kalman gain
        A = Pk_pri*H';
        B = (H*Pk_pri*H'+Ruwb);
        Kk = A/B;
        yd = ruwb;
        % update X and P
        Xk_post = Xk_pri+Kk*(yk-yd);
        Pk_post = (eye(length(X0))-Kk*H)*Pk_pri;
        X_old = Xk_post;
        P_old = Pk_post;
    else
        % State prediction
        Xk_pri = Fk*X_old;
        Pk_pri = Fk*P_old*Fk' + Q;
        pseudo = ps{index_g,2};
        pseudo_rate = ps_rate{index_g,2};
        satpos = Satposition{index_g,2};
        uwb = UWB_Range(i,2:5)';
        yk = [pseudo;pseudo_rate;uwb];
        % Calculate the Jacobian
        ri = Xk_pri(1:3)';
        cgnss = repmat(ri,length(satpos(:,1)),1);
        rgnss = sqrt(sum((satpos - cgnss).^2, 2));
        HGn = -(satpos-cgnss)./repmat(rgnss,1,3);
        vi = Xk_pri(4:6)';
        ai = Xk_pri(7:9)';
        tdk = Xk_pri(end);
        dr = ri - vi*tdk - 0.5*ai*tdk^2;
        cuwb = repmat(dr,length(UWB_base_ecef(:,1)),1);
        ruwb = sqrt(sum((UWB_base_ecef - cuwb).^2, 2));
        HUm = -(UWB_base_ecef - cuwb)./repmat(ruwb,1,3);
        tem = vi+ai*tdk;
        Htdk = -HUm*tem';
        H = zeros(2*n+m,length(X0));
        H(1:n,1:3) = HGn;
        H(1:n,10) = ones(n,1);
        H((n+1):2*n,4:6) = HGn;
        H((n+1):2*n,11) = ones(n,1);
        H((2*n+1):(2*n+m),1:3) = HUm;
        H((2*n+1):(2*n+m),4:6) = -HUm*tdk;
        H((2*n+1):(2*n+m),7:9) = -0.5*HUm*tdk^2;
        H((2*n+1):(2*n+m),12) = Htdk;
        % Calculate the Kalman gain
        A = Pk_pri*H';
        B = (H*Pk_pri*H'+R);
        Kk = A/B;
        tuk = Xk_pri(10);
        fuk = Xk_pri(11);
        pd = rgnss + tuk;
        cv = repmat(vi,length(satpos(:,1)),1);
        trans_rece = (satpos-cgnss)./rgnss;
        prd = HGn*vi' + fuk;
        rd = ruwb;
        yd = [pd;prd;rd];
        % update X and P
        Xk_post = Xk_pri+Kk*(yk-yd);
        Pk_post = (eye(length(X0))-Kk*H)*Pk_pri;
        X_old = Xk_post;
        P_old = Pk_post;
        
        % Storage
        speed(index_g) = norm(Xk_post(4:6));
        [r(index_g,1),r(index_g,2),r(index_g,3)] = ecef2enu(Xk_post(1),Xk_post(2),Xk_post(3),recPos0(1),recPos0(2),recPos0(3),wgs84);
        v(index_g,:) = ((R_ecef2enu')*Xk_post(4:6))';
        tu(index_g) = Xk_post(10);
        fu(index_g) = Xk_post(11);
        td(index_g) = Xk_post(end);
        index_g = index_g + 1;
    end
end

save("ekf_td.mat","time","td");
save("ekf_pos.mat","r");






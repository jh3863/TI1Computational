% TI1_01_PS02
% =========================================================================
% Author: ...
% Date: 231119
% Version: 1.0 231119 JH Initial Release
%            -
%
% Source:
%
% Description:
%       - Computational 
%
% Required Input:
%       - No path definitions needed. Set path to folder containing this
%       script and the plots subfolder
%
% Output:
%
% Improvements:
%       - move VFI (& plotting) into separate scripts/functions
%       -
%
%=========================================================================


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%% 00a SETUP %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define model parameters
beta = 0.99;
z    = 1.00;
alpha= 0.70;

% define depreciation rate
delta= 0.5;

% define computation parameters
tol = 10^-6;
gNk   = 250;
iterMax = 10000;

% define the utility function
u =@(c) log(c)

% Create k grid
kbar = [0.1 (z./(delta))^(1/(1-alpha)) ];
gridK= [kbar(1):(kbar(2)-kbar(1))/(gNk-1):kbar(2)]';
gridZ= [1 0.5]'; gNz   = size(gridZ,1);

% definex C for a given k and k' using the budget constraint
cGivenKK = @(k,kprime,z) z.*k.^alpha + (1-delta).*k - kprime;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%% 01B Gamma = 1, gNk=250 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gamma = 1.0 ;
Q = [gamma (1-gamma); ...
     gamma (1-gamma); ...
    ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% 01B1 CALCULATE value & policy function %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize value function iteration (VFI)
Vdiff=10; iter=0; 

% define initial value function guess
V = 0.*meshgrid(gridZ,gridK);
Vnew = 0.5*V
% define objects to track v_n at each iteration
VdiffH = nan(1,iterMax);

polKindex = 0.*meshgrid(gridZ,gridK);

% run VFI
tic 
while Vdiff>tol & iter <iterMax
    
    for iZ =1:gNz
        cfeasible = cGivenKK(gridK,gridK',gridZ(iZ))>=0;
        c = cGivenKK(gridK,gridK',gridZ(iZ)) .* cfeasible;
        [VnewCol , index] =max( u(c) + beta * Q(iZ,:) * V' ,[], 2   );  %max( u(c) + beta * V * Q(iZ,:)',[], 2   );
        Vnew(:,iZ) = VnewCol;
        polKindex(:,iZ) = index;
    end
    %Vdiff = sum((Vnew- V).^2,[1 2]); V= Vnew; iter = iter+1;
    Vdiff = max(abs(Vnew- V),[],[1 2]); V= Vnew; iter = iter+1;

    % store history
    VdiffH(iter) = Vdiff./tol; %Vhistory(:,iter) = V;
end
VFIps02_time = toc;
%VdiffH(1000:iter+1)
iter

load('dta/PS01.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% 01B2 PLOT value & policy function %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PS02.Q1b.V = V;
PS02.Q1b.polK = gridK(polKindex(:,:));

figure; set(gcf,'position',[400,800,1500,600]);
p=subplot(1,2,1)
plot(gridK, PS01.V, 'LineWidth',4); grid on; hold on; ylabel('k_{t+1}'); xlabel('k_{t}');
plot(gridK, V(:,1),'--', 'LineWidth',4); grid on; 
title('Value Function V(k_t) '); p.FontSize =25;
legend([{'V^d PS02'} {'V^d PS01'} ],'Location','southeast');

p = subplot(1,2,2);
plot(gridK, gridK(polKindex(:,1)) , 'LineWidth',4); grid on; hold on; ylabel('k_{t+1}'); xlabel('k_{t}');
plot(gridK, PS01.polKindex,'--' , 'LineWidth',2);
title('Policy Function k_{t+1}(k_t) '); p.FontSize =25;
legend([{'k_{t+1} PS02'} {'k_{t+1} PS01'} ],'Location','southeast');

saveas(gcf,'plots/PS02B.png');

max(abs([PS01.V - V(:,1)])) % 0.0062
max(abs([PS01.polKindex -   gridK(polKindex(:,1))] )) % 0.0062

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%% 01C Stochastic DP: Gamma = 0.5, gNk=250  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gamma = 0.5 ;
Q = [gamma (1-gamma); ...
     gamma (1-gamma); ...
    ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% 01C CALCULATE value & policy function %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize value function iteration (VFI)
Vdiff=10; iter=0; 

% define initial value function guess
V = 0.*meshgrid(gridZ,gridK);
Vnew = 0.5*V
% define objects to track v_n at each iteration
VdiffH = nan(1,iterMax);

polKindex = 0.*meshgrid(gridZ,gridK);

% run VFI
tic 
while Vdiff>tol & iter <iterMax
    
    for iZ =1:gNz
        cfeasible = cGivenKK(gridK,gridK',gridZ(iZ))>=0;
        c = cGivenKK(gridK,gridK',gridZ(iZ)) .* cfeasible;
        [VnewCol , index] =max( u(c) + beta * Q(iZ,:) * V' ,[], 2   );  %max( u(c) + beta * V * Q(iZ,:)',[], 2   );
        Vnew(:,iZ) = VnewCol;
        polKindex(:,iZ) = index;
    end
    %Vdiff = sum((Vnew- V).^2,[1 2]); V= Vnew; iter = iter+1;
    Vdiff = max(abs(Vnew- V),[],[1 2]); V= Vnew; iter = iter+1;

    % store history
    VdiffH(iter) = Vdiff./tol; %Vhistory(:,iter) = V;
end
VFIps02_time = toc
%VdiffH(1000:iter+1)
iter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% 01B2 PLOT value & policy function %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PS02.Q1c.V = V;
PS02.Q1c.Vcomb = gamma.*V(:,1) + (1-gamma).*V(:,2);
PS02.Q1c.Vcomb
PS02.Q1c.polK = gridK(polKindex(:,:));

figure; set(gcf,'position',[400,800,1500,600]);
p=subplot(1,2,1)
plot(gridK, PS02.Q1c.Vcomb,'b-', 'LineWidth',4); grid on; hold on; ylabel('k_{t+1}'); xlabel('k_{t}');
plot(gridK, PS02.Q1b.V(:,1),'r--', 'LineWidth',4); grid on; 
title('Value Function V(k_t) '); p.FontSize =25;
legend([{'V^s'} {'V^d'}],'Location','southeast');

p = subplot(1,2,2);
plot(gridK, PS02.Q1c.polK(:,1),'b-' , 'LineWidth',4); grid on; hold on; ylabel('k_{t+1}'); xlabel('k_{t}');
plot(gridK, PS02.Q1c.polK(:,2),'b-' , 'LineWidth',4); grid on; hold on; ylabel('k_{t+1}'); xlabel('k_{t}');

plot(gridK, PS02.Q1b.polK(:,1) ,'r--' , 'LineWidth',2);
plot(gridK, PS02.Q1b.polK(:,2) ,'r--' , 'LineWidth',2);
title('Policy Function k_{t+1}(k_t) '); p.FontSize =25;
legend([{'V^s(1.0)'} {'V^s(0.5)'} {'V^d(1.0)'} {'V^d(0.5)'}],'Location','southeast');

saveas(gcf,'plots/PS02C.png');

[PS02.Q1c.polK(:,1)  PS02.Q1b.polK(:,1)]
[PS02.Q1c.polK(:,2)  PS02.Q1b.polK(:,2)]
max(abs([PS01.V - V(:,1)])) % 0.0062
max(abs([PS01.polKindex -   gridK(polKindex(:,1))] )) % 0.0062

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%% 01D Simulate Path %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set seed
rng(7);

k0 = 5;
Ntimesteps = 300;
Npaths = 25;
zdraws    = rand(Ntimesteps,Npaths);
zdrawsBin = zdraws < 0.5;

kpath = nan(Ntimesteps+1,Npaths);
kindexpath = nan(Ntimesteps+1,Npaths);
kpath(1,:) = k0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% 01D1 Simluate paths %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% simulate for all paths in parallel
for iTime = 1:Ntimesteps
    [~, kindex ]= min(abs(gridK-kpath(iTime,:)));
    kpath(iTime+1,:)      = PS02.Q1c.polK( sub2ind([gNk gNz], kindex  , zdrawsBin(iTime,:)+1));
end

%
[kpath [zdrawsBin(1,:); zdrawsBin]]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% 01D2 PLOT simluated paths %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; set(gcf,'position',[400,800,1500,600]);
p=subplot(1,2,1)
plot([1:50+1],kpath(1:51,1), 'LineWidth',4); grid on; hold on; ylabel('k_{t+1}'); xlabel('t');
title('Single simulated path'); p.FontSize =25;

p = subplot(1,2,2);
plot([1:Ntimesteps+1],kpath(:,:)); grid on; ylabel('k_{t+1}'); xlabel('t');
title('Many simualted paths'); p.FontSize =25;

saveas(gcf,'plots/PS02D.png');

% mean after 50 time periods
mean(kpath(50:end,:))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%% 01E Marko Shocks %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%% 01F Tauchen Method for AR(1) approximation %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gNz = 5; mu=2.5; rho=0.95; sigma=0.150; m=2.5;
[gridZ, Q ]= tauchen(gNz,mu,rho,sigma,m);
gridZ

%gridZ = gridZ


% define wider capital grid
gNk   = 500;
kbar = [0.01 (max(gridZ)./(delta))^(1/(1-alpha))/6 ];
gridK= [kbar(1):(kbar(2)-kbar(1))/(gNk-1):kbar(2)]';

%cGivenKK = @(k,kprime,z) z.*k.^alpha + (1-delta).*k - kprime;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% 01F CALCULATE value & policy function %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize value function iteration (VFI)
Vdiff=10; iter=0; 

% define initial value function guess
V = 0.*meshgrid(gridZ,gridK);
Vnew = 0.5*V
% define objects to track v_n at each iteration
VdiffH = nan(1,iterMax);

polKindex = 0.*meshgrid(gridZ,gridK);

% run VFI
tic 
while Vdiff>tol & iter <iterMax
    
    for iZ =1:gNz
        cfeasible = cGivenKK(gridK,gridK',gridZ(iZ))>=0;
        c = cGivenKK(gridK,gridK',gridZ(iZ)) .* cfeasible;
        [VnewCol , index] =max( u(c) + beta * Q(iZ,:) * V' ,[], 2   );  %max( u(c) + beta * V * Q(iZ,:)',[], 2   );
        Vnew(:,iZ) = VnewCol;
        polKindex(:,iZ) = index;
    end
    %Vdiff = sum((Vnew- V).^2,[1 2]); V= Vnew; iter = iter+1;
    Vdiff = max(abs(Vnew- V),[],[1 2]); V= Vnew; iter = iter+1;

    % store history
    VdiffH(iter) = Vdiff./tol; %Vhistory(:,iter) = V;
end
VFIps02_time = toc
%VdiffH(1000:iter+1)
iter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% 01F2 PLOT value & policy function %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PS02.Q1f.V = V;
PS02.Q1f.polK = gridK(polKindex(:,:));

% set seed
rng(7);

k0 = 5;
Ntimesteps = 300;
Npaths = 25;
zdraws    = rand(Ntimesteps,Npaths);
zdrawsBin = zdraws < 0.5;

kpath = nan(Ntimesteps+1,Npaths);
kindexpath = nan(Ntimesteps+1,Npaths);
kpath(1,:) = k0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% 01D1 Simluate paths %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% simulate for all paths in parallel
for iTime = 1:Ntimesteps
    [~, kindex ]= min(abs(gridK-kpath(iTime,:)));
    kpath(iTime+1,:)      = PS02.Q1f.polK( sub2ind([gNk gNz], kindex  , zdrawsBin(iTime,:)+1));
end

%
[kpath [zdrawsBin(1,:); zdrawsBin]]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% 01D2 PLOT simluated paths %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; set(gcf,'position',[400,800,1500,600]);
p=subplot(1,2,1)
plot([1:50+1],kpath(1:51,1), 'LineWidth',4); grid on; hold on; ylabel('k_{t+1}'); xlabel('t');
title('Single simulated path'); p.FontSize =25;

p = subplot(1,2,2);
plot([1:Ntimesteps+1],kpath(:,1:3)); grid on; ylabel('k_{t+1}'); xlabel('t');
title('Many simualted paths'); p.FontSize =25;

saveas(gcf,'plots/PS02F.png');


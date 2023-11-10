% TI1_00_main
% =========================================================================
% Author: Jacob Hartwig (jhartwig@uchicago.edu)
% Date: 231109
% Version: 1.0 231109 JH Initial Release
%            -
%
% Source:
%
% Description:
%       - setup workspace environment
%
% Required Input:
%
% Output:
%
% Improvements:
%       -
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
delta= 1;

% define computation parameters
tol = 10^-6;
gN   = 500;
iterMax = 10000;
Vhistory = nan(gN,iterMax);

% define the utility function
u =@(c) log(c)

% Create k grid
kbar = [0.1 (z./(delta))^(1/(1-alpha)) ];
gridK= [kbar(1):(kbar(2)-kbar(1))/(gN-1):kbar(2)]'

% definex C for a given k and k' using the budget constraint
cGivenKK = @(k,kprime) z.*k.^alpha + (1-delta).*k - kprime

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%% 01A Delta = 1, gN=500 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define constants as derived in Q1
A = 1/(1-beta).*(log(z) + log(1-alpha*beta)+alpha * beta/ (1-alpha * beta) *log(alpha * beta));
B = alpha / (1-alpha * beta);

% define theoretical value & policy function function
Vtrue = @(k) A + B * log(k)
Ptrue = @(k) (beta*B*z*k.^alpha)./(beta*B+1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% 01A CALCULATE value & policy function %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize value function iteration (VFI)
Vdiff=10; iter=0; 

% define initial value function guess
V = 0.* gridK;

% define objects to track v_n at each iteration
Vhistory(:,1) = V; VdiffH = nan(1,iterMax);

% run VFI
while Vdiff>tol & iter <iterMax
    cfeasible = cGivenKK(gridK,gridK')>=0;
    c = cGivenKK(gridK,gridK') .* cfeasible;
    [Vnew , index] = max( u(c) + beta* V',[], 2   );
    Vdiff = sum((Vnew- V).^2); V= Vnew; iter = iter+1;

    % store history
    VdiffH(iter) = Vdiff./tol; Vhistory(:,iter) = V;
end
%VdiffH(1000:iter+1)
iter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% 01A PLOT value & policy function %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; set(gcf,'position',[400,800,1500,600]);
p =subplot(1,2,1);
plot(gridK, Vtrue(gridK),'--', 'LineWidth',8); grid on; hold on; ylabel('k_{t+1}'); xlabel('k_{t}');
plot(gridK, Vhistory(:,[iter]), 'LineWidth',4); title('Value Function V(k_t) '); p.FontSize =30;


p = subplot(1,2,2);
plot(gridK, Ptrue(gridK),'--', 'LineWidth',8); grid on; hold on; ylabel('k_{t+1}'); xlabel('k_{t}');
plot(gridK, gridK(index) , 'LineWidth',4); title('Policy Function k_{t+1}(k_t) '); p.FontSize =30;
legend([{'Theoretial'} {'VFI'} ],'Location','southeast');

saveas(gcf,'plots/Q2a.png');

% [Vtrue(gridK)   Vhistory(:,[iter])]
% [Ptrue(gridK)   gridK(index)]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%% 01B Delta = 0.5, gN=500 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gN   = 2000; Vhistory = nan(gN,iterMax);
% define depreciation rate
delta= 0.5;

% update k grid
kbar = [0.1 (z./(delta))^(1/(1-alpha)) ];
gridK= [kbar(1):(kbar(2)-kbar(1))/(gN-1):kbar(2)]'

% define C for a given k and k' using the budget constraint
cGivenKK = @(k,kprime) z.*k.^alpha + (1-delta).*k - kprime

% theoretical steady state values
theoryK_steadystate = (alpha*z*beta/(1-beta+beta*delta))^(1/(1-alpha));
theoryC_steadystate = z*theoryK_steadystate^alpha-delta*theoryK_steadystate;

%(alpha*z*beta/(1-beta+beta*delta))^(1/(1-alpha));

% initialize value function iteration (VFI)
Vdiff=10; iter=0; 

% define initial value function guess
V = 0.* gridK;

% define objects to track v_n at each iteration
Vhistory(:,1) = V; VdiffH = nan(1,iterMax);

% run VFI
while Vdiff>tol & iter <iterMax
    cfeasible = cGivenKK(gridK,gridK')>=0;
    c = cGivenKK(gridK,gridK') .* cfeasible;
    [Vnew , index] = max( u(c) + beta* V',[], 2   );
    Vdiff = sum((Vnew- V).^2); V= Vnew; iter = iter+1;

    % store history
    VdiffH(iter) = Vdiff./tol; Vhistory(:,iter) = V;
end
VdiffH(iter-5:iter+1)
iter

figure; set(gcf,'position',[400,800,1500,600]);
p =subplot(1,2,1);
plot(gridK, Vhistory(:,[iter]), 'LineWidth',4); grid on; hold on; ylabel('k_{t+1}'); xlabel('k_{t}');
 title('Value Function V(k_t) '); p.FontSize =30;


p = subplot(1,2,2);
plot(gridK, gridK(index) , 'LineWidth',4); grid on; hold on; ylabel('k_{t+1}'); xlabel('k_{t}');
plot(gridK, gridK,'--' , 'LineWidth',2);
title('Policy Function k_{t+1}(k_t) '); p.FontSize =30;
legend([{'VFI'} {'45°'} ],'Location','southeast');

saveas(gcf,'plots/Q2b.png');


[gridK gridK(index) gridK-gridK(index)]

[~,index_steadyK]= min(abs(gridK-gridK(index)));
kplus = gridK(index); K_steadystate = kplus(index_steadyK)
[K_steadystate theoryK_steadystate]

[cGivenKK(K_steadystate,K_steadystate) theoryC_steadystate]

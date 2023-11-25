% TI1_00_main
% =========================================================================
% Author: ...
% Date: 231109
% Version: 1.0 231109 JH Initial Release
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

% define depreciation rate
delta= 0.5;

% define C for a given k and k' using the budget constraint
cGivenKK = @(k,kprime) z.*k.^alpha + (1-delta).*k - kprime

% theoretical steady state values
theoryK_steadystate = (alpha*z*beta/(1-beta+beta*delta))^(1/(1-alpha));
theoryC_steadystate = z*theoryK_steadystate^alpha-delta*theoryK_steadystate;


% update k grid
gN   = 250; Vhistory = nan(gN,iterMax);
kbar = [0.1 (z./(delta))^(1/(1-alpha)) ];
gridK= [kbar(1):(kbar(2)-kbar(1))/(gN-1):kbar(2)]'


% initialize value function iteration (VFI)
Vdiff=10; iter=0; 

% define initial value function guess
V = 0.* gridK;

% define objects to track v_n at each iteration
Vhistory(:,1) = V; VdiffH = nan(1,iterMax);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% 01B CALCULATE value & policy function%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
while Vdiff>tol & iter <iterMax
    cfeasible = cGivenKK(gridK,gridK')>=0;
    c = cGivenKK(gridK,gridK') .* cfeasible;

    % new value function & index of optimal capital choice (ie policy funtion)
    [Vnew , index] = max( u(c) + beta* V',[], 2   );
    Vdiff = sum((Vnew- V).^2); V= Vnew; iter = iter+1;

    % store history
    VdiffH(iter) = Vdiff./tol; Vhistory(:,iter) = V;
end
toc
VdiffH(iter-5:iter+1)
iter


PS01.polKindex = gridK(index)
PS01.V = Vhistory(:,[iter]);

save('dta/PS01.mat',"PS01")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% 01B PLOT value & policy function %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% 01B COMPUTE steady states %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[gridK gridK(index) gridK-gridK(index)]

[~,index_steadyK]= min(abs(gridK-gridK(index)));
kplus = gridK(index); K_steadystate = kplus(index_steadyK);
K_ss_diff = K_steadystate-theoryK_steadystate;
[K_steadystate theoryK_steadystate  K_ss_diff K_ss_diff/theoryK_steadystate]

C_ss_diff = cGivenKK(K_steadystate,K_steadystate)-theoryC_steadystate;
[cGivenKK(K_steadystate,K_steadystate) theoryC_steadystate C_ss_diff C_ss_diff/theoryC_steadystate]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%% 01C Delta = 0.5,Loop over grid sizes %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define grid sizes
gNgrid = [50, 100, 500, 1000, 2000];
gNgridLength = size(gNgrid,2);

% create arrays to store run time and steady states
gNgridTime = nan(1,gNgridLength);
gNgridKss  = nan(1,gNgridLength);
gNgridCss  = nan(1,gNgridLength);

% loop over different grid sizes
for iNgrid=1:gNgridLength
    gN   = gNgrid(iNgrid); 

    % update k grid
    Vhistory = nan(gN,iterMax);
    kbar = [0.1 (z./(delta))^(1/(1-alpha)) ];
    gridK= [kbar(1):(kbar(2)-kbar(1))/(gN-1):kbar(2)]';
    
    
    % initialize value function iteration (VFI)
    Vdiff=10; iter=0; 
    
    % define initial value function guess
    V = 0.* gridK;
    
    % define objects to track v_n at each iteration
    Vhistory(:,1) = V; VdiffH = nan(1,iterMax);
    
    % run VFI
    tic
    while Vdiff>tol & iter <iterMax
        cfeasible = cGivenKK(gridK,gridK')>=0;
        c = cGivenKK(gridK,gridK') .* cfeasible;
        [Vnew , index] = max( u(c) + beta* V',[], 2   );
        Vdiff = sum((Vnew- V).^2); V= Vnew; iter = iter+1;
    
        % store history
        VdiffH(iter) = Vdiff./tol; Vhistory(:,iter) = V;
    end
    gNgridTime(iNgrid) = toc

    % extract steady state values
    [~,index_steadyK]= min(abs(gridK-gridK(index)));
    kplus = gridK(index); gNgridKss(iNgrid) = kplus(index_steadyK);

    gNgridCss(iNgrid) = cGivenKK(gNgridKss(iNgrid),gNgridKss(iNgrid));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% 01C PLOT run time & steady state deviations %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; set(gcf,'position',[400,800,1500,600]);
p =subplot(1,3,1);
%plot(gNgrid, gNgridTime,'O-', 'LineWidth',4); grid on; hold on; ylabel('time'); xlabel('grid size');
plot(gNgrid, gNgridTime,'O-', 'LineWidth',4); grid on; hold on; ylabel('time (in seconds)'); xlabel('grid size');
title('A) Computational Time'); p.FontSize =20;


p = subplot(1,3,2);
plot(gNgrid, gNgridKss , 'bO-', 'LineWidth',4); grid on; hold on; ylabel('k_{ss} & c_{ss}'); xlabel('grid size');
plot(gNgrid, ones(1,gNgridLength)*theoryK_steadystate , 'k--', 'LineWidth',2); 
plot(gNgrid, gNgridCss , 'rO-', 'LineWidth',4); grid on; hold on;
plot(gNgrid, ones(1,gNgridLength)*theoryC_steadystate , 'k--', 'LineWidth',2); 
title('B) Steady state comparison'); p.FontSize =20;
%legend([{'k_{ss}'} {''} {'c_{ss}'} {''}],'Location','southeast');

p = subplot(1,3,3);
plot(gNgrid, abs(gNgridKss-theoryK_steadystate), 'bO-', 'LineWidth',4); grid on; hold on; ylabel('|x_{ss}-x_{ss}^T|'); xlabel('grid size');
plot(gNgrid, abs(gNgridCss-theoryC_steadystate), 'rO-', 'LineWidth',4); grid on; hold on;
title('C) Steady state absolute difference'); p.FontSize =20;
legend([{'k_{ss}'} {'c_{ss}'}],'Location','southeast');

saveas(gcf,'plots/Q2c.png');

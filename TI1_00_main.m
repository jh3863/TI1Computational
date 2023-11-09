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

tol = 10^-6;

beta = 0.99;
z    = 1.00;
alpha= 0.70;

delta= 1;
gN   = 500;
iterMax = 10000;
Vhistory = nan(gN,iterMax);


u =@(c) log(c)

% Create k grid
kbar = [0.1 (z./(delta))^(1/(1-alpha)) ];
gridK= [kbar(1):(kbar(2)-kbar(1))/(gN-1):kbar(2)]'



cGivenKK = @(k,kprime) z.*k.^alpha + (1-delta).*k - kprime



%kprime = @(k,c) z.*k.^alpha + (1-delta).*k - c;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%% 01A Delta = 1, gN=500 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = 1/(1-beta).*(log(z) + log(1-alpha*beta)+alpha * beta/ (1-alpha * beta) *log(alpha * beta));
B = alpha / (1-alpha * beta);

Vtrue = @(k) A + B * log(k)
Ptrue = @(k) (beta*B*z*k.^alpha)./(beta*B+1)

%kprime(0.5,0)
%u(cGivenKK(gridK,gridK'))

V = 0.* gridK;
Vhistory(:,1) = V;
Vdiff=10; iter=0; VdiffH = nan(1,iterMax);
while Vdiff>tol & iter <iterMax
    cfeasible = cGivenKK(gridK,gridK')>=0;
    c = cGivenKK(gridK,gridK') .* cfeasible;
    [Vnew , index] = max( u(c) + beta* V',[], 2   );
    Vdiff = sum((Vnew- V).^2); V= Vnew; iter = iter+1;

    % store history
    VdiffH(iter) = Vdiff./tol; Vhistory(:,iter) = V;
end
VdiffH(1000:iter+1)
iter

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





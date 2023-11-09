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
%% %%% 00a SETUP %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = 1/(1-beta).*(log(z) + log(1-alpha*beta)+alpha * beta/ (1-alpha * beta) *log(alpha * beta));
B = alpha / (1-alpha * beta);


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


figure
subplot(1,2,1);

figure
plot(gridK, A + B*log(gridK), 'LineWidth',4); grid on; hold on;
plot(gridK, Vhistory(:,[iter]), 'LineWidth',2); grid on; hold on;

[A + B*log(gridK)   Vhistory(:,[iter])]
%plot(gridK, Vhistory(:,[1:100:1000 iter]), 'LineWidth',2); grid on; hold on;

title('Value Function')

subplot(1,2,2);
plot(gridK, gridK(index) , 'LineWidth',2); grid on; title('K` Policy Function')




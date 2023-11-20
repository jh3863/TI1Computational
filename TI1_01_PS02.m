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
gN   = 250;
iterMax = 10000;

% define the utility function
u =@(c) log(c)

% Create k grid
kbar = [0.1 (z./(delta))^(1/(1-alpha)) ];
gridK= [kbar(1):(kbar(2)-kbar(1))/(gN-1):kbar(2)]';
gridZ= [1 0.5]';

% definex C for a given k and k' using the budget constraint
cGivenKK = @(k,kprime) z.*k.^alpha + (1-delta).*k - kprime

gamma = 0.5 ;
Q = [gamma (1-gamma); ...
     gamma (1-gamma); ...
    ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%% 01B Gamma = 1, gN=250 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% 01A CALCULATE value & policy function %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize value function iteration (VFI)
Vdiff=10; iter=0; 

% define initial value function guess
V = 0.*meshgrid(gridZ,gridK);

% define objects to track v_n at each iteration
VdiffH = nan(1,iterMax);

% run VFI
while Vdiff>tol & iter <iterMax
    cfeasible = cGivenKK(gridK,gridK')>=0;
    c = cGivenKK(gridK,gridK') .* cfeasible;
    [Vnew , index] = max( u(c) + beta* Q * V',[], 2   );

    Vdiff = sum((Vnew- V).^2); V= Vnew; iter = iter+1;

    % store history
    VdiffH(iter) = Vdiff./tol; Vhistory(:,iter) = V;
end
%VdiffH(1000:iter+1)
iter


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%% 01C Stochastic DP: Gamma = 0.5, gN=250  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%% 01D Simulate Path %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%% 01E Marko Shocks %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%% 01F Tauchen Method for AR(1) approximation %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


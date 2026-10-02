close all; clear all;

%% --- Base parameters (north) ---
Lambda = [300; 50; 659.8; 27; 36.1; 0.32];

mu = zeros(length(Lambda),1);
b_hat_zL = [73.4; 44.5; 87.8; 165.4; 15; 1963.3];
b_hat_zN = [3.7; 34.2; 0.35; 277.3; 0.36; 750].*[1; 1; 3.59; 1; 1.03; 1.25];
 
base_c = 0.1815;
b_hat_L = b_hat_zL./[base_c; base_c; base_c; base_c; 1-base_c; base_c];
b_hat_N = b_hat_zN./[base_c; base_c; base_c; base_c; 1-base_c; base_c];
bA = 1037.5*1.08;

betaL = [0.91; 0.57; 0.56; 0.06; 0.024];
betaN = [0.91; 0.57; 0.56; 0.06; 0.024];
alpha = [1;1;1;1;1];

k = [2; 2; 2; 2; 1];
r = 2.52e6; 
K = 5;
sL = 0.43; sN = 0.55; sA = 0.84;
n_hosts = length(Lambda);

%% Initital Conditions (north)
L0 = 115000;
N0 = 11500;
A0 = 3450;
H(:,1) = [40; 30; 59.2; 6; 27; 0.46];
% %% --- Base parameters (south) ---
% Lambda = [273.5; 50; 174.9; 27; 1278.1; 0.21];
% 
% mu = zeros(length(Lambda),1);
% b_hat_zL = [3.1; 4.5; 19; 16.54; 25.2; 196.3];
% b_hat_zN = [0.25; 3.42; 0.035; 27.73; 2.5; 75.0].*[1.24; 1; 0.83; 1; 1.19; 1.00];
% 
% base_c = 0.0204;
% b_hat_L = b_hat_zL./[base_c; base_c; base_c; base_c; 1-base_c; base_c];
% b_hat_N = b_hat_zN./[base_c; base_c; base_c; base_c; 1-base_c; base_c];
% bA = 1037.5*1.00;
% 
% betaL = [0.91; 0.57; 0.56; 0.06; 0.24];
% betaN = [0.91; 0.57; 0.56; 0.06; 0.24];
% alpha = [1;1;1;1;1];
% 
% sL = 0.43; sN = 0.55; sA = 0.84;
% k = [2; 2; 2; 2; 1];
% n_hosts = length(Lambda);
% r = 2.52e6; K=5;
% %% Initital Conditions (south)
% L0 = 115000;
% N0 = 11500;
% A0 = 3450;
% H(:,1) = [27.9; 30; 30; 6; 445.6; 0.40];

%% --- Compute mu and z values ---

for i = 1:n_hosts-1
    mu(i) = k(i)*log(Lambda(i)/(k(i)*H(i,1))+1);
end
mu(n_hosts) = 2*log((Lambda(n_hosts)+sqrt(Lambda(n_hosts)^2+4*H(n_hosts,1)^2)) / (2*H(n_hosts,1)));

SLsum0 = 0; SNsum0 = 0; 
for i = 1:n_hosts-1
    SLsum0 = SLsum0 + b_hat_zL(i)*(H(i,1)+Lambda(i)/k(i));
    SNsum0 = SNsum0 + b_hat_zN(i)*H(i,1);
end
SLsum0 = SLsum0 + b_hat_zL(n_hosts)*(H(n_hosts,1)*exp(-mu(n_hosts)/2)+Lambda(n_hosts));
SNsum0 = SNsum0 + b_hat_zN(n_hosts)*H(n_hosts,1)*exp(-mu(n_hosts)/2);
SAsum0 = bA*H(n_hosts,1);

zL = -L0*log(1-0.1/sL)/SLsum0;
zN = -N0*log(1-0.1/sN)/SNsum0;
zA = -A0*log(1-0.3/sA)/SAsum0;

bL = zeros(n_hosts,1); bN = zeros(n_hosts,1);
idx = [1:(n_hosts-2) n_hosts];
for i = idx
    bL(i) = b_hat_L(i) * base_c;
    bN(i) = b_hat_N(i) * base_c;
end
bL(n_hosts-1) = b_hat_L(n_hosts-1) * (1-base_c);
bN(n_hosts-1) = b_hat_N(n_hosts-1) * (1-base_c);

SLsum = 0; SNsum = 0; 
for i = 1:n_hosts-1
    SLsum = SLsum+bL(i)*(H(i)+Lambda(i)/k(i));
    SNsum = SNsum+bN(i)*H(i);
end
SLsum = SLsum+bL(n_hosts)*(H(n_hosts)*exp(-mu(n_hosts)/2)+Lambda(n_hosts));
SNsum = SNsum+bN(n_hosts)*H(n_hosts)*exp(-mu(n_hosts)/2);
SAsum = bA*H(n_hosts);
%% Loop parameter
% c
c_values = [0,base_c,1];
n_c = length(c_values);

% % z
% z_sets = [
%             %zL*0.5, zN*0.5, zA*0.5;
%             zL, zN, zA;
%             %zL*1.5, zN*1.5, zA*1.5
%             ];

% s
% used 0.1, normal, and 1 for south s
%s_sets = [sL*1e-1, sN*1e-1, sA*1e-1;
    %sL, sN, sA;
%    1, 1, 1];
% s_sets = [sL*0.15, sN*0.15, sA*0.15];

% H6
% H6_sets = [0.01; 2];
%% Outer loop

if isempty(gcp('nocreate'))
    parpool;   % start parallel pool only if one doesn't exist
end

for c_idx = 1:n_c
    c = c_values(c_idx); 

% for z_idx = 1:size(z_sets,1)
%     zL = z_sets(z_idx,1);
%     zN = z_sets(z_idx,2);
%     zA = z_sets(z_idx,3);
% 
%     c=base_c;

% % vary sL, sN, sA
% for s_idx = 1:size(s_sets,1)
%     sL = s_sets(s_idx,1);
%     sN = s_sets(s_idx,2);
%     sA = s_sets(s_idx,3);
% 
%       c=base_c;

% for H6_idx = 1:size(H6_sets,1)
%     H(n_hosts,1) = H6_sets(H6_idx,1);
% 
%     c=base_c;
%% Loop values
r_values = linspace(0,3e6,50);
K_values = linspace(1,200,50);
n_r = length(r_values);
n_K = length(K_values);

%% Storage
stability_map = nan(n_K,n_r);
stable_fp1 = cell(n_K, n_r);
stable_fp2 = cell(n_K, n_r);
unstable_fp1 = cell(n_K, n_r);
unstable_fp2 = cell(n_K, n_r);

unique_positive_fp1 = cell(n_K,n_r);  % stores unique postive fixed point from system 1 for each c
nonmatching_fp2 = cell(n_K,n_r);      % stores fixed points from system 2 that do not match the system 1 solutions
tol_match = 1e-6;  % tolerance for defining matching fixed points

unique_solutions1 = []; % will hold all unique soltuions to first gen system
unique_solutions2 = []; % will hold all unique solutions to second gen system

options = optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-5, 'TolX', 1e-5); % 


L_range = [0, 0.5e1, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e10];
N_range = [0, 0.5e1, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e10]; 
A_range = [0, 0.5e1, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e10];


[L0_grid,N0_grid,A0_grid] = ndgrid(L_range,N_range,A_range);

x0s = [L0_grid(:), N0_grid(:), A0_grid(:)];
n_guesses = size(x0s,1);

if isempty(gcp('nocreate'))
    parpool;   % start parallel pool only if one doesn't exist
end

% Loop over values
for iK = 1:n_K
    for ir = 1:n_r
        iK
        ir
        r = r_values(ir);
        K = K_values(iK);
        params = {r, K, sL, sN, sA, H(1:n_hosts),zL, zN, zA, Lambda, SLsum, SNsum, SAsum};

        solutions1 = []; % store all solutions to system 1 for this c
        solutions2 = []; % store all solutions to system 2 for this c

        parfor q = 1:n_guesses
                    x0 = x0s(q,:)';
            
                    [x_sol1, fval1, exitflag1] = fsolve(@(x) mySystem(x,params),x0,options);
                    [x_sol2, fval2, exitflag2] = fsolve(@(x) mySystem2(x,params),x0,options);
                
    
                    % Only keep valid solutions
                        if exitflag1 > 0  && all(x_sol1>=0) && norm(fval1) < 1e-6 % exitflag>0 means fsovle converged; all(x_sol>=0) keeps only non-negative solutions, norm(fval)<1e-6 ensures small residual
                            solutions1 = [solutions1; x_sol1(:)']; % stores solultions  
                        end
                        if exitflag2 > 0  && all(x_sol2>=0) && norm(fval2) < 1e-6       
                            solutions2 = [solutions2; x_sol2(:)'];
                        end
        end      
%%
    unique_solutions1 = unique(round(solutions1, 3), 'rows');      
    unique_solutions2 = unique(round(solutions2, 3), 'rows');
        
    % --- Compare all system 1 solutions with system 2 solutions (no filtering) ---

    % If system 1 has solutions, use them as references
    if ~isempty(unique_solutions1)
        ref_fp = unique_solutions1;   % <--- no filtering, take all solutions
    else
        ref_fp = [];                  % no reference solutions
    end

    % Now compare system 2 to system 1
    if ~isempty(unique_solutions2)

        if isempty(ref_fp)
            % Nothing to compare: keep all system 2 solutions
            nonmatching_fp2{iK,ir} = unique_solutions2;

        else
            % Compute differences between each system 2 solution and each system 1 solution
            % Using pdist2 so all pairwise distances are computed
            D = pdist2(unique_solutions2, ref_fp);

            % Keep system 2 points that do NOT match ANY system 1 point
            nonmatch_mask = min(D, [], 2) > 1e-2;

            nonmatching_fp2{iK,ir} = unique_solutions2(nonmatch_mask, :);

        end

    else
        nonmatching_fp2{iK,ir} = [];
    end

    stable1 = []; unstable1 = []; stable2 = []; unstable2 = [];
    % System 1 stability
    J_fun = makeJacobian(params);
    for s1 = 1:size(unique_solutions1,1)
        x = unique_solutions1(s1,:)';

        J = J_fun(x);
        eigvals = eig(J);
        if max(abs(eigvals))<1
            stable1 = [stable1; x'];
        else
            unstable1 = [unstable1; x'];
        end
    end

    % System 2 stability (ONLY nonmatching points)
    these_fp2 = nonmatching_fp2{iK,ir};       % <- filtered list
    for s = 1:size(these_fp2,1)

    % Use this point as the representative
    x1 = these_fp2(s,:)';

    % Generate the orbit in the correct order
    x2 = tickMap(x1,params);

    % Compute Jacobians
    J1 = J_fun(x1);
    J2 = J_fun(x2);

    % Jacobian of one trip around the cycle
    M = J2*J1;

    eigvals = eig(M);

    if max(abs(eigvals)) < 1
        stable2 = [stable2; x1'];
    else
        unstable2 = [unstable2; x1'];
    end
end

    stable_fp1{iK,ir}   = stable1;
    unstable_fp1{iK,ir} = unstable1;
    stable_fp2{iK,ir}   = stable2;
    unstable_fp2{iK,ir} = unstable2;
    

hasStable1 = ~isempty(stable1);
hasStable2 = ~isempty(stable2);
hasUnstable1 = ~isempty(unstable1);

stableExtinct = false;
stableCoexistFP = false;
stablesync2 = false;
stableasync2 = false;

% check for stable fixed point values
if hasStable1
    for kk = 1:size(stable1,1)

        x = stable1(kk,:);

        if all(x==0)
            stableExtinct = true;
        else
            stableCoexistFP = true;
        end

    end
end

if hasStable2
    if any(stable2(:)==0)
        stablesync2 = true;
    end
    if all(stable2(:)>0)
        stableasync2 = true;
    end

end
% --- assign colors ---
code = 0;

if stableExtinct
    code = code + 1;
end
if stableCoexistFP
    code = code + 2;
end
if stablesync2
    code = code + 4;
end
if stableasync2
    code = code + 8;
end

stability_map(iK,ir) = code;

    end
end

%% --- Plot results (2D) ---
S = 0.5*sL*sN*sA;
codes = [4 6 10 2 1];

labels = {'SC','FP+SC','FP+AC','FP','EE'};
plot_data = nan(size(stability_map));

for kk = 1:length(codes)
    plot_data(stability_map==codes(kk)) = kk;
end

cmap = [0.8 0.2 0.2 % red - stable sync cycle
        0.9 0.55 0.1 % green - stable fp, stable sync cycle
        0 0.6 0.55 % cyan - stable fp, stable async cycle
        0.2 0.4 0.75 % blue - stable fp 
        0.8 0.8 0.8 % gray - DNE
        ];

imagesc(r_values*S,K_values,plot_data)
set(gca,'YDir','normal')
hold on
r_boundary = K_values;
plot(r_boundary, K_values, 'k--', 'LineWidth', 1);

colormap(cmap)
clim([0.5 length(codes)+0.5])

colorbar('Ticks',1:length(codes), ...
         'TickLabels',labels)

xlabel('$rS$', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$K$', 'Interpreter', 'latex', 'FontSize', 14)

end
delete(gcp('nocreate'));
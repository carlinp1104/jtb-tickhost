close all;
clear all;

%% --- Parameter Setup ---

%% --- Base parameters (north) ---
Lambda = [300; 50; 659.8; 27; 36.1; 0.32];

mu = zeros(length(Lambda),1);
b_hat_zL = [73.4; 44.5; 87.8; 165.4; 28.9; 1963.3];
b_hat_zN = [3.7; 34.2; 0.35; 277.3; 0.36; 750].*[1; 1; 3.59; 1; 1.07; 1.25];

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
% b_hat_zL = [3.1; 4.5; 34.5; 16.54; 25.2; 196.3];
% b_hat_zN = [0.25; 3.42; 0.035; 27.73; 2.5; 75.0].*[1.24; 1; 0.90; 1; 1.19; 1.00];
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

%%
% range of c values
c_values = linspace(0,1,150);

% range of r values
% r_values = [linspace(0,2.5e6,200) linspace(2.5e6, 1e7,100)];

n_c = length(c_values);

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


%%
unique_solutions1 = []; % will hold all unique soltuions to first gen system
unique_solutions2 = []; % will hold all unique solutions to second gen system

opts = optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-8, 'TolX', 1e-8); % 

unique_positive_fp1 = cell(n_c,1);  % stores unique postive fixed point from system 1 for each c
nonmatching_fp2 = cell(n_c,1);      % stores fixed points from system 2 that do not match the system 1 solutions
tol_match = 1e-6;  % tolerance for defining matching fixed points
regime = strings(n_c,1);

% Range of initial guesses to try
L_range = [0, 0.5e1, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e10];
N_range = [0, 0.5e1, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e10]; 
A_range = [0, 0.5e1, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e10];

if isempty(gcp('nocreate'))
    parpool;   % start parallel pool only if one doesn't exist
end

% Loop over c values
for ii = 1:n_c
    ii
        c = c_values(ii);
        %c = base_c;

        bL = zeros(n_hosts,1); bN = zeros(n_hosts,1);
        idx = [1:(n_hosts-2) n_hosts];
        for i = idx
            bL(i) = b_hat_L(i) * c;
            bN(i) = b_hat_N(i) * c;
        end
        bL(n_hosts-1) = b_hat_L(n_hosts-1) * (1-c);
        bN(n_hosts-1) = b_hat_N(n_hosts-1) * (1-c);
    
        SLsum = 0; SNsum = 0; 
        for i = 1:n_hosts-1
            SLsum = SLsum+bL(i)*(H(i)+Lambda(i)/k(i));
            SNsum = SNsum+bN(i)*H(i);
        end
        SLsum = SLsum+bL(n_hosts)*(H(n_hosts)*exp(-mu(n_hosts)/2)+Lambda(n_hosts));
        SNsum = SNsum+bN(n_hosts)*H(n_hosts)*exp(-mu(n_hosts)/2);
        SAsum = bA*H(n_hosts);


        params = {r, K, sL, sN, sA, H(1:n_hosts),zL, zN, zA, Lambda, SLsum, SNsum, SAsum};

    solutions1 = []; % store all solutions to system 1 for this c
    solutions2 = []; % store all solutions to system 2 for this c

    parfor l = 1:length(L_range)
        for n = 1:length(N_range)
            for a = 1:length(A_range)

                x0 = [L_range(l); N_range(n); A_range(a)];

                [x_sol1, fval1, exitflag1] = fsolve(@(x) mySystem(x, params), x0, opts);
                [x_sol2, fval2, exitflag2] = fsolve(@(x) mySystem2(x, params), x0, opts); % finds fixed points of system 2
               

                % Only keep valid solutions
                if exitflag1 > 0 && all(x_sol1 >= 0) && norm(fval1) < 1e-6 % exitflag>0 means fsovle converged; all(x_sol>=0) keeps only non-negative solutions, norm(fval)<1e-6 ensures small residual
                    solutions1 = [solutions1; x_sol1(:)']; % stores solultions  
                end
                if exitflag2 > 0 && all(x_sol2 >= 0) && norm(fval2) < 1e-6
                    solutions2 = [solutions2; x_sol2(:)'];
                    x1 = x_sol2
                    x2 = tickMap(x1, params);
                    solutions2 = [solutions2; x2(:)']; 
                end
            end
        end
    end      

    unique_solutions1 = unique(round(solutions1, 2), 'rows');      
    unique_solutions2 = unique(round(solutions2, 2), 'rows');
        
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
            nonmatching_fp2{ii} = unique_solutions2;

        else
            % Compute differences between each system 2 solution and each system 1 solution
            % Using pdist2 so all pairwise distances are computed
            D = pdist2(unique_solutions2, ref_fp);

            % Keep system 2 points that do NOT match ANY system 1 point
            nonmatch_mask = min(D, [], 2) > 1e-2;

            nonmatching_fp2{ii} = unique_solutions2(nonmatch_mask, :);

        end

    else
        nonmatching_fp2{ii} = [];
    end
%%
    stable1 = []; unstable1 = [];
    stable2 = []; unstable2 = [];

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
    these_fp2 = nonmatching_fp2{ii};       % <- filtered list
    J_fun = makeJacobian(params);
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

    stable_fp1{ii}   = stable1;
    unstable_fp1{ii} = unstable1;
    stable_fp2{ii}   = stable2;
    unstable_fp2{ii} = unstable2;

end

%% --- Plot results (2D) ---
S = 0.5*sL*sN*sA;
varNames = {'Larvae (larvae/ha)', 'Nymphs (nymphs/ha)', 'Adults (adults/ha)'}; % variable labels
nVars = numel(varNames);

for v = 1:nVars
    figure; hold on;

    hStable = plot(nan, nan, 'b.', 'MarkerFaceColor','b','DisplayName', 'Stable');
    hUnstable = plot(nan, nan, 'r.', 'MarkerFaceColor','r','DisplayName', 'Unstable');

    xlabel('c', 'FontSize', 14);
    ylabel(varNames{v}, 'FontSize', 14);
    grid on;

    for ii = 1:n_c
        c = c_values(ii);

            % Plot unstable system 2 fps
            if ~isempty(unstable_fp2{ii})
                plot(c*ones(size(unstable_fp2{ii},1),1), unstable_fp2{ii}(:,v), 'r.', 'MarkerFaceColor', 'r');
            end

            % Plot unstable system 1 fps
            if ~isempty(unstable_fp1{ii})
                plot(c*ones(size(unstable_fp1{ii},1),1), unstable_fp1{ii}(:,v), 'r.', 'MarkerFaceColor', 'r');
            end      

            % Plot stable system 2 fps
            if ~isempty(stable_fp2{ii})
                plot(c*ones(size(stable_fp2{ii},1),1),stable_fp2{ii}(:,v), 'b.', 'MarkerFaceColor', 'b');
            end

            % Plot stable system 1 fps
            if ~isempty(stable_fp1{ii})
                plot(c*ones(size(stable_fp1{ii},1),1),stable_fp1{ii}(:,v), 'b.', 'MarkerFaceColor', 'b');
            end

    end
end


%% --- Plot results (3D) ---

figure; hold on; grid on;
xlabel('c', 'FontSize', 14);
ylabel('L', 'FontSize', 14);
zlabel('N', 'FontSize', 14);
view(3);
hStable3 = plot3(nan,nan,nan,'b.','MarkerFaceColor','b','DisplayName','Stable');
hUnstable3 = plot3(nan,nan,nan,'r.','MarkerFaceColor','r','DisplayName','Unstable');

for ii = 1:n_c
    c = c_values(ii);

    % Plot system 1 stable fp
    if ~isempty(stable_fp1{ii})
        plot3(c*ones(size(stable_fp1{ii},1),1), stable_fp1{ii}(:,1), stable_fp1{ii}(:,2), 'b.', 'MarkerFaceColor', 'b');
    end

    % Plot system 1 unstable fp
    if ~isempty(unstable_fp1{ii})
        plot3(c*ones(size(unstable_fp1{ii},1),1), unstable_fp1{ii}(:,1), unstable_fp1{ii}(:,2), 'r.', 'MarkerFaceColor', 'r');
    end

    % Plot system 2 stable fp
    if ~isempty(stable_fp2{ii})
        plot3(c*ones(size(stable_fp2{ii},1),1), stable_fp2{ii}(:,1), stable_fp2{ii}(:,2), 'b.', 'MarkerFaceColor', 'b');
    end

    % Plot system 2 unstable fp
    if ~isempty(unstable_fp2{ii})
        plot3(c*ones(size(unstable_fp2{ii},1),1), unstable_fp2{ii}(:,1), unstable_fp2{ii}(:,2), 'r.', 'MarkerFaceColor', 'r');
    end
end

legend([hStable3 hUnstable3], {'Stable','Unstable'}, 'FontSize', 14);
%%%%%%%%%%%%%%%%%%%%%%%%%%%
delete(gcp('nocreate'));
%% Run this after running DRIVER

% Plotter for cycles
n = length(c_values);

fp_branch_stable = nan(n,3);
fp_branch_unstable = nan(n,3);
fp_branch_zeros = nan(n,3);

lower_branch_stable = nan(n,3);
upper_branch_stable = nan(n,3);

lower_branch_unstable_sync = nan(n,3);
upper_branch_unstable_sync = nan(n,3);
lower_branch_unstable_async = nan(n,3);
upper_branch_unstable_async = nan(n,3);

for i = 1:n
    if ~isempty(stable_fp1{i})
        fp_branch_stable(i,:) = stable_fp1{i}(1,:);
    end
     if ~isempty(unstable_fp1{i})
         if size(unstable_fp1{i},1)==1
             fp_branch_zeros(i,:) = unstable_fp1{i}(1,:);
         elseif size(unstable_fp1{i},1)==2
             fp_branch_zeros(i,:) = unstable_fp1{i}(1,:);
             fp_branch_unstable(i,:) = unstable_fp1{i}(2,:);             
         end
    end

    if ~isempty(stable_fp2{i})
        if size(stable_fp2{i},1) == 2
        lower_branch_stable(i,:) = stable_fp2{i}(1,:);
        upper_branch_stable(i,:) = stable_fp2{i}(2,:);
        else
            warning('Unexpected number of rows in cell %d',i);
        end
    end

    if ~isempty(unstable_fp2{i})

        if size(unstable_fp2{i},1) ==2
            if all(unstable_fp2{i}>0)
                lower_branch_unstable_async(i,:) = unstable_fp2{i}(1,:);
                upper_branch_unstable_async(i,:) = unstable_fp2{i}(2,:);
            else
                lower_branch_unstable_sync(i,:) = unstable_fp2{i}(1,:);
                upper_branch_unstable_sync(i,:) = unstable_fp2{i}(2,:);
            end

        elseif size(unstable_fp2{i},1) == 4
            lower_branch_unstable_sync(i,:) = unstable_fp2{i}(1,:);
            lower_branch_unstable_async(i,:) = unstable_fp2{i}(2,:);
            upper_branch_unstable_async(i,:) = unstable_fp2{i}(3,:);
            upper_branch_unstable_sync(i,:) = unstable_fp2{i}(4,:);

        else
            warning('Unexpected number of rows in cell %d',i);
        end
    end
end

%% Plotting
% Linear scale
varNames = {'Larvae (larvae/ha)','Nymphs (nymphs/ha)','Adults (adults/ha)'};
for col = 1:3        % 1=L, 2=N, 3=A

figure; hold on;
xlabel('$c$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel(varNames{col},'FontSize',14);
hUnstable = plot(nan, nan, '--r', 'LineWidth', 1.5);
hStable = plot(nan, nan, 'b', 'LineWidth', 1.5);
hold on

y1 = lower_branch_unstable_sync(:,col);
mask = ~isnan(y1);
plot(c_values(mask), y1(mask),'r--','LineWidth',3)
plot(c_values(mask), y1(mask),'w','LineWidth',1)

y1 = upper_branch_unstable_sync(:,col);
mask = ~isnan(y1);
plot(c_values(mask), y1(mask),'r--','LineWidth',3)
plot(c_values(mask), y1(mask),'w','LineWidth',1)

y1 = lower_branch_unstable_async(:,col);
mask = ~isnan(y1);
plot(c_values(mask), y1(mask),'r--','LineWidth',3)
plot(c_values(mask), y1(mask),'w','LineWidth',1)

y1 = upper_branch_unstable_async(:,col);
mask = ~isnan(y1);
plot(c_values(mask), y1(mask),'r--','LineWidth',3)
plot(c_values(mask), y1(mask),'w','LineWidth',1)

y1 = lower_branch_stable(:,col);
mask = ~isnan(y1);
plot(c_values(mask), y1(mask),'b','LineWidth',3)
plot(c_values(mask), y1(mask),'w','LineWidth',1)

y1 = upper_branch_stable(:,col);
mask = ~isnan(y1);
plot(c_values(mask), y1(mask),'b','LineWidth',3)
plot(c_values(mask), y1(mask),'w','LineWidth',1)

plot(c_values, fp_branch_stable(:,col),'b','LineWidth',1.5);
plot(c_values, fp_branch_unstable(:,col),'r--','LineWidth',1.5);
%plot(c_values, fp_branch_zeros(:,col),'r--','LineWidth',1.5);

legend([hStable, hUnstable], {'Stable', 'Unstable'}, 'FontSize',14,'Location','best');
end

% asinh scale
a = 1e4;
varNames = {'Larvae (larvae/ha)','Nymphs (nymphs/ha)','Adults (adults/ha)'};
for col = 1:3        % 1=L, 2=N, 3=A

figure; hold on;
xlabel('$c$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel(varNames{col},'FontSize',14);
hUnstable = plot(nan, nan, '--r', 'LineWidth', 1.5);
hStable = plot(nan, nan, 'b', 'LineWidth', 1.5);
hold on

plot(c_values, asinh(fp_branch_zeros(:,col)/a),'r--','LineWidth',1.5);

y1 = lower_branch_unstable_sync(:,col);
mask = ~isnan(y1);
plot(c_values(mask), asinh(y1(mask)/a),'r--','LineWidth',3)
plot(c_values(mask), asinh(y1(mask)/a),'w','LineWidth',1)

y1 = upper_branch_unstable_sync(:,col);
mask = ~isnan(y1);
plot(c_values(mask), asinh(y1(mask)/a),'r--','LineWidth',3)
plot(c_values(mask), asinh(y1(mask)/a),'w','LineWidth',1)

y1 = lower_branch_unstable_async(:,col);
mask = ~isnan(y1);
plot(c_values(mask), asinh(y1(mask)/a),'r--','LineWidth',3)
plot(c_values(mask), asinh(y1(mask)/a),'w','LineWidth',1)

y1 = upper_branch_unstable_async(:,col);
mask = ~isnan(y1);
plot(c_values(mask), asinh(y1(mask)/a),'r--','LineWidth',3)
plot(c_values(mask), asinh(y1(mask)/a),'w','LineWidth',1)

y1 = lower_branch_stable(:,col);
mask = ~isnan(y1);
plot(c_values(mask), asinh(y1(mask)/a),'b','LineWidth',3)
plot(c_values(mask), asinh(y1(mask)/a),'w','LineWidth',1)

y1 = upper_branch_stable(:,col);
mask = ~isnan(y1);
plot(c_values(mask), asinh(y1(mask)/a),'b','LineWidth',3)
plot(c_values(mask), asinh(y1(mask)/a),'w','LineWidth',1)

plot(c_values, asinh(fp_branch_stable(:,col)/a),'b','LineWidth',1.5);
plot(c_values, asinh(fp_branch_unstable(:,col)/a),'r--','LineWidth',1.5);

legend([hStable, hUnstable], {'Stable', 'Unstable'}, 'FontSize',14,'Location','best');
end

% 3D
figure; hold on; grid on;
xlabel('$c$', 'Interpreter','latex', 'FontSize', 14);
ylabel('Larvae (larvae/ha)', 'FontSize', 14);
zlabel('Adults (adults/ha)', 'FontSize', 14);
view(3);
hStable3 = plot3(nan,nan,nan,'b','LineWidth',1.5);
hUnstable3 = plot3(nan,nan,nan,'--r','LineWidth',1.5);

plot3(c_values, fp_branch_zeros(:,1),fp_branch_zeros(:,3),'r--','LineWidth',1.5);

y1 = lower_branch_unstable_sync(:,1);
y2 = lower_branch_unstable_sync(:,3);
mask = ~isnan(y1);
plot3(c_values(mask), y1(mask),y2(mask),'r--','LineWidth',3)
plot3(c_values(mask), y1(mask),y2(mask),'w','LineWidth',1)

y1 = upper_branch_unstable_sync(:,1);
y2 = upper_branch_unstable_sync(:,3);
mask = ~isnan(y1);
plot3(c_values(mask), y1(mask),y2(mask),'r--','LineWidth',3)
plot3(c_values(mask), y1(mask),y2(mask),'w','LineWidth',1)

y1 = lower_branch_unstable_async(:,1);
y2 = lower_branch_unstable_async(:,3);
mask = ~isnan(y1);
plot3(c_values(mask), y1(mask),y2(mask),'r--','LineWidth',3)
plot3(c_values(mask), y1(mask),y2(mask),'w','LineWidth',1)

y1 = upper_branch_unstable_async(:,1);
y2 = upper_branch_unstable_async(:,3);
mask = ~isnan(y1);
plot3(c_values(mask), y1(mask),y2(mask),'r--','LineWidth',3)
plot3(c_values(mask), y1(mask),y2(mask),'w','LineWidth',1)

y1 = lower_branch_stable(:,1);
y2 = lower_branch_stable(:,3);
mask = ~isnan(y1);
plot3(c_values(mask), y1(mask),y2(mask),'b','LineWidth',3)
plot3(c_values(mask), y1(mask),y2(mask),'w','LineWidth',1)

y1 = upper_branch_stable(:,1);
y2 = upper_branch_stable(:,3);
mask = ~isnan(y1);
plot3(c_values(mask), y1(mask),y2(mask),'b','LineWidth',3)
plot3(c_values(mask), y1(mask),y2(mask),'w','LineWidth',1)

plot3(c_values, fp_branch_stable(:,1),fp_branch_stable(:,3),'b','LineWidth',1.5);
plot3(c_values, fp_branch_unstable(:,1),fp_branch_unstable(:,3),'r--','LineWidth',1.5);

legend([hStable3 hUnstable3], {'Stable','Unstable'}, 'FontSize', 14, 'Location', 'best');

% 3D
figure; hold on; grid on;
xlabel('$c$', 'Interpreter','latex', 'FontSize', 14);
ylabel('Larvae (larvae/ha)', 'FontSize', 14);
zlabel('Nymphs (nymphs/ha)', 'FontSize', 14);
view(3);
hStable3 = plot3(nan,nan,nan,'b','LineWidth',1.5);
hUnstable3 = plot3(nan,nan,nan,'--r','LineWidth',1.5);

plot3(c_values, fp_branch_zeros(:,1),fp_branch_zeros(:,2),'r--','LineWidth',1.5);

y1 = lower_branch_unstable_sync(:,1);
y2 = lower_branch_unstable_sync(:,2);
mask = ~isnan(y1);
plot3(c_values(mask), y1(mask),y2(mask),'r--','LineWidth',3)
plot3(c_values(mask), y1(mask),y2(mask),'w','LineWidth',1)

y1 = upper_branch_unstable_sync(:,1);
y2 = upper_branch_unstable_sync(:,2);
mask = ~isnan(y1);
plot3(c_values(mask), y1(mask),y2(mask),'r--','LineWidth',3)
plot3(c_values(mask), y1(mask),y2(mask),'w','LineWidth',1)

y1 = lower_branch_unstable_async(:,1);
y2 = lower_branch_unstable_async(:,2);
mask = ~isnan(y1);
plot3(c_values(mask), y1(mask),y2(mask),'r--','LineWidth',3)
plot3(c_values(mask), y1(mask),y2(mask),'w','LineWidth',1)

y1 = upper_branch_unstable_async(:,1);
y2 = upper_branch_unstable_async(:,2);
mask = ~isnan(y1);
plot3(c_values(mask), y1(mask),y2(mask),'r--','LineWidth',3)
plot3(c_values(mask), y1(mask),y2(mask),'w','LineWidth',1)

y1 = lower_branch_stable(:,1);
y2 = lower_branch_stable(:,2);
mask = ~isnan(y1);
plot3(c_values(mask), y1(mask),y2(mask),'b','LineWidth',3)
plot3(c_values(mask), y1(mask),y2(mask),'w','LineWidth',1)

y1 = upper_branch_stable(:,1);
y2 = upper_branch_stable(:,2);
mask = ~isnan(y1);
plot3(c_values(mask), y1(mask),y2(mask),'b','LineWidth',3)
plot3(c_values(mask), y1(mask),y2(mask),'w','LineWidth',1)

plot3(c_values, fp_branch_stable(:,1),fp_branch_stable(:,2),'b','LineWidth',1.5);
plot3(c_values, fp_branch_unstable(:,1),fp_branch_unstable(:,2),'r--','LineWidth',1.5);

legend([hStable3 hUnstable3], {'Stable','Unstable'}, 'FontSize', 14, 'Location', 'best');

%% To correct labels on asinh scale figures
%% Run this after running (1) DRIVER and (2) plotter with the figure needing to be fixed open
% % for c_values
% ax = gca;
% % %% Larvae
% % yt=[0 0.5320 1.0639 1.5959 2.1278 2.6598 3.1917 3.7237 4.2556 4.7876];
% % ylabels={'0', '0.06', '0.1','0.2','0.4', '0.7','1.2','2.1','3.5','6'};
% 
% % %% Nymphs
% % yt=[0 0.4550 0.9099 1.3649 1.8198 2.2748 2.7297 3.1847 3.6396 4.0946];
% % ylabels={'0', '0.05', '0.1','0.2','0.3', '0.5','0.8','1.2','1.9','3'};
% 
% %% Adults
% yt=[0 0.4347 0.894 1.3041 1.7388 2.1736 2.6083 3.0430 3.4777 3.9124];
% ylabels={'0', '0.05', '0.1','0.2','0.3', '0.4','0.7','1','1.6','2.5'};
% 
% yticks(yt)
% yticklabels([])
% 
% xlim_vals = xlim;
% 
% for i = 1:length(yt)
%     y_pos = yt(i);
% 
%     text(xlim_vals(1) - 0.015*(xlim_vals(2)-xlim_vals(1)), y_pos, ylabels{i}, ...
%         'HorizontalAlignment','right', ...
%         'VerticalAlignment','middle', ...
%         'Clipping','off');
% end
% 
% text(0, 1.04, '\times 10^5', 'Units','normalized')
% 
% % ylh = ylabel('Larvae (larvae/ha)', 'FontSize', 14);
% % ylh = ylabel('Nymphs (nymphs/ha)', 'FontSize', 14);
% ylh = ylabel('Adults (adults/ha)', 'FontSize', 14);
% 
% ylh.Units = 'normalized';
% ylh.Position(1) = -0.08;
% 
% % ylim([0 5])
% % ylim([0 4.5])
% ylim([0 4])

% % for r_values
% ax = gca;
% 
% % --- X data ---
% xt = [0 0.17074 0.19706 0.5 1 1.1466 1.5 2]*1e5;
% xlabels = {'0','B_2','B_3','0.5','1','B_4','1.5','2'};
% xticks(xt)
% 
% % --- Y data ---
% % % Larvae
% % yt = [0 0.7222    1.4444    2.1667    2.8889    3.6111    4.3333    5.0556    5.7778 6.5];
% % ylabels = {'0', '0.08', '0.2','0.4','0.9', '1.8','3.8','7.8','16.1','33.3'};
% % yticks(yt)
% %Nymphs&Adults
% yt = [0 0.6111 1.2222 1.8333 2.4444 3.0556 3.6667 4.2778 4.8889 5.5];
% ylabels = {'0', '0.06', '0.1', '0.3', '0.6','1','2','3.6','6.6','12.2'};
% yticks(yt)
% 
% % Turn off default labels
% ax.XTickLabel = [];
% ax.YTickLabel = [];
% 
% % --- Custom X labels ---
% for i = 1:length(xt)
%     x_pos = xt(i);
% 
%     if strcmp(xlabels{i}, 'B_2')
%         x_pos = x_pos - 0.025e5;
%     elseif strcmp(xlabels{i}, 'B_3')
%         x_pos = x_pos + 0.025e5;
%     end
% 
%     text(x_pos, ax.YLim(1)-0.015*(ax.YLim(2)-ax.YLim(1)), xlabels{i}, ...
%         'HorizontalAlignment','center', ...
%         'VerticalAlignment','top');
% end
% 
% % --- Custom Y labels ---
% for i = 1:length(yt)
%     y_pos = yt(i);
% 
%     text(ax.XLim(1)-0.015*(ax.XLim(2)-ax.XLim(1)), y_pos, ylabels{i}, ...
%         'HorizontalAlignment','right', ...
%         'VerticalAlignment','middle');
% end
% 
% % --- Exponent notes ---
% text(0.93, -0.09, '\times 10^5', 'Units','normalized') %x
% text(0, 1.04, '\times 10^5', 'Units','normalized') %y
% 
% % --- X label ---
% xlh = xlabel('$\mathcal{R}_d$', 'FontSize', 14, 'Interpreter', 'latex');
% xlh.Position(2) = xlh.Position(2) - 0.265;
% 
% % --- Y label ---
% % ylh = ylabel('Larvae (larvae/ha)', 'FontSize', 14);
% % ylh = ylabel('Nymphs (nymphs/ha)', 'FontSize', 14);
% ylh = ylabel('Adults (adults/ha)', 'FontSize', 14);
% 
% ylh.Units = 'normalized';   % key step
% ylh.Position(1) = ylh.Position(1) - 0.055;
% 
% ylim([0 5.5])
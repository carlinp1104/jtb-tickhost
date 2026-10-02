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
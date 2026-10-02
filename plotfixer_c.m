%% Run this after running (1) DRIVER and (2) plotter

ax = gca;
% %% Larvae
% yt=[0 0.5320 1.0639 1.5959 2.1278 2.6598 3.1917 3.7237 4.2556 4.7876];
% ylabels={'0', '0.06', '0.1','0.2','0.4', '0.7','1.2','2.1','3.5','6'};

% %% Nymphs
% yt=[0 0.4550 0.9099 1.3649 1.8198 2.2748 2.7297 3.1847 3.6396 4.0946];
% ylabels={'0', '0.05', '0.1','0.2','0.3', '0.5','0.8','1.2','1.9','3'};

%% Adults
yt=[0 0.4347 0.894 1.3041 1.7388 2.1736 2.6083 3.0430 3.4777 3.9124];
ylabels={'0', '0.05', '0.1','0.2','0.3', '0.4','0.7','1','1.6','2.5'};

yticks(yt)
yticklabels([])

xlim_vals = xlim;

for i = 1:length(yt)
    y_pos = yt(i);

    text(xlim_vals(1) - 0.015*(xlim_vals(2)-xlim_vals(1)), y_pos, ylabels{i}, ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','middle', ...
        'Clipping','off');
end

text(0, 1.04, '\times 10^5', 'Units','normalized')

% ylh = ylabel('Larvae (larvae/ha)', 'FontSize', 14);
% ylh = ylabel('Nymphs (nymphs/ha)', 'FontSize', 14);
ylh = ylabel('Adults (adults/ha)', 'FontSize', 14);

ylh.Units = 'normalized';
ylh.Position(1) = -0.08;

% ylim([0 5])
% ylim([0 4.5])
ylim([0 4])
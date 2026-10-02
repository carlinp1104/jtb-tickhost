%% Run this after running (1) DRIVER and (2) plotter

ax = gca;

% --- X data ---
xt = [0 0.17074 0.19706 0.5 1 1.1466 1.5 2]*1e5;
xlabels = {'0','B_2','B_3','0.5','1','B_4','1.5','2'};
xticks(xt)

% --- Y data ---
% % Larvae
% yt = [0 0.7222    1.4444    2.1667    2.8889    3.6111    4.3333    5.0556    5.7778 6.5];
% ylabels = {'0', '0.08', '0.2','0.4','0.9', '1.8','3.8','7.8','16.1','33.3'};
% yticks(yt)
%Nymphs&Adults
yt = [0 0.6111 1.2222 1.8333 2.4444 3.0556 3.6667 4.2778 4.8889 5.5];
ylabels = {'0', '0.06', '0.1', '0.3', '0.6','1','2','3.6','6.6','12.2'};
yticks(yt)

% Turn off default labels
ax.XTickLabel = [];
ax.YTickLabel = [];

% --- Custom X labels ---
for i = 1:length(xt)
    x_pos = xt(i);

    if strcmp(xlabels{i}, 'B_2')
        x_pos = x_pos - 0.025e5;
    elseif strcmp(xlabels{i}, 'B_3')
        x_pos = x_pos + 0.025e5;
    end

    text(x_pos, ax.YLim(1)-0.015*(ax.YLim(2)-ax.YLim(1)), xlabels{i}, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','top');
end

% --- Custom Y labels ---
for i = 1:length(yt)
    y_pos = yt(i);

    text(ax.XLim(1)-0.015*(ax.XLim(2)-ax.XLim(1)), y_pos, ylabels{i}, ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','middle');
end

% --- Exponent notes ---
text(0.93, -0.09, '\times 10^5', 'Units','normalized') %x
text(0, 1.04, '\times 10^5', 'Units','normalized') %y

% --- X label ---
xlh = xlabel('$\mathcal{R}_d$', 'FontSize', 14, 'Interpreter', 'latex');
xlh.Position(2) = xlh.Position(2) - 0.265;

% --- Y label ---
% ylh = ylabel('Larvae (larvae/ha)', 'FontSize', 14);
% ylh = ylabel('Nymphs (nymphs/ha)', 'FontSize', 14);
ylh = ylabel('Adults (adults/ha)', 'FontSize', 14);

ylh.Units = 'normalized';   % key step
ylh.Position(1) = ylh.Position(1) - 0.055;

ylim([0 5.5])
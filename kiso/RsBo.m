close all;
Rs = readmatrix("Rs.xlsx");
xRs = Rs(:,1);
yRs = Rs(:,2);

Bo = readmatrix("Bo.xlsx");
xBo = Bo(:,1);
yBo = Bo(:,2);

figure;
axRs = subplot(1,2,1);
axBo = subplot(1,2,2);



hold(axRs, "on");
plot(axRs, xRs, yRs, 'LineWidth', 2);
xlabel(axRs, 'Pressure [psia]');
ylabel(axRs, 'Rs [scf/bbl]');
ylim(axRs, [0, 2500]);
xlim(axRs, [0, 9000]);
xticks(axRs, 0:1000:9000);
grid(axRs, "on");
axRs.FontSize = 20;
hold off;

hold(axBo, "on");
plot(axBo, xBo, yBo, 'LineWidth', 2);
xlabel(axBo, 'Pressure [psia]');
ylabel(axBo, 'Bo, RC/SC');
ylim(axBo, [0.8, 2.4]);
xlim(axBo, [0, 9000]);
xticks(axBo, 0:1000:9000);
grid(axBo, "on");
axBo.FontSize = 20;
hold off;
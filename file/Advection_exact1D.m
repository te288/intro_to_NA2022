%% Input Parameters
c = 1;
t = 0;
x = linspace(-3,10, 200);
ax = subplot(1,1,1);
plot(ax, x, exp(-(x-c*t).^2)); % Initial plot
xlim(ax, [-3,10]);ylim(ax, [-0.1,1.1]);
name = title(['t = ',num2str(t),'[s]']);
hold(ax ,"on")


for t = 0:1:15 % time data
    name = title(['t = ',num2str(t,'%04.1f'),'[s]']);
    plot(ax, x, exp(-(x-c*t).^2));
    hold(ax ,"on")
    % filename = sprintf('plotEXACT%04.1f.png', t); % file name 
    %saveas(fig, filename); % save figure as png
end
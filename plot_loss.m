function plot_loss(losses, dt, nt, imin, icrit)
    figure;
    min_ = min(losses);
    max_ = max(losses);
    hold on
    plot((1:1:nt)*dt/1000, losses, 'LineWidth',2)
    patch([imin imin icrit icrit]*dt/1000, [min_ max_ max_ min_], [0.5,0.5,0.5], 'FaceAlpha', 0.4, 'EdgeColor', 'none')
    hold off
    legend('Training', 'Error','Location', 'best')
    xlabel('Time (s)')
    ylabel('Error $e(t)$', 'Interpreter', 'LaTeX')
    title('Error curve')
    xlim([0,nt*dt/1000]);
    set(gca, 'YScale', 'log')
end
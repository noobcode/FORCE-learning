function plot_loss(losses, dt, nt, imin, icrit)
    figure;
    min_ = min(losses);
    max_ = max(losses);
    
    %patch([imin imin icrit icrit]*dt/1000, [min_ max_ max_ min_], [0.5,0.5,0.5], 'FaceAlpha', 0.4, 'EdgeColor', 'none')
    hold on
    plot((1:1:nt)*dt/1000, losses, 'LineWidth',1)
    xlabel('Time (s)')
    ylabel('Error $e(t)$', 'Interpreter', 'LaTeX')
    %legend('Training', 'Error','Location', 'best')
    legend('Error','Location', 'best')
    title('Error Curve')
    xlim([0,nt*dt/1000]);
    ylim([min_, max_]);
    set(gca, 'YScale', 'log')
    hold off
end
function plot_approximant_vs_target(zx, current, dt, nt, imin, icrit)
        figure;
        min_ = min(min(zx),min(current));
        max_ = max(max(zx),max(current));
        patch([imin imin icrit icrit]*dt/1000, [min_ max_ max_ min_], [0.5,0.5,0.5], 'FaceAlpha', 0.4, 'EdgeColor', 'none')
        hold on
        plot(dt*(1:1:nt)/1000, zx,'k','LineWidth',2)
        plot(dt*(1:1:nt)/1000, current,'b--','LineWidth',2)
        xlabel('Time (s)')
        ylabel('$\hat{x}(t)$','Interpreter','LaTeX')
        legend('Training', 'Target Signal', 'Approximant','Location', 'best')
        title('Target vs Approximant')
        xlim([0,nt*dt/1000]);
        hold off
end
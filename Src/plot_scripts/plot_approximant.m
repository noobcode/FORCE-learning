function plot_approximant(current, dt, nt, imin, icrit, icrit_2)
    figure;
    
    min_ = min(current);
    max_ = max(current);
    
    patch([imin imin icrit icrit]*dt/1000, [min_ max_ max_ min_], [0.5,0.5,0.5], 'FaceAlpha', 0.4, 'EdgeColor', 'none')
    hold on
    
    if icrit_2 ~= -1
        patch([icrit_2 icrit_2 nt nt]*dt/1000, [min_ max_ max_ min_], [0.9,0.1,0.1], 'FaceAlpha', 0.4, 'EdgeColor', 'none')
    end
    
    plot(dt*(1:1:nt)/1000, current,'b--','LineWidth',2)
    xlabel('Time (s)')
    ylabel('$\hat{x}(t)$','Interpreter','LaTeX')
    title('Approximant')
    xlim([0,nt*dt/1000]);
    ylim([min_, max_]);
    
    if icrit_2 ~= -1
        legend('Training', 'Post Control', 'Approximant','Location', 'best')
    else
        legend('Training', 'Approximant','Location', 'best')
    end
    
    hold off

end
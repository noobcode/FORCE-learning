function plot_approximant_vs_target(zx, current, dt, nt, imin, icrit, icrit_2)
    figure;
    min_ = min(min(zx),min(current));
    max_ = max(max(zx),max(current));
    
    patch([imin imin icrit icrit]*dt/1000, [min_ max_ max_ min_], [0.5,0.5,0.5], 'FaceAlpha', 0.4, 'EdgeColor', 'none')
    hold on
    
    if icrit_2 ~= -1
        patch([icrit_2 icrit_2 nt nt]*dt/1000, [min_ max_ max_ min_], [0.9,0.1,0.1], 'FaceAlpha', 0.4, 'EdgeColor', 'none')
    end
    
    % pad target signal with itself
    if length(zx) < nt
        times = round(nt/length(zx));
        temp = zx;
        for i=1:times-1
           zx = [zx,temp];
        end
    end
    
    plot(dt*(1:1:nt)/1000, zx,'k','LineWidth',2)
    plot(dt*(1:1:nt)/1000, current,'b--','LineWidth',2)
    xlabel('Time (s)')
    ylabel('$\hat{x}(t)$','Interpreter','LaTeX')
    title('Target vs Approximant')
    xlim([0,nt*dt/1000]);
    ylim([min_, max_]);
    
    if icrit_2 ~= -1
        legend('Training', 'Post Control', 'Target Signal', 'Approximant','Location', 'best')
    else
        legend('Training', 'Target Signal', 'Approximant','Location', 'best')
    end
    
    hold off
end
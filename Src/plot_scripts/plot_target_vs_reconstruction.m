function plot_target_vs_reconstruction(original_zx, reconstructed_zx, dt, nt, imin, icrit)
    figure;
    min_ = min(min(original_zx), min(reconstructed_zx));
    max_ = max(max(original_zx), max(reconstructed_zx));
    hold on
    patch([imin imin icrit icrit]*dt/1000, [min_ max_ max_ min_], [0.5,0.5,0.5], 'FaceAlpha', 0.4, 'EdgeColor', 'none')
    plot(dt*(1:1:nt)/1000, original_zx,'k','LineWidth',2)
    plot(dt*(1:1:nt)/1000, reconstructed_zx,'b--','LineWidth',2)
    xlabel('Time (s)')
    ylabel('$\hat{x}(t)$','Interpreter','LaTeX')
    legend('Training', 'Target Signal', 'Reconstructed approximant','Location', 'best')
    xlim([0,nt*dt/1000]);
    hold off
end
function [] = plot_multidim_approx_vs_target(zx, current, dt, nt, imin, icrit, icrit_2, target_length)
    %% padding
    padding_done = false;
    [n_joints, n_samples] = size(zx);
    
    temp = zx;
    % pad target signal with itself
    if n_samples < nt
        padding_done = true;
        times = round(nt/n_samples);
        for i=1:times-1
           zx = [zx,temp];
        end
    end

    tt = (1:1:nt)*dt;

    %% one plot
    figure; 
    hold on
    plot(tt, zx,'-','LineWidth',1)
    plot(tt, current','--','LineWidth',1)
    hold off
    
    %% vertical sub-plots
    target_length = sum(target_length);
    len = round(target_length/dt);
    
    figure;
    for j = 1:1:n_joints
        max_y = max(max(zx(j,:)), max(current(j,end-len+1:end)));  % lascio cosi?
        min_y = min(min(zx(j,:)), min(current(j,end-len+1:end)));
        
        subplot(n_joints,1,j)
        hold on
        if imin ~= -1 && icrit ~= -1
            patch([imin imin icrit icrit]*dt, [min_y max_y max_y min_y], [0.4,0.4,0.4], 'FaceAlpha', 0.3, 'EdgeColor', 'none')
        end
        if icrit_2 ~= -1
            patch([icrit_2 icrit_2 nt nt]*dt, [min_y max_y max_y min_y], [0.1,0.5,0.1], 'FaceAlpha', 0.3, 'EdgeColor', 'none')
        end
        plot(tt, zx(j,:),'k','LineWidth',1)
        plot(tt, current(j,:),'r--','LineWidth',1)
        hold off
        
        if j < n_joints
            set(gca,'XTick',[]);
        end
        
        set(gca,'YTick',[]);
        ylabel("\theta_" + num2str(j))
        
        ylim([min_y,max_y])
        %ylim([-2,2])
        xlim([1,nt]*dt)
        
        if j == n_joints
            if icrit_2 ~= -1
                legend({'training', 'control', 'target', 'approx'}, 'Orientation','horizontal')
            else
                legend({'training', 'target', 'approx'}, 'Orientation','horizontal')
            end
        end
        
    end

    xlabel('Time (ms)');
    set(gcf,'Position',[50 50 600 1000]);
    suplabel('Joint Trajectories', 't');
    [a, h1] = suplabel('Joint Angle (rad)', 'y');
    set(h1);

    %% grid sub-plot
    % funziona solo se si è fatto il padding
    if padding_done == true
    target_length = sum(target_length);
    
    len = round(target_length/dt);
    %last_zx = zx(:, end-len+1:end);
    last_current = current(:, end-len+1:end);
    tt_1_len = (1:1:len)*dt/1000;
    
    figure;
    for j = 1:1:n_joints
        max_y = max(max(temp(j,:)), max(last_current(j,:)));
        min_y = min(min(temp(j,:)), min(last_current(j,:)));
        
        subplot(2,4,j)
        hold on
        plot(tt_1_len, temp(j,:), 'k','LineWidth',1)
        plot(tt_1_len, last_current(j,:), 'r--','LineWidth',1)
        hold off
        ylim([min_y, max_y])
        xlim([1, len]*dt/1000)
        
        if j < 5
            set(gca,'XTick',[]);
        end
        
        if j == 7
            legend('target', 'approx')
        end
        title("Joint " + num2str(j)) 
    end
    set(gcf,'Position',[100 100 1000 550]);
    [a, h1] = suplabel('Time (s)', 'x');
    [a, h2] = suplabel('Joint Angle (rad)', 'y');
    set(h1); set(h2);
    end
    
end
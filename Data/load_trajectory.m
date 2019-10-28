function [new_positions, positions, time_stamps, target_length] = load_trajectory(file_name, dt, do_plot)
    
    data_table = readtable(file_name);
    
    positions = data_table(1:end,2:8);
    time_stamps = data_table(1:end,end);
    
    % convert to array
    positions = positions{:,:};
    time_stamps = time_stamps{:,:};
    
    target_length = time_stamps(end); % milliseconds
    tt = (1:1:round(target_length/dt))*dt;
    
    [n_samples, n_joints] = size(positions);
    
    n_samples
    
    new_positions = zeros(n_joints, length(tt));
    
    for i = 1:n_joints
        p = polyfit(time_stamps, positions(:,i), 3);
        new_positions(i,:) = polyval(p, tt);
    end

    %% plot
    
    
    if do_plot
        figure;
        for i = 1:n_joints
            plot(tt, new_positions(i,:), '-o'), hold on
        end
        hold off

        figure
        plot(time_stamps, positions, '-o')
        legend('\theta_1','\theta_2','\theta_3','\theta_4','\theta_5','\theta_6','\theta_7','Location', 'best')
        xlabel('Time (ms)')
        ylabel('Joint Angle (radians)')
        title('Joints Trajectory')
        xlim([0, target_length])
        hold off
    end
end
function [zx1_new, zx2_new, zx3_new, zx4_new, target_length] = load_four_trajectories_same_length(dt)
    [~, zx1, t_stamps1, T1] = load_trajectory('trajectory_1.csv', dt, false);
    [~, zx2, t_stamps2, T2] = load_trajectory('trajectory_2.csv', dt, false);
    [~, zx3, t_stamps3, T3] = load_trajectory('trajectory_3.csv', dt, false);
    [~, zx4, t_stamps4, T4] = load_trajectory('trajectory_4.csv', dt, false);
    
    target_length = max([T1 T2 T3 T4]);
    
    tt = (1:1:round(target_length/dt))*dt;

    [~, n_joints] = size(zx1);

    zx1_new = zeros(n_joints, length(tt));  zx2_new = zeros(n_joints, length(tt));
    zx3_new = zeros(n_joints, length(tt));  zx4_new = zeros(n_joints, length(tt));

    for i = 1:n_joints
        p1 = polyfit(t_stamps1, zx1(:,i), 3);
        p2 = polyfit(t_stamps2, zx2(:,i), 3);
        p3 = polyfit(t_stamps3, zx3(:,i), 3);
        p4 = polyfit(t_stamps4, zx4(:,i), 3);

        zx1_new(i,:) = polyval(p1, tt);
        zx2_new(i,:) = polyval(p2, tt);
        zx3_new(i,:) = polyval(p3, tt);
        zx4_new(i,:) = polyval(p4, tt);
    end

    %{
    figure;
    plot(tt, zx1_new);
    figure;
    plot(tt, zx2_new);
    figure;
    plot(tt, zx3_new);
    figure;
    plot(tt, zx4_new);

    size(zx1_new), size(zx2_new),size(zx3_new),size(zx4_new)
    %}
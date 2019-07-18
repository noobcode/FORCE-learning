clear all
close all

d_values = [-50 0 50 100 150 200];
legend_values = ["0.3 Hz" "0.5 Hz" "0.8 Hz" "1 Hz" "5 Hz" "7 Hz" "12 Hz"];

% mean error Hz
error_0_3 = [0.0035 0.8354 1.0363 1.0002 1.3219 0.9555];
error_0_5 = [0.0005 0.9633 1.1828 1.1485 1.1839 1.4324];
error_0_8 = [0.0004 0.0314 0.8236 0.2271 0.6201 1.0904];
error_1 =   [0.0002 0.0042 0.1604 0.0707 0.2803 0.8043];
error_5 =   [0.0008 0.0024 0.0031 0.0080 0.0568 0.2092];
error_7 =   [0.0031 0.0029 0.0017 0.0159 0.0922 0.3423];
error_12 =  [0.0274 0.0110 0.0077 0.0199 0.1143 0.2831];

% std error Hz
error_std_0_3 = [0.0022 0.2985 0.3147 0.3012 0.3440 0.1196];
error_std_0_5 = [0.0004 0.6157 1.0435 0.2579 0.2102 0.3537];
error_std_0_8 = [0.0007 0.0383 0.3235 0.1961 0.3683 0.4142];
error_std_1 =   [0.0002 0.0050 0.2359 0.0513 0.2162 0.3749];
error_std_5 =   [0.0011 0.0028 0.0042 0.0130 0.0450 0.0986];
error_std_7 =   [0.0031 0.0030 0.0011 0.0121 0.0546 0.2388];
error_std_12 =  [0.0175 0.0148 0.0069 0.0242 0.1300 0.2311];

% mean firing rate Hz
firing_rate_0_3 = [127.4213 66.9317 49.9101 40.6669 35.3160 10.8206];
firing_rate_0_5 = [126.7475 78.6007 55.8020 42.9709 35.7839 28.7686];
firing_rate_0_8 = [124.3889 75.4338 49.2399 41.6946 31.2717 27.1882];
firing_rate_1 =   [122.9335 75.5090 53.8603 42.1748 33.6764 26.4775];
firing_rate_5 =   [111.8739 74.2237 58.6097 49.1684 42.1641 36.7864];
firing_rate_7 =   [110.0842 72.8004 58.0457 49.0444 42.4460 36.9730];
firing_rate_12 =  [107.4061 70.2030 56.3157 48.1619 42.3296 37.7564];

% std firing rate Hz
firing_rate_std_0_3 = [3.6725 5.6981 5.4620 3.1070 8.9538 3.3611];
firing_rate_std_0_5 = [3.4055 16.2739 16.1879 3.2579 1.9988 7.0668];
firing_rate_std_0_8 = [3.2347 2.3926 5.6293 0.7785 1.8470 0.8322];
firing_rate_std_1 =   [3.1984 1.6012 1.2385 0.8422 1.2541 1.2853];
firing_rate_std_5 =   [3.1431 1.6655 1.1252 0.9106 0.5901 0.9040];
firing_rate_std_7 =   [3.1920 1.5903 1.1428 0.9931 0.9463 1.2213];
firing_rate_std_12 =  [0.0175 0.0148 0.0069 0.0242 0.1300 0.2311];

%% plot 
figure
hold on
errorbar(d_values, error_0_3, error_std_0_3, 'LineWidth',1)
errorbar(d_values, error_0_5, error_std_0_5, 'LineWidth',1)
errorbar(d_values, error_0_8, error_std_0_8, 'LineWidth',1)
errorbar(d_values, error_1, error_std_1, 'LineWidth',1)
errorbar(d_values, error_5, error_std_5, 'LineWidth',1)
errorbar(d_values, error_7, error_std_7, 'LineWidth',1)
errorbar(d_values, error_12, error_std_12, 'LineWidth',1)
xlabel('d')
ylabel('average error')
set(gca, 'YScale', 'log')
legend(legend_values, 'Location', 'best')
hold off

figure
hold on
errorbar(d_values, firing_rate_0_3, firing_rate_std_0_3, 'LineWidth',1)
errorbar(d_values, firing_rate_0_5, firing_rate_std_0_5, 'LineWidth',1)
errorbar(d_values, firing_rate_0_8, firing_rate_std_0_8, 'LineWidth',1)
errorbar(d_values, firing_rate_1, firing_rate_std_1, 'LineWidth',1)
errorbar(d_values, firing_rate_5, firing_rate_std_5, 'LineWidth',1)
errorbar(d_values, firing_rate_7, firing_rate_std_5, 'LineWidth',1)
errorbar(d_values, firing_rate_12, firing_rate_std_12, 'LineWidth',1)
xlabel('d')
ylabel('average firing rate')
legend(legend_values, 'Location', 'best')
hold off
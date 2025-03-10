% 清除所有变量和关闭所有图形
clear all;
close all;

% OFDM仿真参数
carrier_count = 200; % 子载波数
symbol_count = 100;  % 总符号数
ifft_length = 512;  % IFFT长度
CS_length = 20;     % 循环后缀长度
alpha = 1.5 / 32;   % 升余弦窗系数
bit_per_symbol = 4; % 调制方式决定
SNR_range = 0:1:20; % 信噪比范围
mult_path_am = [1 0.2 0.1]; % 多径幅度
%符号周期长度由保护间隔和有用数据持续时间长度（符号周期）组成，
%一般选择符号周期长度是保护间隔长度的5倍，即保护间隔的长度为有效符号周期的1/4。
mutt_path_time = [0 20 50]; % 多径时延  均方根值为31.09  保护间隔的时间长度为时延扩展均方根值的2-4倍
CP_lengths = [16, 32, 64, 128, 256]; % 不同的循环前缀长度

% 初始化误码率存储数组
error_rates_diff_CP = cell(length(CP_lengths), 1);

% 对每个CP长度进行仿真
for cp_idx = 1:length(CP_lengths)
    CP_length = CP_lengths(cp_idx);
    error_rates = [];
    
    for snr = SNR_range
        SNR = snr;
        
        % 产生随机序列
        bit_length = carrier_count * symbol_count * bit_per_symbol;
        bit_sequence = round(rand(1, bit_length))'; % 列向量
        
        % 子载波调制
        carrier_position = 29:228;
        conj_position = 485:-1:286;
        bit_moded = qammod(bit_sequence, 16, 'InputType', 'bit');
        
        % 串并转换
        ifft_position = zeros(ifft_length, symbol_count);
        bit_moded = reshape(bit_moded, carrier_count, symbol_count);
        ifft_position(carrier_position, :) = bit_moded;
        ifft_position(conj_position, :) = conj(bit_moded);
        signal_time = ifft(ifft_position, ifft_length);
        
        % 加循环前缀和后缀
        signal_time_C = [signal_time(end-CP_length+1:end, :); signal_time];
        signal_time_C = [signal_time_C; signal_time_C(1:CS_length, :)];

        % 加窗
        signal_window = signal_time_C .* repmat(rcoswindow(alpha, size(signal_time_C, 1)), 1, symbol_count);
        
        % 发送信号，多径信道
        signal_Tx = reshape(signal_window, 1, []);
        signal_origin = reshape(signal_time_C, 1, []);
        path2 = 0.2 * [zeros(1, 20) signal_Tx(1:end-20)];
        path3 = 0.1 * [zeros(1, 50) signal_Tx(1:end-50)];
        signal_Tx_mult = signal_Tx + path2 + path3;
        
        % 加AWGN
        signal_power_sig = var(signal_Tx);
        signal_power_mut = var(signal_Tx_mult);
        SNR_linear = 10^(SNR / 10);
        noise_power_mut = signal_power_mut / SNR_linear;
        noise_power_sig = signal_power_sig / SNR_linear;
        noise_sig = randn(size(signal_Tx)) * sqrt(noise_power_sig);
        noise_mut = randn(size(signal_Tx_mult)) * sqrt(noise_power_mut);
        Rx_data_sig = signal_Tx + noise_sig;
        Rx_data_mut = signal_Tx_mult + noise_mut;
        
        % 串并转换
        Rx_data_mut = reshape(Rx_data_mut, ifft_length + CS_length + CP_length, []);
        Rx_data_sig = reshape(Rx_data_sig, ifft_length + CS_length + CP_length, []);
        
        % 去循环前缀和后缀
        Rx_data_sig(1:CP_length, :) = [];
        Rx_data_sig(end-CS_length+1:end, :) = [];
        Rx_data_mut(1:CP_length, :) = [];
        Rx_data_mut(end-CS_length+1:end, :) = [];
        
        % FFT
        fft_sig = fft(Rx_data_sig);
        fft_mut = fft(Rx_data_mut);
        
        % 降采样
        data_sig = fft_sig(carrier_position, :);
        data_mut = fft_mut(carrier_position, :);
        
        % 逆映射
        bit_demod_sig = reshape(qamdemod(data_sig, 16, 'OutputType', 'bit'), [], 1);
        bit_demod_mut = reshape(qamdemod(data_mut, 16, 'OutputType', 'bit'), [], 1);
        
        % 计算误码率
        error_bit_sig = sum(bit_demod_sig ~= bit_sequence);
        error_bit_mut = sum(bit_demod_mut ~= bit_sequence);
        error_rate_sig = error_bit_sig / bit_length;
        error_rate_mut = error_bit_mut / bit_length;
        error_rates = [error_rates; error_rate_sig, error_rate_mut];
    end
    
    % 存储当前CP长度下的误码率
    error_rates_diff_CP{cp_idx} = error_rates;
end

% 绘制不同CP长度下的误码率与信噪比曲线
figure;
hold on;
colors = lines(length(CP_lengths)); % 选择不同的颜色来区分不同的CP长度
for i = 1:length(CP_lengths)
    % 绘制单径情况下的BER
    plot(SNR_range, error_rates_diff_CP{i}(:, 1), 'Color', colors(i, :), 'LineWidth', 1.5, 'DisplayName', ['CP = ', num2str(CP_lengths(i))]);
    
    % 绘制多径情况下的BER
    plot(SNR_range, error_rates_diff_CP{i}(:, 2), 'Color', colors(i, :), 'LineStyle', '--', 'LineWidth', 1.5, 'DisplayName', ['CP = ', num2str(CP_lengths(i)), ' (Multipath)']);
end
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs SNR for Different CP Lengths');
legend show;
grid on;
hold off;

% 升余弦窗函数
function window = rcoswindow(alpha, bit_length)
    warning off;
    window = zeros(1, bit_length / 2);
    t = 1:bit_length / 2;
    T = bit_length / (2 * (1 + alpha));
    window(t) = 0.5 * (1 - sin(pi / (2 * alpha * T) * (t - T)));
    window(1:(1 - alpha) * T) = 1;
    window = [fliplr(window) window]';
end
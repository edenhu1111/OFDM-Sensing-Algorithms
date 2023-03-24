%% Sensing Channel Generation
% Description: Generating a simulated time shifted and doppler modulaed
% signal.
% TxSignal_cp: transmit signal
%
% range: target range. When it is a row vector, it indicated that there are
% more than one target;
%
% velocity:target relative velocity. Its size should be the same as "range"
%
% SNR: Signal-to-noise Ratio
%
% Author: Yunbo HU (SIMIT, UCAS)
% GitHub: https://github.com/edenhu1111
%% Code
function RxSignal = sensingSignalGen(TxSignal_cp,range,velocity,SNR)
    global c0 lambda M delta_f
    delay = round(2 * range / c0 * delta_f * M);
    h_gain = exp(1j*2*pi*rand(size(delay)));    
    doppler = 2*velocity/lambda;
    max_delay = max(delay,[],'all');
    RxSignal = zeros(size(TxSignal_cp,1)+max_delay,1);
    d = zeros(size(TxSignal_cp));
    for p = 1:length(delay)
        ii = 0:length(d)-1;
        d = exp(1j*2*pi*doppler(p)*ii'/(delta_f*M));
        RxSignal = RxSignal + h_gain(p) * ...
            [zeros(delay(p),1);...
            TxSignal_cp .* d...
            ;zeros(max_delay-delay(p),1)];
    end
    RxSignal = RxSignal + 10^(-SNR/10)*(randn(size(RxSignal)) + 1j*randn(size(RxSignal)))/sqrt(2);
end


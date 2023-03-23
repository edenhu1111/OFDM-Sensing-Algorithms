%% Cycle Cross-Correlation based sensing algorithm
% RxSignal: Received sensing signal (column vector, only in time domain)
%
% TxSignal_cp: Transmitted signal (only in time domain)
% 
% mildM: number of data per sub-block
%
% Qbar: number of overlapped data between adjacent sub-blocks
%
% mildQ: length of VCP(Virtual Cyclic Prefix)
% Code Author: Yunbo HU
    % Reference: K. Wu, J. A. Zhang, X. Huang, and Y. J. Guo, 
    %¡®Integrating Low-Complexity and Flexible Sensing Into Communication Systems¡¯, 
    %  IEEE Journal on Selected Areas in Communications

function [r_cc,RDM] = cccSensing(RxSignal,TxSignal_cp,mildM,Qbar,mildQ)
% regrouping original received signal
% adjcent groups overlap Qbar points
mildN = floor((length(TxSignal_cp)-Qbar-mildQ)/(mildM - Qbar));
Rx_sub = zeros(mildM,mildN);
Tx_sub = zeros(size(Rx_sub));
for ii = 0:size(Rx_sub,2)-1 
   Rx_sub(:,ii+1) =  RxSignal(ii*(mildM-Qbar)+1 : ii*(mildM-Qbar)+mildM ,:);
   Tx_sub(:,ii+1) =  TxSignal_cp(ii*(mildM-Qbar)+1 : ii*(mildM-Qbar)+mildM,:);
end


% add VCP
for ii = 0:size(Rx_sub,2)-1 
   Rx_sub(1:mildQ,ii+1) =  Rx_sub(1:mildQ,ii+1) + RxSignal(ii*(mildM-Qbar)+mildM+1 : ii*(mildM-Qbar)+mildM+mildQ);
end
% Cross Correlation
r_cc = ifft( fft(Rx_sub) .* conj(fft(Tx_sub)) );
RDM = fft(r_cc.',10*mildN);
end


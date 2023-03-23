%% MUSIC for OFDM sensing
% CIM: Input channel information matrix(pre-processed by known transmitted symbols)
% k:target number
% Author: Yunbo HU(SIMIT, UCAS)
function [P_music_range,P_music_velo] = MUSICforOFDMsensing(CIM,k)
global M N c0 delta_f lambda Ts
%% range estimation
R_range = CIM*CIM'/N;
[V,D]=eig(R_range);
[~,ind_D] = sort(diag(D),'descend');
U_n = V(:,ind_D(k+1:end));

delay = linspace(0,2*100/c0*delta_f,M);
ii = 0:M-1;
A = exp(-1j*2*pi*kron(ii',delay));
P_music_range = zeros(size(A,2),1);
for jj = 1:size(A,2)
    P_music_range(jj) = 1/(A(:,jj)'*(U_n*U_n')*A(:,jj));
end


%% Velocity Estimation
R_dop = CIM.'*conj(CIM)/M;
[V,D]=eig(R_dop);
[~,ind_D] = sort(diag(D),'descend');
U_n = V(:,ind_D(k+1:end));

doppler = linspace(0,2*100/lambda*Ts,M);
ii = 0:N-1;
A = exp(1j*2*pi*kron(ii',doppler));
P_music_velo = zeros(size(A,2),1);
for jj = 1:size(A,2)
    P_music_velo(jj) = 1/(A(:,jj)'*(U_n*U_n')*A(:,jj));
end


end


%% ESPRIT for OFDM sensing
% CIM: Input channel information matrix(pre-processed by known transmitted symbols)
% k:target number
% Author: Yunbo HU(SIMIT, UCAS)
% GitHub: https://github.com/edenhu1111
function [range,velocity] = ESPRITforOFDMsensing(CIM,k)
global lambda delta_f c0 Ts
[M,N] = size(CIM);
%% range estimation
z = [CIM(1:M-1,:);CIM(2:M  ,:)];
R_zz = z*z'/N;
[U,~,~] = svd(R_zz);
Es = U(:,1:k);
Esx = Es(1:M-1,:);
Esy = Es(M:end,:);

EE = [Esx,Esy];
[F,~,~] = svd(EE'*EE);
F = F(:,end-k+1:end);
F1 = F(1:k,:);
F2 = F(k+1:2*k,:);
psi = -F1*inv(F2);
[~,D] = eig(psi);


phi = angle(diag(D));
phi(phi>0) = phi(phi>0) - 2*pi;
tau = -phi/(2*pi*delta_f);
range = tau*c0/2;

%% doppler estimation
z = [CIM(:,1:N-1),CIM(:,2:N)];
R_zz = z.'*conj(z)/M;
[U,~,~] = svd(R_zz);
Es = U(:,1:k);
Esx = Es(1:N-1,:);
Esy = Es(N:end,:);

EE = [Esx,Esy];
[F,~,~] = svd(EE'*EE);
F = F(:,end-k+1:end);
F1 = F(1:k,:);
F2 = F(k+1:2*k,:);
psi = -F1*inv(F2);
[~,D] = eig(psi);

phi = angle(diag(D));
phi(phi<0) = phi(phi<0) + 2*pi;
doppler = phi/(2*pi*Ts);
velocity = doppler*lambda/2;

end


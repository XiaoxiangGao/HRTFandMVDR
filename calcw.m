% This file gives an example which calculates the HRTF of a rigid sphere. 
% The result is then used to compute binaural MVDR beamformer coefficients.
% Refer to 'Dual-Channel Speech Enhancement by Superdirective Beamforming'.

clear
clc

%% Set parameters

betamin     = 0.1;
thetamin    = 150 * pi/180;        % 150 degrees
FS          = 16e3;
NFFT        = 512;
f           = 0 : FS/NFFT : FS-FS/NFFT;
omega       = 2*pi*f;
c           = 343;
a           = 8.75e-2;
OMEGA0      = c/a;
NF          = length(omega);
EPSI        = 0;


%% Compute cross-spectral-density matrix
% In pratice, only one element is needed

theta       = -pi : 5*pi/180 : pi-5*pi/180;
midthe      = [pi+theta(1:37), theta(38:end)-pi];
NTH         = length(theta);
gamma_mod   = zeros(1,NTH);
tau_mod     = zeros(1,NTH);
D_L         = zeros(NTH,NF);
D_R         = zeros(NTH,NF);

for nth = 1:NTH
    % compute D_R
    gamma_mod(nth) = (1+betamin/2) + (1-betamin/2) * cos(theta(nth)*pi/thetamin);
    if (0 - EPSI)<= abs(theta(nth)) < (pi/2 + EPSI)
        tau_mod(nth) = -a * cos(theta(nth)) / c;
    else
        tau_mod(nth) = a * (abs(theta(nth)) - pi/2) / c;
    end
    
    for k = 1:NF
        D_R(nth,k) = (1 + 1i*gamma_mod(nth)*omega(k)/(2*OMEGA0)) * exp(-1i*omega(k)*tau_mod(nth)) / (1+1i*omega(k)/(2*OMEGA0));
    end
    
    % compute D_L
    midgamma = (1+betamin/2) + (1-betamin/2) * cos(midthe(nth)*pi/thetamin);
    if (0 - EPSI)<= abs(midthe(nth)) < (pi/2 + EPSI)
        midtau = -a * cos(midthe(nth)) / c;
    else
        midtau = a * (abs(midthe(nth)) - pi/2) / c;
    end
    for k = 1:NF
        D_L(nth,k) = (1 + 1i*(midgamma*omega(k)/(2*OMEGA0))) * exp(-1i*omega(k)*midtau) / (1+1i*omega(k)/(2*OMEGA0));
    end
end


phi         = zeros(1,NF);
numerator   = zeros(1,NF);
denominator = zeros(1,NF);
for k = 1:NF
    middenom_1  = 0;
    middenom_2  = 0;
    for nth = 1:NTH
        numerator(k)    = numerator(k) + D_L(nth,k) * conj(D_R(nth,k));
        middenom_1      = middenom_1 + abs(D_L(nth,k))^2;
        middenom_2      = middenom_2 + abs(D_R(nth,k))^2;
    end
    denominator(k)  = sqrt( middenom_1 * middenom_2 );
    phi(k)          = numerator(k) / denominator(k);
end


%% Compute weigth vector w in the front direction, theta_s = 0

mus     = 0;
w       = zeros(2, NFFT/2+1);
wout    = zeros((NFFT/2+1)*2,2);

for k = 1:(NFFT/2+1)
    phimat  = [1, abs(phi(k)); abs(phi(k)), 1];
    INVphi  = inv(phimat);
    Dvec    = [D_L(59,k); D_R(59,k)];
    w(:,k)  = (INVphi + mus*eye(2)) * Dvec / ( Dvec' * (INVphi + mus*eye(2)) * Dvec );
    wout(k,:) = [real(w(1,k)) imag(w(1,k))];
    wout(NFFT/2+1+k,:) = [real(w(2,k)) imag(w(2,k))];
end

w(:,1) = [0,0];
wout(1,:) = [0 0];
wout(258,:) = [0 0];
%% plot

% plot(abs(w(1,:)))

save('w_far.mat','w')
% save wfront.txt -ascii wout
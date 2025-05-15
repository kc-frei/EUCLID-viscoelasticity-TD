function [N,Tsim,uc_th,time,dt] = TD_sim_setting_v03(target_mat)

tauGuK = union(target_mat.tauG,target_mat.tauK,'sorted');
fmin = .1*1/tauGuK(end); %0.001 was the best   % actually we do not know these vcalues. Optimal value 0.1 or even 0.01*(1/tau(end)) must be used.
fmax =     1/tauGuK(1);     % actually we do not know these vcalues. Optimal value 5*(1/tau(1));  or even 10* The larger the better you get E_inf
% omemin = 2*pi*fmin;
% omemax = 2*pi*fmax;
% Tmax = 2*pi/omemin;
% Tmin = 2*pi/omemax;
Tmax = 1/fmin;
Tmin = 1/fmax;
numT  = 1; % 5 was the best!!!!
Tsim = numT*Tmax
dt = 10*Tmin%/10; % %/5 was veru good
fs = 1/dt;

fNyq = fs/2;

N = Tsim/dt;
N = N - mod(N,2);

time = (0:N-1)*dt;
Tsim = N/fs;
df = 1/Tsim;
dome = 2*pi*df;

freqvec  = [0:N-1]*df;
freqvec = freqvec(2:N/2+1);
num_exc_ome_aux = 120;
%fgen_gen = logspace(floor(log10(fmin)),ceil(log10(fmax)),num_exc_ome);
fgen_gen = df*unique(round(logspace(0,ceil(log10(ceil(fmax/df))),num_exc_ome_aux)));
num_exc_ome = length(fgen_gen);
for i = 1:num_exc_ome
    [minval,imin] = min(abs(fgen_gen(i) - freqvec))
    fgen_gen2(i) = freqvec(imin)
    indvec(i) = imin;
end
omevec = 2*pi*fgen_gen2;

% Nyquist check
if ((fgen_gen2(end) - fNyq) > 1e-6) == 1
    disp('Abort since maximum bandlimit for generarting the signl is larger than Nyquist fr.')
    return
end
% if fgen_gen2(1) > fmin || fgen_gen2(end) <= fmax
%     disp('Warning frequency band exacited inconsistent with material rel times')
% end
Amp0 = 0.05;
Amp = Amp0*ones(1,length(fgen_gen2));
S0  = Amp0^2/2/dome;
sig2 = S0*2*pi*(freqvec(end)-freqvec(1))
sig  = sqrt(sig2)

for i = 1:num_exc_ome
    Uc_th(:,i) =  Amp(i)*cos(2*pi*fgen_gen2(i)*time + pi/3);
end
%uc_th = sum(Uc_th,2);
uc_th = Amp0*ones(size(time))';

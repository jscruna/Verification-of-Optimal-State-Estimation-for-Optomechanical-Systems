%% Here we introduce the information about the optomechanical system
clear("all")
sam=100; % number of samples to consider the Monte Carlo simulation
rng shuffle
% We generate different seeds to avoid dependance between random variables
rng(1)
rand1=randn(1,sam);
rng(2)
rand2=randn(1,sam);
rng(3)
rand3=randn(1,sam);
rng(4)
rand4=randn(1,sam);
rng(5)
rand5=randn(1,sam);
rng(6)
rand6=randn(1,sam);
rng(7)
rand7=randn(1,sam);

%Optomechanical parameters

mass = 7.71e-6 + 0.01*rand1*1e-6; % mass of the mechanical oscillator with error
mass_mean = 7.71e-6; %mean value

kappa_p=2*pi*(1.64+0.02*rand2)*1e6; % optical loss of the optical cavity with error
kappa_p_mean=2*pi*1.64*1e6; %mean value

g=-2*pi*(3.2+0.2*rand3)*1e4; %optomechanical coupling with error
g_mean=-2*pi*3.2*1e4; %mean value

omega_m=2*pi*(280+7*rand4); %resonance frequency of the mechanical oscillator
omega_m_mean=2*pi*280; %mean value

Detune=0.0292+0.0004*rand5; % detunning of the optical cavity
Detune_mean=0.0292; %mean value

Temp=(11+2*rand6)*1e-3; %Temperature of the mechanical system
Temp_mean=11e-3;

n_th=(8+2*rand7)*1e5; %Mean phonon number occupation of the mechanical system
n_th_mean=8e5;

gamma_m = 1.1*2*pi;%7.04; %mechanical damping
gamma_m_mean = 1.1*2*pi;

finesse = 1849;
beta = 0.631;
N_th = 19;
hbar = 6.626e-34 / (2 * pi); % Planck's constant
k_B = 1.380648e-23; % Boltzman's constant
c = 299.792458e6; % Speed of light
lambda = 1064e-9; % Wavelenth
omega_laser = c / lambda * 2 * pi;
P_in = 30e-3 * 0.8453; % Initial intensity of the laser
L = 0.098; % Length of the cavity
eta = 0.92; % Efficiency

kappa_in_over_all = 0.2247; % normalized optical loss using detuning

x_zpf = sqrt(hbar ./ 2 ./mass ./omega_m); % Zero-point fluctuation for position with error
x_zpf_mean=sqrt(hbar/2/mass_mean/omega_m_mean); % Mean value
%p_zpf = sqrt(mean(mass) *mean(omega_m) * hbar / 2); % Zero-point fluctuation for momentum with error
p_zpf_mean=sqrt(mass_mean*omega_m_mean*hbar/2); % Mean Value

%Wiener Filter parameters

n_x = 2 * gamma_m * (2 * n_th + 1) + 16 * g.^2 * (2 * N_th + 1) ./ ((1 + 4 * Detune.^2).*kappa_p);% Total phonon number occupation in the mechanical system including back-action
n_x_mean= 2*gamma_m*(2*n_th_mean+1)+16*g_mean^2*(2*N_th +1 )/((1+4*Detune_mean^2)*kappa_p_mean);
lambda_x = 64 * g.^2 .* eta.* Detune.^2 ./((2 * eta * N_th + 1).* kappa_p.*(1 + 4 * Detune.^2).^2);% Measurement rate of the mechanical (capacity to resolve the zero-point fluctuation of the mechanical system)
lambda_x_mean=64*g_mean^2*eta*Detune_mean^2/((2*eta*N_th + 1)*kappa_p_mean*(1 + 4*Detune_mean^2)^2);
sigma_x = (-32 * g.^2 .* eta.* Detune) .* (2 * N_th + 1)./ ((1 + 4.* Detune.^2).^2.* kappa_p * (2 * eta * N_th + 1));% It's correlated with the degree of squeezing in the mechanical system
sigma_x_mean= (-32 *g_mean^2*eta*Detune_mean)*(2*N_th + 1)/((1 + 4*Detune_mean^2)^2*kappa_p_mean*(2*eta*N_th+1));
omega_x = sqrt(sqrt(omega_m.^4 + 2 * sigma_x.* omega_m.^3 + n_x .* lambda_x.* omega_m.^2)); % Modified resonance frequency due to the measurement
omega_x_mean= sqrt(sqrt(omega_m_mean^4 + 2 * sigma_x_mean* omega_m_mean^3 + n_x_mean* lambda_x_mean* omega_m_mean^2));
gamma_x = sqrt(gamma_m.^2 - 2 * omega_m.* (omega_m + sigma_x) + 2 * omega_x.^2); % Modified mechanical damping due to measurement
gamma_x_mean = sqrt(gamma_m^2 - 2*omega_m_mean*(omega_m_mean + sigma_x_mean) + 2 * omega_x_mean^2);
M = (2*N_th*eta+1);

%% Theoretical values of the variances for prediction and retrodiction

V11_aux=(gamma_x_mean-gamma_m)/lambda_x_mean; % Variance in position of the mechanical system after considering Bayes' theorem for a causal optomechanical system
VE11_aux=(gamma_x_mean+gamma_m)/lambda_x_mean; % Variance in position of the mechanical system after considering Bayes' theorem for an anti-causal optomechanical system

V22_aux=(gamma_x_mean-gamma_m)/(lambda_x_mean*2*omega_m_mean^2)*(2*omega_m_mean*(omega_m_mean + sigma_x_mean) + gamma_x_mean*(gamma_x_mean-gamma_m));% Variance in momentum of the mechanical system after considering Bayes' theorem for an optomechanical system
VE22_aux=(gamma_x_mean+gamma_m)/(lambda_x_mean*2*omega_m_mean^2)*(2*omega_m_mean*(omega_m_mean+ sigma_x_mean) + gamma_x_mean*(gamma_x_mean+gamma_m));% Variance in momentum of the mechanical system after considering Bayes' theorem for an optomechanical system

V12_aux=(gamma_x_mean-gamma_m)^2/(2*omega_m_mean*lambda_x_mean);% Covariance of the mechanical system after considering Bayes' theorem for an optomechanical system
VE12_aux=-(gamma_x_mean+gamma_m)^2/(2*omega_m_mean*lambda_x_mean); %Covariance of the mechanical system after considering Bayes' theorem for an optomechanical system


%%
%measurement data 11mK

load('your_experimental_data.mat');

%% PSD generation for a fix value (mean value) of prediction and retrodiction (Fig. 3)
%Wiener filter position for prediction (using mean value)
Hq_pred_fix=@(omega) (1/sqrt(lambda_x_mean*M))*((omega_x_mean^2-omega_m_mean^2)-i*omega*(gamma_x_mean-gamma_m))./((omega_x_mean^2-omega.^2)+i*gamma_x_mean*omega);
%Wiener filter position for retrodiction (using mean value)
Hq_retr_fix=@(omega) (1/sqrt(lambda_x_mean*M))*((omega_x_mean^2-omega_m_mean^2)+i*omega*(gamma_x_mean+gamma_m))./conj((omega_x_mean^2-omega.^2)+i*gamma_x_mean*omega);
%Wiener Filter momentum for prediction
Hp_pred_fix=@(omega) 1/(sqrt(lambda_x_mean*M)*omega_m_mean)*(-(gamma_x_mean-gamma_m)*omega_m_mean^2-i*omega*(omega_x_mean^2-omega_m_mean^2+(gamma_x_mean-gamma_m)^2-(gamma_x_mean-gamma_m)*gamma_x_mean))./((omega_x_mean^2-omega.^2)+i*gamma_x_mean*omega);
%Wiener Filter momentum for retrodiction
Hp_retr_fix=@(omega) 1/(sqrt(lambda_x_mean*M)*omega_m_mean)*((gamma_x_mean + gamma_m)*omega_m_mean^2-i*omega*(omega_x_mean^2-omega_m_mean^2+(gamma_x_mean+gamma_m)^2-(gamma_x_mean+gamma_m)*gamma_x_mean))./conj((omega_x_mean^2-omega.^2)+i*gamma_x_mean*omega);

%% Verification process under different prediction and retrodiction filters (Fig. 4)

delete(gcp('nocreate'))
%sample=7;
fs=1e6;% frequency sample of the data
resol=linspace(1,120,5);
%initialization of the variances
var_test=zeros(length(resol),sam,1);  % Variable to store the information about the position when we test information about verification process
var_test_p=zeros(length(resol),sam,1); % Variable to store the information about the momentum when we test the information about verification process
var_test_qp=zeros(length(resol),sam,1); % Variable to store the information about the cross-variance
parpool(5)
parfor j=1:sam % loop for the Monte Carlo simulation considering the error of the experimental data
rPD2=rPD;
fmin=1;
%causal notch
for n_notch=1:1:40
    [b,a] = butter(1,[(n_notch+1)*50-fmin/2 (n_notch+1)*50+fmin/2]/(fs/2),'stop');
    rPD2=filter(b,a,rPD2);
end
cal=-pi*c*mass(j)/finesse/cos(beta)*(1-kappa_in_over_all)*(omega_m(j)^2*0.73*500); % Calibration factor to convert units (from Volts to Meters)
rPD2=rPD2/cal; % Applying convertion factor
rPD2=rPD2/x_zpf(j); % Renormalizing using the zero-point fluctuation in position
obs2=rPD2*sqrt(lambda_x(j)*M); % X quadrature ( outcome from the photodetector)

%Wiener Filter for position
Hq_pred=@(omega) (1/sqrt(lambda_x(j)*M))*((omega_x(j)^2-omega_m(j)^2)-i*omega*(gamma_x(j)-gamma_m))./((omega_x(j)^2-omega.^2)+i*gamma_x(j)*omega);
Hq_retr=@(omega) (1/sqrt(lambda_x(j)*M))*((omega_x(j)^2-omega_m(j)^2)+i*omega*(gamma_x(j)+gamma_m))./conj((omega_x(j)^2-omega.^2)+i*gamma_x(j)*omega);
%Wiener Filter for momentum
Hp_pred=@(omega) 1/(sqrt(lambda_x(j)*M)*omega_m(j))*(-(gamma_x(j)-gamma_m)*omega_m(j)^2-i*omega*(omega_x(j)^2-omega_m(j)^2+(gamma_x(j)-gamma_m)^2-(gamma_x(j)-gamma_m)*gamma_x(j)))./((omega_x(j)^2-omega.^2)+i*gamma_x(j)*omega);
Hp_retr=@(omega) 1/(sqrt(lambda_x(j)*M)*omega_m(j))*((gamma_x(j) + gamma_m)*omega_m(j)^2-i*omega*(omega_x(j)^2-omega_m(j)^2+(gamma_x(j)+gamma_m)^2-(gamma_x(j)+gamma_m)*gamma_x(j)))./conj((omega_x(j)^2-omega.^2)+i*gamma_x(j)*omega);

%resol=linspace(1,700,70);
%resol=linspace(1,100,10);
fresol=zeros(1,length(resol)); % Range for the frequency resolution

for mm=1:5% PUT THE SIZE OF THE FREQ. RESOLUTION SIZE

a='red';
fmin=resol(mm);
fresol(mm)=fmin;
[ Pxx5, Freq5, n ] = Pxx_fun( obs2, fs, fmin, a );
%position verification
var_test(mm,j,1)=bandpower(Hq_retr(2*pi*Freq5).*conj(Hq_retr(2*pi*Freq5)).*Pxx5-Hq_pred(2*pi*Freq5).*conj(Hq_pred(2*pi*Freq5)).*Pxx5,Freq5,[200 fs/1000],'psd');
%momentum verification
var_test_p(mm,j,1)=bandpower(Hp_retr(2*pi*Freq5).*conj(Hp_retr(2*pi*Freq5)).*Pxx5-Hp_pred(2*pi*Freq5).*conj(Hp_pred(2*pi*Freq5)).*Pxx5,Freq5,[200 fs/1000],'psd');
%covariance verification
var_test_qp(mm,j,1)=bandpower(real(Hq_retr(2*pi*Freq5).*conj(Hp_retr(2*pi*Freq5)).*Pxx5-Hq_pred(2*pi*Freq5).*conj(Hp_pred(2*pi*Freq5)).*Pxx5),Freq5,[200 fs/1000],'psd');
end

end

%% Setting the intervals
%fresol=zeros(1,length(resol));
x=resol;
y_mean=zeros(1,length(resol));
y_p_mean=zeros(1,length(resol));
y_qp_mean=zeros(1,length(resol));
y_upper=zeros(1,length(resol));
y_lower=zeros(1,length(resol));
y_upper_p=zeros(1,length(resol));
y_lower_p=zeros(1,length(resol));
y_upper_qp=zeros(1,length(resol));
y_lower_qp=zeros(1,length(resol));
for kk=1:length(resol)
y_upper(kk)=mean(var_test(kk,:))+std(var_test(kk,:));
y_lower(kk)=mean(var_test(kk,:))-std(var_test(kk,:));
y_mean(kk)=mean(var_test(kk,:));
y_upper_p(kk)=mean(var_test_p(kk,:))+std(var_test_p(kk,:));
y_lower_p(kk)=mean(var_test_p(kk,:))-std(var_test_p(kk,:));
y_p_mean(kk)=mean(var_test_p(kk,:));
y_upper_qp(kk)=mean(var_test_qp(kk,:))+std(var_test_qp(kk,:));
y_lower_qp(kk)=mean(var_test_qp(kk,:))-std(var_test_qp(kk,:));
y_qp_mean(kk)=mean(var_test_qp(kk,:));
end

%% Plot using verification process and raw data vs. the theoretical value
figure()

fill([x, fliplr(x)], [y_upper, fliplr(y_lower)], [0.5 0.1 1], 'EdgeColor', 'none');hold on;
plot(x,(2*V11_aux)*ones(1,length(resol)),'LineWidth',2,'LineStyle','--',...
    'Color',[1 0.411764705882353 0.16078431372549])
%legend({'Experiment','Theory'},'Orientation','horizontal')

fill([x, fliplr(x)], [y_upper_p, fliplr(y_lower_p)], [0.5 0.1 1], 'EdgeColor', 'none');hold on;

plot(x,(2*V22_aux)*ones(1,length(resol)),'LineWidth',2,'LineStyle','--',...
    'Color',[1 0.411764705882353 0.16078431372549]);grid on;
legend({'Experiment','theory'},'Orientation','horizontal')
xticks([1 10 20 30 40 50 60 70 80 90 100 110 120])
xticklabels({'1','10','20','30','40','50','60','70','80','90','100','110','120'})
xlim([1 120])
ylim([0 2.2e4])
ylabel('Verification Variance','FontSize',12,'FontName','Times New Roman')
xlabel('Frequency resolution (Hz)','FontSize',12,'FontName','Times New Roman')

%%Here we define the function for the power spectral density
function [ Pyy_ave, f, n ] = Pxx_fun( posx, fs, fmin, a )
overlap=0.5;
N=ceil(fs/fmin);
n=pow2(nextpow2(N));
%num_windows=floor((fs*T/n-1)/(1-overlap))+1;
[Pyy_ave,f]= pwelch(posx,hanning(N),[],n ,fs,'onesided' );
%figure()
%loglog(f,Pyy_ave,'color',a)
end

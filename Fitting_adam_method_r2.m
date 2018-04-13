% 17.06.16 - agh load existing transfer function data, fit using simple VE
% model
% r2: fit Loic's multiharmonic data

clear all, close all, warning off
addpath('~/Documents/postdoc/analysis')
addpath('../Step4_LielegOptimization_161207')


%datdir = '/Volumes/Hossein//PostdocJournal/Optical trap/Selected for Paper/';
%datdir = '~/Documents/McGill/data/Hossein/Selected for Paper/';
%flist{1} = {'WT/1/v1_comprehensive_AGH37.mat','WT/2/v1_comprehensive_AGH30.mat','WT/3/3_C1_Jan15_P4.mat','WT/4/C1_Feb3_P3.mat','WT/C1_Dec1_P2.mat'};
%flist{2} = {'K255E/1/v1_comprehensive_AGH19.mat','K255E/2/v1_comprehensive_AGH30.mat','K255E/3/v1_comprehensive_AGH26.mat','K255E/4/4_C193_dec15_p4.mat'...
%    'K255E/7/C193_Jan27_p3.mat','K255E/9/C193_March2_P3.mat'};
%Wfrf=50e-3;

datdir = '~/Documents/McGill/data/loic-chaubet/selected_FINAL/';
flist{1} = {'Nov_24/C1/C1_P6_p2_Nov_24_2017_0.4amp_X1.0X_17freq0_amp_corrected.mat','Nov_27/C1/C1_P7_p1_Nov_27_2017_0.3amp_X1.0X_17freq_amp_corrected.mat',...
    'Nov_27/C1/C1_P7_p1_Nov_27_2017_0.3amp_X1.2X_17freq_amp_corrected.mat'};
flist{2} = {'Nov_24/C193/C193_P5_p2_Nov_24_2017_0.4amp_X1.0X_17freq0_amp_corrected.mat','Nov_24/C193/C193_P5_p2_Nov_24_2017_0.4amp_X1.2X_17freq0_amp_corrected.mat',...
    'Nov_27/C193/p1/C193_P6_p1_Nov_27_2017_0.3amp_X1.0X_17freq_amp_corrected.mat','Nov_27/C193/p2/C193_P6_p2_Nov_27_2017_0.4amp_X1.0X_17freq_amp_corrected.mat',...
    'Nov_27/C193/p2/C193_P6_p2_Nov_27_2017_0.4amp_X1.2X_17freq_amp_corrected.mat','Nov_29_and_Dec_12/C193/BB_C193_P6_p2_Nov_29_2017_0.4_1.0x_17freq_FINAL_best.mat'};
Wfrf=0.5e-3;
fmin=400;
fmax=8000;

zve=1;
zph=0;

filename = fullfile(datdir,flist{2}{6});
load(filename,...
      'Fexc','H','FP','PYm','Pyy','Pyy2','fsamp','Rbead','Zbead','C','kT');

i = sqrt(-1);

%Wfrf=50e-3;%[10e-3, 100e-3] *(numel(Fexc)./numel(fdata1));

jph = find(Fexc<100);
kph = round(mean(angle(H(jph)))/pi); %wrap phase to zero
%H = H.*10.*exp(-i*pi).*exp(i*0.000048*2*pi.*Fexc); %correct for polarity of AOD and QPD; phase lag between acquisition channels
H = H.*exp(i*kph*pi).*exp(i*0.00005*2*pi.*Fexc); %c
%H = H;

jkeep = find(abs(H) - mean(abs(H)) < std(abs(H)));
H = H(jkeep);
Fexc = Fexc(jkeep);
C = C(jkeep);

Mfrf = 20*log10(abs(H));
PHfrf = (180/pi).*unwrap(angle(H));
PYm = 2*PYm;
Hall = H; %record for G* calculation
Fexcall=Fexc;
Fexc_PS_reject=[Fexc,56.5,85,92,132,162.5,169.5,337,899,1079,1347,2441]; %spikes in PS to remove 
Call = C;


% plot spectra ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%filename='C1_Dec1_P4';
%close all;
%falias=5e3;%6e3;%20e3;

figure, subplot(311), semilogx(Fexc,Mfrf,'o','linewidth',2)
ylabel('Mag, dB'), set(gca,'Xlim',[0.01 5e3],'fontsize',14,'fontweight','bold');
subplot(312), semilogx(Fexc,PHfrf,'o','linewidth',2)
ylabel('Phase, deg'), xlabel('f, Hz'), set(gca,'Xlim',[0.01 5e3],'fontsize',14,'fontweight','bold');
hfrf=gcf;

subplot(313), semilogx(Fexc,C,'o','linewidth',2)
ylabel('Coherence'), xlabel('f, Hz'), set(gca,'fontsize',14,'fontweight','bold')
%disp('Click on minimum acctable coherence')
%[chert,cohthreshh]=ginput(1);
cohthreshh = 0.9;
disp('Click on minimum acceptable frequency')
[freq_low,chert]=ginput(1);
disp('Click on maximum acceptable frequency')
[freq_high, chert]=ginput(1);


find_coh=find(C>cohthreshh&Fexc>freq_low&Fexc<freq_high);
Fexc=Fexc(find_coh);
H=H(find_coh);
PHfrf=PHfrf(find_coh)%*180/pi;  % to convert from radian to degree
C=C(find_coh);
figure(hfrf), subplot(311);
Mfrf=Mfrf(find_coh);
Pyy=Pyy(find_coh);
Pyy2=Pyy2(find_coh);

% finish data preparation%%%%%
%close all;
figure(hfrf), hold on, subplot(311), semilogx(Fexc,Mfrf,'o','markerfacecolor','b','linewidth',2)
ylabel('Mag, dB'), set(gca,'Xlim',[0.01 5e3]);
subplot(312), hold on, semilogx(Fexc,PHfrf,'o','markerfacecolor','b','linewidth',2)
ylabel('Phase, deg'), xlabel('f, Hz'), set(gca,'Xlim',[0.01 5e3]);
subplot(313), hold on, semilogx(Fexc,C,'o','markerfacecolor','b','linewidth',2)
ylabel('Coherence'), xlabel('f, Hz')


if zve==0|zve==2,
    ig.fz_frf=0;
    ig.Pz=0;
elseif zve==1,
    disp('Click on zero corner frequeny (Hz)')
    [ig.fz_frf,ig.Pz]=ginput(1);
end
disp('Click on pole corner frequency (Hz)')
[ig.fp_frf,ig.MinfdB]=ginput(1);

figure, loglog(FP,PYm,'b','linewidth',2), hps=gcf;
xlabel('f (Hz)'), ylabel('Power Spectrum (V^2*s)')

% disp('Click on corner frequency (Hz) and steady-state amplitude ($V^2$)')
% [ig.fp_ps,ig.Pss]=ginput(1);
%ig.Dfit=ig.Pss*pi^2*(400+ig.fp_frf^2);


% fit frequency response and power spectrum  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ig.Minf=10^(ig.MinfdB/20);
ig.Gamma0=9.42e-6; %[pN*s/nm]
%disp('Must fix: Calculate Gamma0(T)')

% find slope of power spectrum at high F %%%

disp('Click on min and max freq. to calculate slope. The max frequency will also be use as maximum limit of fitting')
[slope_psmin,ymin]=ginput(1);
[slope_psmax,ymax]=ginput(1);  %%% hoddrin changed this

disp('Click on corner frequency of QPD response' )
[falias,yfalias]=ginput(1);%6e3;%20e3;



jps=find(FP>=slope_psmin & FP<=slope_psmax);
pslope=polyfit(log(2*pi.*FP(jps)),log(PYm(jps)),1);
ig.Phi=-pslope(1);
hold on, loglog(FP,exp(polyval(pslope,log(2*pi.*FP))),'r')


ig.Pss=exp(polyval(pslope,log(2*pi*ig.fp_frf)));
ig.Dfit=ig.Pss*pi^2*((3*ig.fp_frf)^2);

Hth=tf(ig.Minf*[1 2*pi*ig.fz_frf],[1 2*pi*ig.fp_frf]);
[Hm,Hph,Hw]=bode(Hth);

% figure, subplot(311), semilogx(Fexc,Mfrf,'o','linewidth',2)
% ylabel('Mag, dB'), set(gca,'Xlim',[0.01 5e3]);
% subplot(312), semilogx(Fexc,PHfrf,'o','linewidth',2)
% ylabel('Phase, deg'), xlabel('f, Hz'), set(gca,'Xlim',[0.01 5e3]);
% subplot(313), semilogx(Fexc,C,'o','linewidth',2)
% ylabel('Coherence'), xlabel('f, Hz')
% hfrf=gcf;

% figure(hfrf), subplot(311), hold on,
% semilogx(Hw./(2*pi),20.*log10(Hm(:)),'g')
% subplot(312), hold on,
% semilogx(Hw./(2*pi),Hph(:),'g')

ig.Beta=1000;%1/ig.Minf;
ig.GammaR=5;%(kT/(ig.Beta^2*ig.Dfit))/ig.Gamma0;


ig.Alpha=min([abs(pslope(1))-1,1]);

wc=2*pi*ig.fp_frf;

ig.k_cyt1=1e-4;%0.01;%max([0.001,sqrt(4*kT*ig.GammaR*ig.Gamma0/ig.Pss)-ig.k_trap-ig.k_cyt0]);
ig.k_cyt0=ig.fz_frf*2*pi*ig.GammaR*ig.Gamma0;
ig.k_trap=ig.fp_frf*2*pi*ig.GammaR*ig.Gamma0-ig.k_cyt0%0.05;
ig.keq0=ig.k_trap+ig.k_cyt0+ig.k_cyt1*(ig.fp_frf)^ig.Alpha;

    %   Beta    GammaR    Alpha    k_trap   k_cyt0   k_cyt1 mass    nu];
lb=[0.5.*ig.Beta   0.5*ig.GammaR     ig.Alpha-0.1        0.1.*ig.k_trap         0.1*ig.k_cyt0    0      1e-1    0];
ub=[1.5*ig.Beta     100*ig.GammaR     ig.Alpha+0.1        10.*ig.k_trap         10.*ig.k_cyt0    0.1    1e1     1.01];

%Beta    GammaR    Alpha    k_trap   k_cyt0   k_cyt1 mass    nu];
%lb=[0.1.*ig.Beta   0.01*ig.GammaR    0.5       0.005     0.001        0.001       1e-1    1e-6];
%ub=[10*ig.Beta     10*ig.GammaR     0.9       0.2   10.*ig.k_cyt0    1e-1    1e1    1.01];
%ub=[10*ig.Beta     10*ig.GammaR     ig.Alpha+0.05       0.5   10.*ig.k_cyt0    1e-1    1e1    1.01];


%disp([' - fit using ',num2str(numel(u10)),' parameters -'])
%options=optimset('TolFun',1e-8,'TolX',1e-6,'Display','on','MaxFunEvals',1e14,'MaxIter',1e4);%,...
%'LargeScale','off','LevenbergMarquardt','on','DiffMinChange',0.1);%,'DiffMinChange',5);%,'MaxFunEvals',1e3,'Display','iter');
%fmin=20;
%fmax=5000;
%fmax=3000;
%jfit1=find((FP>fmin & FP<fmax) & ~(FP>fnl1&FP<fnh1)& ~(FP>fnl2&FP<fnh2)& ~(FP>fnl3&FP<fnh3)& PYm>0 & isnan(PYm)~=1);
%jfit1=find((FP>fmin & FP<fmax) & PYm>0 & isnan(PYm)~=1)
jfit1 = find(FP>fmin & FP<fmax & isnan(PYm)~=1);%& FP~=Fexcall);
for ke = 1:numel(Fexc_PS_reject),
   je = find(abs(FP(jfit1)-Fexc_PS_reject(ke))>0.01*Fexc_PS_reject(ke));
   jfit1=jfit1(je);
end

fdata1=FP(jfit1);%(Fexc+10):10:7e3;
ydata1=cumtrapz(FP(jfit1),PYm(jfit1));%interp1(fx,IPSx,fdata,'pchip')
ydat10=ydata1;

figure, semilogx(fdata1,ydata1,'b','linewidth',2), hips=gcf; %,fdata1,IPth1,'g','linewidth',1
xlabel('f (Hz)'), ylabel('$\int P_{th}$ ($V^2$)')
axis tight

%Beta    GammaR    Alpha    k_trap   k_cyt0   k_cyt1 mass    nu
u0=[ig.Beta, ig.GammaR, ig.Alpha, ig.k_trap, ig.k_cyt0, ig.k_cyt1, 10, 1];

options=optimset('TolFun',1e-9,'TolX',1e-9,'Display','iter','MaxFunEvals',1e5,'MaxIter',2e3);%,'DiffMinChange',1e-3);
%options=optimset(options,'Largescale','on');
%options2=optimset(options,'Algorithm',{'levenberg-marquardt',0.01})
%options=optimset(options,'LargeScale','off','DiffMinChange',1e0,'DiffMaxChange',1e1);
%options2=optimset('Algorithm',{'levenberg-marquardt',0.01},'TolFun',1e-9,'TolX',1e-9,'Display','iter','MaxFunEvals',1e5,'MaxIter',2e3);%,'DiffMinChange',5e-3,'DiffMaxChange',1e2,'MaxIter',5e3,'Largescale','off');
options2=optimset('TolFun',1e-9,'TolX',1e-9,'Display','iter','MaxFunEvals',1e5,'MaxIter',2e3);%,'DiffMinChange',1e-6);
%options2=optimset(options,'LargeScale','off','DiffMinChange',1e-3,'DiffMaxChange',1,'TolFun',1e-10,'TolX',1e-9);

%



ff=1:10:2e4;
if zve==0,
%     jve=[1,2,4,7,8];
%     u0=u0(jve);
%     lb=lb(jve);
%     ub=ub(jve);
%     %u0(4)=0; lb(4)=0; ub(4)=1e-10;
%     xdata=[Fexc,Fexc,fdata1'];
%     ydata=[zeros(size(H)),zeros(size(H)),zeros(size(ydat10'))];
%     W=Wfrf.*C(:)';%.^25;
%     lb=[0.1.*ig.Beta   0.5*ig.GammaR        0   1e-6    0];
%     ub=[10*ig.Beta     100*ig.GammaR        1   1e3     1.01];
%     [u2i,resnorm,residual]=lsqcurvefit(@theor_trap_frf_nove_r6,u0,xdata,ydata,lb,ub,options,abs(H(:))',PHfrf,W,ydat10',fsamp,Rbead,Zbead,kT,falias);
%     lb2=[0.8*u2i(1)     0.8*u2i(2)        0     0.1    0];
%     ub2=[1.2*u2i(1)     1.2*u2i(2)        1     1e3     1.01];
%     [u2,resnorm,residual]=lsqcurvefit(@theor_trap_frf_nove_r6,u2i,xdata,ydata,lb2,ub2,options2,abs(H(:))',PHfrf,W,ydat10',fsamp,Rbead,Zbead,kT,falias);
% 
%     errfit=theor_trap_frf_nove_r6(u2,xdata,abs(H(:))',PHfrf,W,ydat10',fsamp,Rbead,Zbead,kT,falias);
%     ffrf_fit=Fexc;
%     jfrf=1:numel(Fexc);
%     jph=jfrf+numel(Fexc);
%     jps=(jph(end)+1):numel(errfit);
%     u2f=zeros(1,8);
%     u2f(jve)=u2;
%     %u2f(6)=u2(3);
%     u2=u2f;
%     [mf,pf]=frf_r6(u2,ffrf_fit,Rbead,Zbead,kT,falias);
%     Mfrf_fit=mf;%((errfit(jfrf)./W)+1).*abs(H(:))';
%     PHfrf_fit=pf;%((errfit(jph).*10./W)+1).*PHfrf;
%     IPS_fit=(errfit(jps)+1).*ydat10';
%     PS_fit=ps_r6(u2,ff,fsamp,Rbead,Zbead,kT,falias);
elseif zve==1,
    if zph==0,
        PHfrf=[];
        xdata=[Fexc,fdata1'];
        ydata=[zeros(size(H)),zeros(size(ydat10'))];
    else
%         xdata=[Fexc,Fexc,fdata1'];
%         ydata=[zeros(size(H)),zeros(size(H)),zeros(size(ydat10'))];
    end
    W=Wfrf.*C(:)';%.^25;
    [u2i,resnorm,residual]=lsqcurvefit(@theor_trap_frf_r6,u0,xdata,ydata,lb,ub,options,abs(H(:))',PHfrf,W,ydat10',fsamp,Rbead,Zbead,kT,falias);

    [mf0,phf0]=frf_r6(u0,Fexc,Rbead,Zbead,kT,falias);
    psf0=ps_r6(u0,ff,fsamp,Rbead,Zbead,kT,falias);
    figure(hfrf), subplot(311), hold on,
    semilogx(Fexc,20.*log10(mf0),'k','linewidth',2),set(gca,'fontsize',14,'fontweight','bold')
    subplot(312), hold on,
    semilogx(Fexc,phf0,'k','linewidth',2)%+(180/pi).*unwrap(angle(exp(-j*1e-6.*Hw)))
    figure(hps), hold on
    loglog(ff,psf0,'k','linewidth',2),set(gca,'fontsize',14,'fontweight','bold')
    
    [mf1,phf1]=frf_r6(u2i,Fexc,Rbead,Zbead,kT,falias);
    psf1=ps_r6(u2i,ff,fsamp,Rbead,Zbead,kT,falias);
    figure(hfrf), subplot(311), hold on,
    semilogx(Fexc,20.*log10(mf1),'g','linewidth',2)
    set(gca,'fontsize',14,'fontweight','bold')
    subplot(312), hold on,
    semilogx(Fexc,phf1,'g','linewidth',2)%+(180/pi).*unwrap(angle(exp(-j*1e-6.*Hw)))
    figure(hps), hold on
    loglog(ff,psf1,'g','linewidth',2),set(gca,'fontsize',14,'fontweight','bold')
    %Beta    GammaR    Alpha    k_trap   k_cyt0   k_cyt1 mass    nu
    lb2=[1     0.5   0.5      0.01    1e-10    1e-14    0.01    1e-4];    
    ub2=[2e3    1e5   1        1      10       5e-1     1e2    1.01];
    W2=Wfrf.*C(:)';
    
    [u2,resnorm,residual,exitflag,output,lambda,jacobian]=lsqcurvefit(@theor_trap_frf_r6,u2i,xdata,ydata,lb2,ub2,options2,abs(H(:))',PHfrf,W2,ydat10',fsamp,Rbead,Zbead,kT,falias);
    
    u2_ci=nlparci(u2,residual,'jacobian',full(jacobian)); %nlparci to calculate confidence intervals on the fit
    
    errfit=theor_trap_frf_r6(u2,xdata,abs(H(:))',PHfrf,W,ydat10',fsamp,Rbead,Zbead,kT,falias);
    ffrf_fit=Fexc;
    jfrf=1:numel(Fexc);
    if isempty(PHfrf)==1,
        jph=[];
        %PHfrf_fit=zeros(size(Fexc));
        jps=(jfrf(end)+1):numel(errfit);
    else
        jph=jfrf+numel(Fexc);
        %PHfrf_fit=((errfit(jph).*10./W)+1).*PHfrf;
        jps=(jph(end)+1):numel(errfit);
    end

    [Mfrf_fit,PHfrf_fit]=frf_r6(u2,Fexc,Rbead,Zbead,kT,falias);
    %Mfrf_fit=((errfit(jfrf)./W)+1).*abs(H(:))';
    IPS_fit=(errfit(jps)+1).*ydat10';
    PS_fit=ps_r6(u2,ff,fsamp,Rbead,Zbead,kT,falias);

% elseif zve==2, %k_cyt0=0;
%     if zph==0,
%         PHfrf=[];
%         xdata=[Fexc,fdata1'];
%         ydata=[zeros(size(H)),zeros(size(ydat10'))];
%     else
%         xdata=[Fexc,Fexc,fdata1'];
%         ydata=[zeros(size(H)),zeros(size(H)),zeros(size(ydat10'))];
%     end
%     W=Wfrf.*C(:)';%.^25;
% 
%     jve=[1:4,6:8];
%     u0=u0(jve);
%     ub=ub(jve);
%     lb=lb(jve);
% 
%     [u2i,resnorm,residual]=lsqcurvefit(@theor_trap_frf_nok0_r6,u0,xdata,ydata,lb,ub,options,abs(H(:))',PHfrf,W.^25,ydat10',fsamp,Rbead,Zbead,kT,falias);
% 
%     %     lb2=[0.7*u2i(1)     0.7*u2i(2)      0.8*u2i(3)        0.8*u2i(4)         0.8*u2i(5)      0      0.1    0];
%     %     ub2=[1.3*u2i(1)     1.3*u2i(2)      1.2*u2i(3)        1.2*u2i(4)         1.2*u2i(5)      1      1e3     1.01];
%     lb2=[];%[0.7*u2i(1)     0.7*u2i(2)      0        0.8*u2i(4)        0.8*u2i(5)      0      0.1    0];
%     ub2=[];%[1.3*u2i(1)     1.3*u2i(2)      0.5        1.2*u2i(4)        1.2*u2i(5)      1      1e3     1.01];
%     [u2,resnorm,residual]=lsqcurvefit(@theor_trap_frf_nok0_r6,u2i,xdata,ydata,lb2,ub2,options2,abs(H(:))',PHfrf,W.^25,ydat10',fsamp,Rbead,Zbead,kT,falias);
% 
%     ffrf_fit=Fexc;
%     jfrf=1:numel(Fexc);
%     errfit=theor_trap_frf_nok0_r6(u2,xdata,abs(H(:))',PHfrf,W,ydat10',fsamp,Rbead,Zbead,kT,falias);
%     ffrf_fit=Fexc;
%     jfrf=1:numel(Fexc);
%     if isempty(PHfrf)==1,
%         jph=[];
%         %PHfrf_fit=zeros(size(Fexc));
%         jps=(jfrf(end)+1):numel(errfit);
%     else
%         jph=jfrf+numel(Fexc);
%         %PHfrf_fit=((errfit(jph).*10./W)+1).*PHfrf;
%         jps=(jph(end)+1):numel(errfit);
%     end
%     u2f=zeros(1,8);
%     u2f(jve)=u2;
%     u2=u2f;
%     [Mfrf_fit,PHfrf_fit]=frf_r6(u2,Fexc,Rbead,Zbead,kT,falias);
%     %Mfrf_fit=((errfit(jfrf)./W)+1).*abs(H(:))';
%     IPS_fit=(errfit(jps)+1).*ydat10';
%     PS_fit=ps_r6(u2,ff,fsamp,Rbead,Zbead,kT,falias);
end

Hth=tf(u2(1)*[1 2*pi*u2(3)],[1 2*pi*u2(2)]);
Hpd=tf(2*pi*u2(end),[1 2*pi*u2(end)]);
Hth=series(Hth,Hpd);
Hw=2*pi.*Fexc;
[Hm,Hph]=bode(Hth,Hw);
figure(hfrf), subplot(311), hold on,
set(gca,'fontsize',14,'fontweight','bold')
semilogx(Fexc,20.*log10(Mfrf_fit),'c','linewidth',2)
subplot(312), hold on,
semilogx(Fexc,PHfrf_fit*pi/180,'c','linewidth',2)%+(180/pi).*unwrap(angle(exp(-j*1e-6.*Hw)))
set(gca,'fontsize',14,'fontweight','bold')

figure(hps), hold on
loglog(ff,PS_fit,'c','linewidth',2)

% IPth1=cumtrapz(fdata1,Pth_t0);%excited_power_spectrum(u,fx(jfit),falias,Fexc,Adrive1,ydat0,f0);
% figure(hips), hold on, semilogx(fdata1,IPth1,'c','linewidth',2)
figure(hips), hold on, semilogx(fdata1,IPS_fit,'c','linewidth',2)

P12=20.*log10(Pyy)-20.*log10(Pyy2);
figure, semilogx(Fexc,P12,'bo')
xlabel('Frequency (Hz)'), ylabel('Pyy(Fexc) - Pyy(2*Fexc)')
set(gca,'fontsize',14,'fontweight','bold')
%calculate trap parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%k_trap=2*pi*(fc_frf-fc_zero)*Gamma;
mi=u2i(end-1); %*10^-21 g
nui=u2i(end); %*10^12 nm^2/s
betai=u2i(1);
gammari=u2i(2);
alphai=u2i(3);
k_trapi=u2i(4);
k_cyt0i=u2i(5);
%k_eq0=u2(6);
%w_0=k_eq0/(gammar*9.42e-6);
k_cyt1i=u2i(6);%(k_eq0-k_trap-k_cyt0)/(w_0^alpha);

%disp(['Dfit = ', num2str(Dfit), ' V'])
%disp(['f_c,zero = ', num2str(fc_zero), ' Hz'])
%disp(['f_c,FRF = ', num2str(fc_frf), ' Hz'])
disp(['INTERMEDIATE FIT --------------------'])
disp(['Beta = ', num2str(betai), ' nm/V'])
disp(['Gamma = ',num2str(gammari),' *gamma0'])
disp(['Alpha = ',num2str(alphai)])
disp(['k_trap = ',num2str(k_trapi),' pN/nm'])
disp(['k_cyt,0 = ', num2str(k_cyt0i), ' pN/nm'])
disp(['k_cyt,1 = ',num2str(k_cyt1i),' pN/nm'])
disp(['Beta*k_trap = ',num2str(betai*k_trapi),' pN/V'])

%k_trap=2*pi*(fc_frf-fc_zero)*Gamma;
m=u2(end-1); %*10^-21 g
nu=u2(end); %*10^12 nm^2/s
beta=u2(1);
gammar=u2(2);
alpha=u2(3);
k_trap=u2(4);
k_cyt0=u2(5);
%k_eq0=u2(6);
%w_0=k_eq0/(gammar*9.42e-6);
k_cyt1=u2(6);%(k_eq0-k_trap-k_cyt0)/(w_0^alpha);

disp(['FINAL FIT -------------------------------'])
%disp(['Dfit = ', num2str(Dfit), ' V'])
%disp(['f_c,zero = ', num2str(fc_zero), ' Hz'])
%disp(['f_c,FRF = ', num2str(fc_frf), ' Hz'])
disp(['Beta = ', num2str(beta), ' nm/V'])
disp(['Gamma = ',num2str(gammar),' *gamma0'])
disp(['Alpha = ',num2str(alpha)])
disp(['k_trap = ',num2str(k_trap),' pN/nm'])
disp(['k_cyt,0 = ', num2str(k_cyt0), ' pN/nm'])
disp(['k_cyt,1 = ',num2str(k_cyt1),' pN/nm'])
disp(['Beta*k_trap = ',num2str(beta*k_trap),' pN/V'])

k_cyt_tot=k_cyt0+k_cyt1.*(2*pi.*FP).^alpha;
Gp=k_cyt_tot./(6*pi*Rbead);
figure, loglog(FP,Gp,'linewidth',2)
xlabel('f (Hz)'), ylabel('Gp ($pN/nm^2$)')

ww=2*pi.*FP;
g=k_cyt0 + k_cyt1.*((j.*ww).^alpha)./gamma(alpha)...
    + gammar*ig.Gamma0.*j.*ww;
Gp=g./(6*pi.*Rbead);
v=imag(g)./ww;

figure, hold on
loglog(FP,real(Gp)*1e6,'color','b','linewidth',2)
ylabel('Gp - Storage Modulus')
set(gca,'Xscale','log'), set(gca,'Yscale','log')
set(gca,'fontsize',14,'fontweight','bold')
 %hold on
loglog(FP,imag(Gp)*1e6,'color','b','linewidth',2)
ylabel('Gpp - Loss Modulus')
set(gca,'Xscale','log'), set(gca,'Yscale','log')
set(gca,'fontsize',14,'fontweight','bold')
figure, hold on
loglog(FP,v,'color','b','linewidth',2)
ylabel('viscosity')
set(gca,'Xscale','log'), set(gca,'Yscale','log')
set(gca,'fontsize',14,'fontweight','bold')

% fname_save_results=[fname_save(1:end-10),'r6c_cal.mat'];
% zsr=input('Save calibration results (0/1)? ');
% if zsr==1
%     save(fname_save_results)
% end


% cd(currentFolder)
% 
% fname_save_results=[fname_save(1:end-10),'r6c_cal.mat'];
% zsr=input('Save calibration results also to the files directory(0/1)? ');
% if zsr==1
%     save(fname_save_results)
% end
% 
% fname_save_results=[fname_save(1:end-10),'r6c_cal-frfdat.mat'];
% zsr=input('Save calibration results also to the files directory(0/1)? ');
% if zsr==1
%     save(fname_save_results)
% end
%cd(currentFolder)




zsave = input('Save fit? ');
if zsave == 1,
savename=[filename,'Curvefitting_adam_method',num2str(chxpd),'.mat'];
save(savename)
end
   


% final edition Hossein K. Heris  Nov 28,2016 The AOD and stage coefficents
% edited Dec 2016. The best working version
% are the most accurate

 % Number of FFT points is proportial to the sampling rate. 
clear
close all
clc

datdir=[pwd '/'];
cd(datdir);
fname_save='Result_';

frequency_amplitude=textread([datdir,'frequency_amp_time_inbetween_stageandAOD.txt'],'');
freq_list=[frequency_amplitude(:,1)]';
AODamp=[frequency_amplitude(:,2)]';
l=size(frequency_amplitude,2);
n=size(frequency_amplitude,1);
l
if l==3
    Samp_rate(1,1:n)=50;
    AODorstage(1,1:n)=1;
elseif l==4
    Samp_rate(1,1:n)=50;
    AODorstage=[frequency_amplitude(:,4)]';
elseif l==5
Samp_rate=[frequency_amplitude(:,5)]';
AODorstage=[frequency_amplitude(:,4)]';
end
%% 1 for AOD, 0 for stage
%urheo_fname={'p1hz-p05mhz' 'p2hz-p05mhz' 'p3hz-p05mhz' 'p4hz-p05mhz' 'p5hz-p05mhz' 'p6hz-p05mhz' 'p7hz-p05mhz' 'p8hz-p05mhz' 'p9hz-p05mhz' '1hz-p05mhz' '2hz-p05mhz' '3hz-p05mhz' '4hz-p05mhz' '5hz-p05mhz' '6hz-p05mhz' '7hz-p05mhz' '8hz-p05mhz' '9hz-p05mhz' '10hz-p05mhz' '20hz-p05mhz' '30hz-p05mhz' '40hz-p05mhz' '50hz-p05mhz' '60hz-p05mhz' '70hz-p05mhz' '80hz-p05mhz' '90hz-p05mhz' '100hz-p05mhz' '200hz-p05mhz' '300hz-p05mhz' '400hz-p05mhz' '500hz-p05mhz'};                         

%freq_list=[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 2 3 4 5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];                                                 
%AODamp=[];%0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05  0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.02 0.02 0.02 0.02 0.02];  % these are in volts
%AODamp(1,1:32)=0.05;




%if_freqzero=find(freq_list==0);
%freq_list=freq_list(1:end-1) 
 
 
pwd
currentFolder = pwd
%cd('/Users/heris/Documents/MATLAB/Optical Trap/Matlab_code_frequecny_phase_K255Eproject')

Fexc=freq_list;
fmin=0.001; % range to start fit
fmax=5e3;
Fspike=[5000 6000];

cht=0 % there is no time channel
%pixelsize=105; %????????
zve=1; % 0 no viscoelastic; 1 viscoelastic
zph=0;
FRF_f=[0,5e3]; %  range to fit frequency response
%kd=1
%kf=1;
kpd=1;  % correction factor for pd, leave it one       ????? ask adam what this is?
fmin=3e-3; % it was 3e2 %range to start fit (PS)
fmax=1e4;
Wfrf=0.1   %0.2;%0.05;%0.2  needs to be adjusted. weighting funtion when fitting frequency and power spectrum data . Different number of points in data sets. 
falias=3e5;  % QPD has filtering effect. first order low pass filter
Temp=37;
checkAOD=5; % this is the column where AOD oscillation is written in the text file
chxpd=2; % this is the associed column for QPD data
checkstage=4; % this is channel for excitation of stage
            

Qo=1e-5;% source of vibration at low frequency 
f0=0;%1;%20;%5; %[Hz] building vibrations (0 if none)

Rbead=250; %Radius of bead (nm)
Zbead=2500; %height of bead (nm)

viscosity0=1e-9; %[pN*s/nm^2]
kT=1.381e-2*(Temp+273.15); %[pN*nm]  %????
kp=0;         %%%%?????     
            
            
 % Excitation parameters
 %cx=cx_list(kf);%20.37*(fm_list(kf)/16); %[nm/V] aod sensitivity, 20.37*(Freq. multiplier/16)
           
cx=-2*1000; %%*pixelsize in nm;    %%% ??? this must change
cxstage=20000; % this is voltage to nm conversion for stage (10 V 200)
Kf=length(freq_list);          
% if ~isempty(if_freqzero)
%     Kf=Kf+1

% end


for kps=1:Kf;
 kps 
datk=textread([datdir,horzcat(num2str(kps),'.txt')],'','headerlines',8);

            %datk=datk(40e3:end,:); disp('exclude first 2 sec of data')
            %datk(:,2:3)=(1/20).*datk(:,2:3); disp('Vpd/20')
            ik=1:numel(datk(:,chxpd));%80e3;%find(caldat(:,1)<=Tmsr);
            fsamp=1e6/(Samp_rate(1,kps));
       if cht==0,
          ctime=ik'./fsamp;
       else
         ctime=caldat(ik,cht);
       end   
       
            lenghtdata=size(datk,1);
            fsamp=1e6/(Samp_rate(1,kps));
            fexc=freq_list(kps);
            Nfft_pre=2^(nextpow2(lenghtdata)-1);%fsamp*2;%fsamp;%floor(fsamp/4);%10e3; % it was 
            Nfft_pre2=floor(Nfft_pre/fsamp)*fsamp;
            if fexc==0
                Nfft= Nfft_pre2;
            else
            Nfft=floor(Nfft_pre2/fexc)*fexc;
            end
            
            Noverlap=floor(Nfft/2);%5e3;
            Nwindow=Nfft;
            Nint=10;
            try
            [psx_k,psdx_k]=power_spectrum([ctime,kpd.*datk(:,chxpd)]',Nwindow,Noverlap,Nfft,Nint,fsamp);
            
            
            FP=psx_k(:,1);
           % PYY(:,kps)=psx_k(:,2);
            fx=FP;
            jfit2=find((fx>fmin & fx<fmax) & abs(fx-Fspike(1))>50 & abs(fx-Fspike(2))>30);
            FP=[];
            FP=fx(jfit2);%(Fexc+10):10:7e3;
            A=psx_k(:,2);
            PYY=[];
            PYY(:,1)=A(jfit2);
            PYm=PYY;
            figure, loglog(FP,PYm,'b','linewidth',2), hps=gcf;
            xlabel('f (Hz)'), ylabel('Power Spectrum ($V^2*s$)')
            catch
                'powererror'
                figure
            end

end
     
     
 selectedset=input('Which power looks good (enter the figure number)?') ;   
 close all;          
 
%  if ~isempty(if_freqzero) 
%     Kf
%     Kf=Kf-1
%  end

 %SS=1
 
 

        for kf=1:Kf-2%Kf%  Kf,
            kf
            fexc=freq_list(kf);
            if fexc==0
                continue
            end
            
            
            feps=fexc/10;   %%%???  allowed variation in excitation 
            
            datk=textread([datdir,horzcat(num2str(kf),'.txt')],'','headerlines',8);  %%%%hossein change this!!! it was kf originally. I added the one. move it
             
            
            [mk,nk]=size(datk);
            if AODorstage(kf)==1
            xtrap=cx.*(datk(:,checkAOD)-75);   
            elseif AODorstage(kf)==0
            xtrap=cxstage.*datk(:,checkstage);
            end
            
            vpd=kpd.*datk(:,chxpd);  %%% ?????????????????????????????????   QPD response
            
           
            lenghtdata=size(datk,1);
            fsamp=1e6/(Samp_rate(1,kf));
            fexc=freq_list(kf);
            %fex_cor=0;
            padcoef=1;
            n=size(xtrap,1)*padcoef;
            matzero=zeros(n,1);
            xtrap2=[xtrap;matzero];
            vpd2=[vpd;matzero];
            
            
            for i=1:1
                if fexc==0
                    continue
                end
            fexc_cor(kf)=fexc;
            
            Nfft_pre=2^(nextpow2(lenghtdata)-1);%fsamp*2;%fsamp;%floor(fsamp/4);%10e3; % it was 
            Nfft_pre2=floor(Nfft_pre/fsamp)*fsamp;
            Nfft=floor(floor(Nfft_pre2/fexc)*fexc);
            Noverlap=floor(Nfft/2);%5e3;
            %Nwindow=Nfft;
            Nint=10;
            
            try
            [P,F]=spectrum(xtrap,vpd,Nfft,Noverlap,hanning(Nfft),fsamp,0.95);
            catch
                'error used alternative'
               Nfft= Nfft_pre2;
               Noverlap=floor(Nfft/2);
               [P,F]=spectrum(xtrap,vpd,Nfft,Noverlap,hanning(Nfft),fsamp,0.95);
            end
            %[P,F]=spectrum(xtrap,vpd,mk,Noverlap,hanning(Nfft),fsamp,0.99);
            if_range=find(abs(F-fexc)<feps);
            if isempty(if_range)
                if_range=find(abs(F-fexc)<0.1);
            end
            
            if_des=if_range(floor(median(find(P(if_range,5)==max(P(if_range,5))))));
            fexc= F(if_des);
            
            end
            

            Fexc(kf)=F(if_des);
            Pxx(kf)=P(if_des,1);
            Pyy(kf)=P(if_des,2);
            H(kf)=P(if_des,4);
            C(kf)=P(if_des,5);
            
            if2_des=find(abs(F-2*Fexc(kf))==min(abs(F-2*Fexc(kf))));
            Pyy2(kf)=P(if2_des,2);
            
            
         
        end
       
        
        
        Mfrf=20.*log10(abs(H));
        PHfrf=unwrap(angle(H));
       % js=find((datk(:,checkAOD)-75)==4);                    %      js=find(chexc_list==4);  this to crroect the phase shift for stage                             %???????????????????????????????????????????????? Ask adam
       % PHfrf(js)=PHfrf(js)-pi;                              %???????????????????????????/???/????  Ask Adam
       PHfrf =(PHfrf - round(PHfrf./pi).*pi);% ???
        
        % calculate power spectrum %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for kps=selectedset
            datk=textread([datdir,horzcat(num2str(kps),'.txt')],'','headerlines',8);
            
            ik=1:numel(datk(:,chxpd));%80e3;%find(caldat(:,1)<=Tmsr);
            fsamp=1e6/(Samp_rate(1,kps));
       if cht==0,
          ctime=ik'./fsamp;
       else
         ctime=caldat(ik,cht);
       end
       
            
            
            lenghtdata=size(datk,1);
            fsamp=1e6/(Samp_rate(1,kps));
            fexc=freq_list(kps);
            Nfft_pre=2^(nextpow2(lenghtdata)-1);%fsamp*2;%fsamp;%floor(fsamp/4);%10e3; % it was 
            Nfft_pre2=floor(Nfft_pre/fsamp)*fsamp;
           if fexc==0
               Nfft=Nfft_pre2
           else
            Nfft=floor(Nfft_pre2/fexc)*fexc;
           end
            Noverlap=floor(Nfft/2);%5e3;
            Nwindow=Nfft;
            Nint=10;
            
            
       [psx_k,psdx_k]=power_spectrum([ctime,kpd.*datk(:,chxpd)]',Nwindow,Noverlap,Nfft,Nint,fsamp);
            FP=psx_k(:,1);
           % PYY(:,kps)=psx_k(:,2);
            fx=FP;
            jfit2=find((fx>fmin & fx<fmax) & abs(fx-Fspike(1))>50 & abs(fx-Fspike(2))>30);
            FP=[];
            FP=fx(jfit2);%(Fexc+10):10:7e3;
            A=psx_k(:,2);
            PYY=[];
            PYY(:,1)=A(jfit2);
        end
        
        PYm=mean(PYY,2);
%loglog(FP,PYY)
        Fi=FP(1:(end-1));
        Pyi=cumtrapz(FP,PYm);

        %save(fname_save,'Fexc','Mfrf','PHfrf','C','H','FP','PYm','Pyy','Pyy2')
jkp=[];
for kkp=1:(numel(FRF_f)/2);
    jfrf=2*(kkp-1)+[1,2];
    jkkp=find(Fexc>=FRF_f(jfrf(1)) & Fexc<=FRF_f(jfrf(2))&Fexc>0);
    jkp=[jkp,jkkp];
end
%file_list=file_list(jkp);
        %freq_list=freq_list(jkp);
         %cx_list=cx_list(jkp);
              %chexc_list=chexc_list(jkp);
Fexc=Fexc(jkp);
Mfrf=Mfrf(jkp);
PHfrf=PHfrf(jkp);
C=C(jkp);
H=H(jkp);
Pyy=Pyy(jkp);
Pyy2=Pyy2(jkp);

% wrap phase
%ph_off=round(mean(PHfrf./180));
%PHfrf=PHfrf-ph_off*180;


% plot spectra ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
figure, subplot(211), semilogx(Fexc,Mfrf,'o','linewidth',2)
ylabel('Mag, dB'), set(gca,'Xlim',[0.01 5e3]);
subplot(212), semilogx(Fexc,PHfrf,'o','linewidth',2)
ylabel('Phase, deg'), xlabel('f, Hz'), set(gca,'Xlim',[0.01 5e3]);
hfrf=gcf;

% figure, subplot(121), semilogx(Fexc,real(H),'o')
% xlabel('f (Hz)'), ylabel('Real(H)')
% subplot(122), semilogx(Fexc,imag(H),'o')
% xlabel('f (Hz)'), ylabel('Imag(H)')

figure, semilogx(Fexc,C,'o','linewidth',2)
ylabel('Coherence'), xlabel('f, Hz')

figure(hfrf), subplot(211),
% if zve==0|zve==2,
%     ig.fz_frf=0;
%     ig.Pz=0;
% elseif zve==1,
%     disp('Click on zero corner frequeny (Hz)')
%     [ig.fz_frf,ig.Pz]=ginput(1);
% end
% disp('Click on pole corner frequency (Hz)')
% [ig.fp_frf,ig.MinfdB]=ginput(1);

figure, loglog(FP,PYm,'b','linewidth',2), hps=gcf;
xlabel('f (Hz)'), ylabel('Power Spectrum ($V^2*s$)')

%disp('Click on corner frequency (Hz) and steady-state amplitude ($V^2$)')
%[ig.fp_ps,ig.Pss]=ginput(1);
%ig.Dfit=ig.Pss*pi^2*(400+ig.fp_frf^2);


% fit frequency response and power spectrum  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ig.Minf=10^(ig.MinfdB/20);
% ig.Gamma0=9.42e-6; %[pN*s/nm]
% %disp('Must fix: Calculate Gamma0(T)')
% 
% % find slope of power spectrum at high F %%%
% disp('Click on min and max freq. to calculate slope')
% [slope_psmin,ymin]=ginput(1);
% [slope_psmax,ymax]=ginput(1);  %%% hoddrin changed this
% jps=find(FP>=slope_psmin & FP<=slope_psmax);
% pslope=polyfit(log(2*pi.*FP(jps)),log(PYm(jps)),1);
% ig.Phi=-pslope(1);
% hold on, loglog(FP,exp(polyval(pslope,log(2*pi.*FP))),'r')
% 
% ig.Pss=exp(polyval(pslope,log(2*pi*ig.fp_frf)));
% ig.Dfit=ig.Pss*pi^2*((3*ig.fp_frf)^2);

% Hth=tf(ig.Minf*[1 2*pi*ig.fz_frf],[1 2*pi*ig.fp_frf]);
% [Hm,Hph,Hw]=bode(Hth);
%figure(hfrf), subplot(211), hold on,
%semilogx(Hw./(2*pi),20.*log10(Hm(:)),'g')
%subplot(212), hold on,
%semilogx(Hw./(2*pi),Hph(:),'g')

% ig.Beta=1/ig.Minf;
% ig.GammaR=(kT/(ig.Beta^2*ig.Dfit))/ig.Gamma0;
% % ig.Alpha=0.75;%max([0,ig.Phi-1]);%1-ig.Phi/2;
% ig.Alpha=abs(pslope(1))-1;
% 
% wc=2*pi*ig.fp_frf;
% %ig.k_cyt1=max([0.001,(ig.GammaR*ig.Gamma0*wc-ig.k_trap-ig.k_cyt0)/(wc^ig.Alpha)]);
% ig.k_cyt1=5e-3;;%0.01;%max([0.001,sqrt(4*kT*ig.GammaR*ig.Gamma0/ig.Pss)-ig.k_trap-ig.k_cyt0]);
% ig.k_cyt0=ig.fz_frf*2*pi*ig.GammaR*ig.Gamma0;
% ig.k_trap=ig.fp_frf*2*pi*ig.GammaR*ig.Gamma0-ig.k_cyt0;%0.05;
% ig.keq0=ig.k_trap+ig.k_cyt0+ig.k_cyt1*(ig.fp_frf)^ig.Alpha;

%ig.k_trap=0.07;
%ig.GammaR=(ig.k_trap+0.01+ig.k_cyt1)/(2*pi*ig.fp_frf*ig.Gamma0);

% ig.k_cyt1=5e-3;
% ig.k_cyt0=5e-2;
% ig.k_trap=0.15;
% ig.Beta=25;
% ig.GammaR=30;
%ig.Alpha=abs(pslope(1))-1;


%     %   Beta    GammaR    Alpha    k_trap   k_cyt0   k_cyt1 mass    nu];
%     lb=[0.1.*ig.Beta   0.5*ig.GammaR     ig.Alpha-0.01        0         0.1*ig.k_cyt0    0      1e-10    0];
%     ub=[10*ig.Beta     100*ig.GammaR     ig.Alpha+0.01        ig.k_trap+0.1         10.*ig.k_cyt0    0.1    1e10     1.01];

%   Beta    GammaR    Alpha    k_trap   k_cyt0   k_cyt1 mass    nu];
% lb=[0.1.*ig.Beta   0.01*ig.GammaR    0.5       0.005     0.001        0.001       1e-1    1e-6];
% %ub=[10*ig.Beta     10*ig.GammaR     0.9       0.2   10.*ig.k_cyt0    1e-1    1e1    1.01];
% ub=[10*ig.Beta     10*ig.GammaR     ig.Alpha+0.05       0.5   10.*ig.k_cyt0    1e-1    1e1    1.01];
% 
% %ub=[1e3     1e4     5e1     100     2           falias+1500 1e10     1e8];
% % else
% %     u10=[ig.Dfit Qo*1e6 f0];
% %     par_range=[0.2*Dfit 0.1*fc Qo f0];
% %     lb=[0,100,0,0];
% %     ub=[50,2000,10,150];
% % end
%disp([' - fit using ',num2str(numel(u10)),' parameters -'])
%options=optimset('TolFun',1e-12,'TolX',1e-10,'Display','on','MaxFunEvals',1e14,'MaxIter',1e4);%,...
%'LargeScale','off','LevenbergMarquardt','on');%,'DiffMinChange',0.1);%,'DiffMinChange',5);%,'MaxFunEvals',1e3,'Display','iter');
% jfit1=find((FP>fmin & FP<fmax) & PYm>0 & isnan(PYm)~=1);
% fdata1=FP(jfit1);%(Fexc+10):10:7e3;
% ydata1=cumtrapz(FP(jfit1),PYm(jfit1));%interp1(fx,IPSx,fdata,'pchip')
% ydat10=ydata1;
% 
% figure, semilogx(fdata1,ydata1,'b','linewidth',2), hips=gcf; %,fdata1,IPth1,'g','linewidth',1
% xlabel('f (Hz)'), ylabel('$\int P_{th}$ ($V^2$)')
% axis tight
% 
% 
% %u0=[ig.Beta, ig.GammaR, ig.Alpha, ig.k_trap, ig.k_cyt0, ig.k_cyt1, 10, 1]; %[Minf,fcp,fcz,Dfit,phi,falias,mass (*1e-21 g/nm^3),nu (*1e12 nm^2/s)]
% u0=[ig.Beta, ig.GammaR, ig.Alpha, ig.k_trap, ig.k_cyt0, ig.k_cyt1, 10, 1]
% %lb=[];
% %ub=[];
% % options=optimset('TolFun',1e-10,'TolX',1e-8,'Display','iter','MaxFunEvals',1e5,'MaxIter',2.5e3);%,'DiffMinChange',1);
% % options=optimset(options,'LargeScale','off','DiffMinChange',1e-1,'DiffMaxChange',10);
% % %options2=options;%optimset(options,'DiffMinChange',1e-3,'DiffMaxChange',1e-1);
% % options2=optimset(options,'LargeScale','off','DiffMinChange',1e-3,'DiffMaxChange',1,'TolFun',1e-10,'TolX',1e-9);

% options=optimset('TolFun',1e-8,'TolX',1e-6,'Display','iter','MaxFunEvals',1e5,'MaxIter',2e3);%,'DiffMinChange',1);
% %options=optimset(options,'Largescale','on');
% options=optimset(options,'Algorithm',{'levenberg-marquardt',0.01})
% %options=optimset(options,'LargeScale','off','DiffMinChange',1e0,'DiffMaxChange',1e1);
% options2=optimset(options,'TolFun',1e-8,'TolX',1e-8,'DiffMinChange',1e-4,'DiffMaxChange',1e2,'MaxIter',5e3,'Largescale','off');
% %options2=optimset(options,'LargeScale','off','DiffMinChange',1e-3,'DiffMaxChange',1,'TolFun',1e-10,'TolX',1e-9);
% 
% 
% Wfrf=Wfrf;%*(numel(Fexc)./numel(fdata1));
% ff=1:10:2e4;
% if zve==0,
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
% elseif zve==1,
%     if zph==0,
%         PHfrf=[];
%         xdata=[Fexc,fdata1'];
%         ydata=[zeros(size(H)),zeros(size(ydat10'))];
%     else
%         xdata=[Fexc,Fexc,fdata1'];
%         ydata=[zeros(size(H)),zeros(size(H)),zeros(size(ydat10'))];
%     end
%     W=Wfrf.*C(:)';%.^25;
%     [u2i,resnorm,residual]=lsqcurvefit(@theor_trap_frf_r6,u0,xdata,ydata,lb,ub,options,abs(H(:))',PHfrf,W,ydat10',fsamp,Rbead,Zbead,kT,falias);
% 
%% 
%     [mf0,phf0]=frf_r6(u0,Fexc,Rbead,Zbead,kT,falias);
%     psf0=ps_r6(u0,ff,fsamp,Rbead,Zbead,kT,falias);
%     figure(hfrf), subplot(211), hold on,
%     semilogx(Fexc,20.*log10(mf0),'k','linewidth',2)
%     subplot(212), hold on,
%     semilogx(Fexc,phf0,'k','linewidth',2)%+(180/pi).*unwrap(angle(exp(-j*1e-6.*Hw)))
%     figure(hps), hold on
%     loglog(ff,psf0,'k','linewidth',2)
%     
%     [mf1,phf1]=frf_r6(u2i,Fexc,Rbead,Zbead,kT,falias);
%     psf1=ps_r6(u2i,ff,fsamp,Rbead,Zbead,kT,falias);
%     figure(hfrf), subplot(211), hold on,
%     semilogx(Fexc,20.*log10(mf1),'g','linewidth',2)
%     subplot(212), hold on,
%     semilogx(Fexc,phf1,'g','linewidth',2)%+(180/pi).*unwrap(angle(exp(-j*1e-6.*Hw)))
%     figure(hps), hold on
%     loglog(ff,psf1,'g','linewidth',2)
% 
%     %     lb2=[0.7*u2i(1)     0.7*u2i(2)      0.8*u2i(3)        0.8*u2i(4)         0.8*u2i(5)      0      0.1    0];
%     %     ub2=[1.3*u2i(1)     1.3*u2i(2)      1.2*u2i(3)        1.2*u2i(4)         1.2*u2i(5)      1      1e3     1.01];
%     %     lb2=[10     1      u2i(3)-0.01        0        0      0      0.1    0];
%     %     ub2=[200    1500   u2i(3)+0.01        1        1      1      1e2    1.01];
%     lb2=[10     0.5   0.5       1e-5        1e-4        1e-6        0.1    1e-4];
%     %ub2=[200    1000   0.9       0.3      0.3      1e-1     1e1    1.01];
%     ub2=[200    1000   0.9      1e-2      0.03      1e-1     1e1    1.01];
%     W2=Wfrf.*C(:)';
    
%     [u2,resnorm,residual,exitflag,output,lambda,jacobian]=lsqcurvefit(@theor_trap_frf_r6,u2i,xdata,ydata,lb2,ub2,options2,abs(H(:))',PHfrf,W2,ydat10',fsamp,Rbead,Zbead,kT,falias);
%     
% %     disp('!!! using nlinfit !!!')
% %     %theor_fun=@(u,fdata) theor_trap_frf_r6(u,fdata,abs(H(:))',PHfrf,W2,ydat10',fsamp,Rbead,Zbead,kT,falias);
% %     %[u2,residual,Jacobian,COVB,mse]=nlinfit(xdata,ydata,theor_fun,u2i,options2);
% %     options2=statset(options2,'Robust','on');
% %      [u2,residual,jacobian,COVB,mse]=nlinfit(xdata',ydata',...
% %         @(u,fdata) theor_trap_frf_r6(u',fdata',abs(H(:))',PHfrf,W2,ydat10',fsamp,Rbead,Zbead,kT,falias)',...
% %         u2i',options2);
% %     
%     u2_ci=nlparci(u2,residual,'jacobian',full(jacobian)); %nlparci to calculate confidence intervals on the fit
%     
%     errfit=theor_trap_frf_r6(u2,xdata,abs(H(:))',PHfrf,W,ydat10',fsamp,Rbead,Zbead,kT,falias);
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
% 
%     [Mfrf_fit,PHfrf_fit]=frf_r6(u2,Fexc,Rbead,Zbead,kT,falias);
%     %Mfrf_fit=((errfit(jfrf)./W)+1).*abs(H(:))';
%     IPS_fit=(errfit(jps)+1).*ydat10';
%     PS_fit=ps_r6(u2,ff,fsamp,Rbead,Zbead,kT,falias);
% 
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
% end

% Hth=tf(u2(1)*[1 2*pi*u2(3)],[1 2*pi*u2(2)]);
% Hpd=tf(2*pi*u2(end),[1 2*pi*u2(end)]);
% Hth=series(Hth,Hpd);
% Hw=2*pi.*Fexc;
% [Hm,Hph]=bode(Hth,Hw);
% figure(hfrf), subplot(211), hold on,
% semilogx(Fexc,20.*log10(Mfrf_fit),'c','linewidth',2)
% subplot(212), hold on,
% semilogx(Fexc,PHfrf_fit,'c','linewidth',2)%+(180/pi).*unwrap(angle(exp(-j*1e-6.*Hw)))
% 
% 
% figure(hps), hold on
% loglog(ff,PS_fit,'c','linewidth',2)
% 
% % IPth1=cumtrapz(fdata1,Pth_t0);%excited_power_spectrum(u,fx(jfit),falias,Fexc,Adrive1,ydat0,f0);
% % figure(hips), hold on, semilogx(fdata1,IPth1,'c','linewidth',2)
% figure(hips), hold on, semilogx(fdata1,IPS_fit,'c','linewidth',2)
% 
% P12=20.*log10(Pyy)-20.*log10(Pyy2);
% figure, semilogx(Fexc,P12,'bo')
% xlabel('Frequency (Hz)'), ylabel('Pyy(Fexc) - Pyy(2*Fexc)')

% calculate trap parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   u1=[Dfit phi] from power spectrum
%   u2=[Minf fc feq] from frequency response
%Dfit=u1(1);
%fc_ps=u1(2);
%fc_frf=u2(2);
%fc_zero=u2(3);
%Beta=1/u2(1);
%Gamma=kT/(Beta^2*Dfit);
%k_cyt=2*pi*fc_zero*Gamma;
%k_trap=2*pi*(fc_frf-fc_zero)*Gamma;
%m=u2(end-1); %*10^-21 g
%nu=u2(end); %*10^12 nm^2/s
% beta=u2(1);
% gammar=u2(2);
% alpha=u2(3);
% k_trap=u2(4);
% k_cyt0=u2(5);
% %k_eq0=u2(6);
% %w_0=k_eq0/(gammar*9.42e-6);
% k_cyt1=u2(6);%(k_eq0-k_trap-k_cyt0)/(w_0^alpha);
% 
% %disp(['Dfit = ', num2str(Dfit), ' V'])
% %disp(['f_c,zero = ', num2str(fc_zero), ' Hz'])
% %disp(['f_c,FRF = ', num2str(fc_frf), ' Hz'])
% disp(['Beta = ', num2str(beta), ' nm/V'])
% disp(['Gamma = ',num2str(gammar),' *gamma0'])
% disp(['Alpha = ',num2str(alpha)])
% disp(['k_trap = ',num2str(k_trap),' pN/nm'])
% disp(['k_cyt,0 = ', num2str(k_cyt0), ' pN/nm'])
% disp(['k_cyt,1 = ',num2str(k_cyt1),' pN/nm'])
% disp(['Beta*k_trap = ',num2str(beta*k_trap),' pN/V'])

% k_cyt_tot=k_cyt0+k_cyt1.*(2*pi.*FP).^alpha;
% Gp=k_cyt_tot./(6*pi*Rbead);
% figure, loglog(FP,Gp,'linewidth',2)
% xlabel('f (Hz)'), ylabel('Gp ($pN/nm^2$)')
% 
% ww=2*pi.*FP;
% g=k_cyt0 + k_cyt1.*((j.*ww).^alpha)./gamma(alpha)...
%     + gammar*ig.Gamma0.*j.*ww;
% Gp=g./(6*pi.*Rbead);
% v=imag(g)./ww;
% 
% figure, hold on
% loglog(FP,real(Gp),'color','b','linewidth',2)
% ylabel('Gpp - Storage Modulus')
% set(gca,'Xscale','log'), set(gca,'Yscale','log')

% figure, hold on
% loglog(FP,imag(Gp),'color','b','linewidth',2)
% ylabel('Gp - Loss Modulus')
% set(gca,'Xscale','log'), set(gca,'Yscale','log')
% 
% figure, hold on
% loglog(FP,v,'color','b','linewidth',2)
% ylabel('viscosity')
% set(gca,'Xscale','log'), set(gca,'Yscale','log')
% 
% 
% fname_save_results=[fname_save(1:end-10),'r6c_cal.mat'];
% zsr=input('Save calibration results (0/1)? ');
% if zsr==1
%     save(fname_save_results)
% end


cd(currentFolder)

% fname_save_results=[fname_save(1:end-10),'r6c_cal.mat'];
% zsr=input('Save calibration results also to the files directory(0/1)? ');
% if zsr==1
%     save(fname_save_results)
% end

% fname_save_results=[fname_save(1:end-10),'r6c_cal-frfdat.mat'];
% zsr=input('Save calibration results also to the files directory(0/1)? ');
% if zsr==1
%     save(fname_save_results)
% end
cd(currentFolder)
savename=horzcat('noloop_2fischerNfftconstant_transferfunctionworkspace_QPDchanel_',num2str(chxpd));
save(savename)
   
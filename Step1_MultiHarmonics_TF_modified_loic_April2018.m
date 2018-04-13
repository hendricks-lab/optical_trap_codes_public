%%% Modified by Loïc Chaubet November 10th 2017 %%%
%/Users/loicchaubet/OneDrive - McGill University/DATA/Hossein/Loic/Blebb/Sept_22_added_Set3/C193_P7/p4_great_analysis/others
%0.4set3_maybe.txt

clear
close all
clc

datdir=[pwd '/'];
cd(datdir);
%%% enter the file name here, make sure the current matlab working folder
%%% is where this file is
datk=textread([datdir,'180s_0.4_CONTROL_not_bad.txt'],'','headerlines',0); % was set to ignore first 8 rows...now changed to ignore 0 (Dec 21st 2017)
fname_save='Result_';
pwd
currentFolder = pwd
 

%%%%% PARAMETERS FOR ADAM'S METHOD FIT (NOT USED IN THIS CODE)  %%%%%%%%%%%
zve=1; % 0 no viscoelastic; 1 viscoelastic
Temp=37;
kT=1.381e-2*(Temp+273.15);
zph = 0;
%if_freqzero=find(freq_list==0);
pixelsize=105; %????????
Wfrf = 0.1;   %0.2;%0.05;%0.2  needs to be adjusted. weighting funtion when fitting frequency and power spectrum data . Different number of points in data sets. 
falias=3e5;  % QPD has filtering effect. first order low pass filter ????????? NEED TO CHECK THIS?
Qo=1e-5;% source of vibration at low frequency 
f0=0;%1;%20;%5; %[Hz] building vibrations (0 if none)
Rbead=250; %Radius of bead (nm)
Zbead=2500; %height of bead (nm)
viscosity0=1e-9; %[pN*s/nm^2]
kp=0;         %%%%????? 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%% Below is the frequencies of excitation exerted by the AOD during 
%%% experiment and the corresponding integer multiple that gives best FRF 
%%% results. The amplitude of the AOD is shown as a reference, but not used 
%%% in the code.

%%%%% SET 1.0xFreq with 17 frequencies (24 November 2017)
%freq_list =         [0.023   0.053   0.17   0.43   1.37    2.37    5      7     19     29     61     87     147    337      899      1347      2441];
% AOD amp (for ref.) [1.4     1.2     0.9    0.8    0.4     0.4    0.29   0.29   0.25   0.24   0.21   0.2    0.13   0.09    0.03      0.02       0.02];
%K_integerOPT =      [3       4       5      19     50      79     182    225    800    1050   2250   2950   4500   13000   40000     53000     99000]; %C193_P5_p2_Nov_24_2017_0.4_1.0x_17freq_FINAL and C1_P6_p2_Nov_24_2017_FINAL
%K_integerOPT =      [2       3       7      22     65     130     240    280   1050   1250   2300   3100   5500   15500   55000     72000    120000];
%K_integerOPT =      [2       3       6      19     53      111     202    234   950   1230   2250   2950   4500   13000   40000     53000     99000];
%K_integerOPT =      [1       2       6      19     53      90     202    234   950   1230   2250   2950   4500   13000   40000     53000     99000];


%%%%% SET 1.2xFreq with 17 frequencies (24 November 2017)
% freq_list =         [0.028   0.064   0.21  0.52   1.64    2.84     6      8.3     23     35     73    103    177    403     1079      1617      2929];
%  AOD amp (for ref.) [1.4     1.2     0.9    0.8    0.4     0.4    0.29   0.29    0.25   0.24   0.21   0.2    0.13   0.09    0.03      0.02      0.02];
%K_integerOPT =      [3       6       16      33     90     150     300    420    1280   1705   3900   5300   9500   24000   65000     88000     200000];
%K_integerOPT = floor(1.1*[3      6       13      33     90     150     300    420    1280   1705   3900   5300   9500   24000   65000     88000     200000]);
%K_integerOPT = floor(2*[2       3       6      19     50      79     182    225    800    1050   2250   2950   4500   13000   40000     53000     99000]);
% HARMONICS at 56 and 131 

%%%%% SET 1.0xFreq with 17 frequencies (29 November 2017 and Dec 12)
%freq_list =         [0.023   0.053   0.17   0.43   1.37    2.37    5      7     19     27     67     87     157    337      899      1347      2441];
%AOD amp (for ref.)  [1.4     1.2     0.9    0.8    0.4     0.4    0.29   0.29   0.25   0.24   0.21   0.2    0.12   0.09    0.03      0.02       0.02];
%K_integerOPT =      [1       2       5      19     50      79     182    225    800    1050   2350   2950   4700   13000   40000     53000     99000];
%K_integerOPT =      [3       4       6      20     50      79     182    225    800    1050   2350   2950   4700   13000   40000     53000     99000];%C193_p2_Nov_29_0.4_1.0x (FINAL)
%K_integerOPT =      floor(1.4*[2       5       6      19     50      79 182    225    800    1050   2350   2950   4700   13000   40000     53000 99000]);%C1_p1_Nov_29_0.3_1.0x (FINAL)
%K_integerOPT =     floor(1.7*[2       3       6      20     50      79     182    225    800    1050   2350   3300   5100   13000   40000     53000     99000]);%C193_P6_p1_Nov_29_2017_0.3_1.0x_17freq_FINAL
%K_integerOPT =     floor(1.9*[2       4       6      20     50      79     182    225    800    1050   2350   3300   5400   13000   40000     53000     99000]);%C193_P6_p1_Nov_29_2017_0.4_1.0x_17freq_FINAL
%K_integerOPT =      [3       5       8      20     50      79     182    225    800    1050   2350   2950   4700   13000   40000     53000     99000];%C1_p1_Dec_12_0.3_1.0x (FINAL), C1_P11_p2_Dec_12_2017_0.4_1.0x(FINAL)
%K_integerOPT =      [3       6       10      20     50      79     182    225    800    1050   2350   2950   4700   13000   40000     53000     99000];
%K_integerOPT =      [3       6       8      20     50      79     182    225    800    1050   2350   2950   4700   13000   40000     53000     99000];%CONTROL_noheater_C193_P11_p4_0.2_1.0x_Dec_12_2017_FINAL
%K_integerOPT =      floor(1.8*[2       4       8      20     50      79     182    225    800    1050   2350   2950   4700   13000   40000     53000     99000]);
%K_integerOPT =       3*[3       5       8      20     50      79     182    225    800    1050   2350   2950   4700   13000   40000     53000     99000];
%K_integerOPT =      [3       5       15      29     75      92     182    225    800    1050   2350   2950   4700   13000   40000     53000     99000]; %Jan 30 C1_P8 p2 
%K_integerOPT =      floor(1.3*[3       6       18      29     75      92     182    255    800    1050   2350   2950   4700   13000   40000     53000     99000]); %Feb 13 C1 CONTROL
%K_integerOPT =      floor(1.4*[4       7       15      29     75      92     182    255    800    1050   2350   2950   4700   13000   40000     53000     99000]);
%K_integerOPT =      [2       5       21      29     75      92     182    225    800    1050   2350   2950   4700   13000   40000     53000     99000]; %Jan 30 C193_P7 p4
%K_integerOPT =      floor(1.3*[2       3       20      29     75      92     182    225    800    1050   2350   2950   4700   13000   40000     53000     99000]);
%K_integerOPT =      floor(0.8*[4       9       15      29     75      92     182    225    800    1050   2350   2950   4700   13000   40000     53000     99000]);


%%%%% SET 1.2xFreq with 17 frequencies (29 November 2017 and Dec 12)
%freq_list =         [0.027   0.063   0.21  0.51   1.63    2.83     6      8.3     23     31     81    103    188    409     1079      1617      2929];
% AOD amp (for ref.) [1.4     1.2     0.9    0.8    0.4     0.4    0.29   0.29    0.25   0.24   0.21   0.2    0.13   0.09    0.03      0.02      0.02];
%K_integerOPT =      [3       6       16      33     90     150     300    420    1280   1705   3900   5300   9500   24000   65000     88000     200000];%C193_p2_Nov_29_0.4_1.2x (FINAL)
%K_integerOPT =      [3       5       17      33     120     210     450 490    1750   2350   6300   6150   11100   26000   66000     100000 200000];%C1_p1_Nov_29_0.3_1.2x (FINAL)
%K_integerOPT =      floor(1.2*[3      5       17      49     120     210     450  540    1750   2750   6600   6800   15100   29000   66000     100000 200000]);
%K_integerOPT = floor(1.4*K_integerOPT);
%K_integerOPT =      [3      7       17      33     165     210     450    490    1750   2350   6300   6150   11100   26000   66000     100000     200000];%C1_p1_Dec_12_0.3_1.2x (FINAL), C1_P11_p2_Dec_12_2017_0.4_1.2x (FINAL)
%K_integerOPT =      [4      7       24      33     165     210     450    490    1750   2350   6300   6150   11100   26000   66000     100000     200000];%%Jan 30 C1_P8 p2

% HARMONICS at 56 and 131, smaller amplitude...

%WATER JAN 9 1.0x 20 FREQ
%freq_list =                   [0.023   0.037   0.053   0.083   0.11    0.17    0.43   1.37    2.37     5       7      19       27      67    87     157     337       899      1347      2441];
%AOD amp (for ref.)            [1.4     1.3     1.2     1.1     1.0     0.9     0.8    0.4     0.4     0.29    0.29    0.25    0.24   0.21    0.2    0.12    0.09      0.03     0.02      0.02];
%K_integerOPT =      floor(1.5*[2       3        5       6       8       10      27     88      93      220     275    890     1200   2750    3350   5800    15000     44000    57000    103000]);

%WATER JAN 9 1.0x 20 FREQ

%%%%% SET 1.0xFreq with 17 frequencies (MARCH 1st 2018!!!!)
freq_list =         [0.023   0.053   0.097  0.17   0.43   1.37    2.37    5      9.1     19     37     87     157    337      899      1347      2441];
%AOD amp (for ref.) [1.4     1.2     0.95   0.9    0.8    0.4     0.4    0.29   0.28   0.25   0.22     0.2    0.12   0.09    0.03      0.02       0.02];
K_integerOPT =      floor(1.5*[3       5       11      16      37     92      96     182    400    800    1250    2950   4900   13000   40000     53000     99000]); 
%K_integerOPT =      floor(2.0*[2       3       6      10      25     92      96     182    400    800    1250    2950   4900   13000   40000     53000     99000]); 


%%% Some other parameters 
coefficent = 2; % used for the power spectrum, may need to be ajusted
Fexc = freq_list;
Fspike = [5000 6000];
Samp_rate = 50; % [us]
cht = 0; % there is no time channel
FRF_f = [0,5e3]; %  range to fit frequency response
kpd = 1;  % correction factor for pd, leave it one       ????? ask adam what this is?
fmin = 3e-3; 
fmax = 1e4;
AODchannelx = 5; % this is the column where input AOD is written in the text file
chxpd = 2; % this is the associed column for X-QPD data
stagechannel = 4; % this is input channel for excitation of stage  
cx = 1750; % conversion factor from AOD MHz to nm (calibrated)
cxstage = 20000; % conversion factor for stage (used if using stage oscillation)
Kf = length(freq_list);          



%%%%%% POWER SPECTRUM SECTION %%%%%
            %datk=datk(40e3:end,:); disp('exclude first 2 sec of data')
            %datk(:,2:3)=(1/20).*datk(:,2:3); disp('Vpd/20')
            ik=1:numel(datk(:,chxpd));%80e3;%find(caldat(:,1)<=Tmsr);
            fsamp=1e6/Samp_rate;
       if cht==0,
          ctime=ik'./fsamp;
       else
         ctime=caldat(ik,cht);
       end   
       
            lengthdata=size(datk,1);
            fsamp=1e6/Samp_rate;
            fexc=0;
            Nfft_pre=2^(nextpow2(lengthdata)-1);%fsamp*2;%fsamp;%floor(fsamp/4);%10e3; % it was 
            Nfft_pre2=floor(Nfft_pre/fsamp)*fsamp;
            if fexc==0
                Nfft= Nfft_pre2*coefficent;
            else
            Nfft=floor(Nfft_pre2/fexc)*fexc*coefficent;
            end
            
            Noverlap=floor(Nfft_pre/2);%5e3;
            Nwindow=floor(Nfft_pre);
            Nint=10;
    
           [psx_k,psdx_k]=power_spectrum([ctime,kpd.*datk(:,chxpd)]',Nwindow,Noverlap,Nfft,Nint,fsamp);
            
            
            FP=psx_k(:,1);
           %PYY(:,kps)=psx_k(:,2);
            fx=FP;
            jfit2=find((fx>fmin & fx<fmax) & abs(fx-Fspike(1))>5 & abs(fx-Fspike(2))>5);
            FP=[];
            FP=fx(jfit2);%(Fexc+10):10:7e3;
            A=psx_k(:,2);
            PYY=[];
            PYY(:,1)=A(jfit2);
            PYm=PYY;
            figure, loglog(FP,PYm,'b','linewidth',2), hps=gcf;
            title('X-QUAD signal (output)')
            xlabel('f (Hz)'), ylabel('Power Spectrum ($V^2*s$)')
      
 
            
%%%%  TRANSFER FUNCTION SECTION %%%%
for kf=1:Kf
    kf
    [mk,nk]=size(datk);
    xtrap=cx.*(datk(:,AODchannelx)-75); 
    %xtrap = cx.*matrix_truncated_column;
    vpd=kpd.*datk(:,chxpd); 
    lengthdata=size(datk,1);
    fsamp=1e6/(Samp_rate);
    fexc=freq_list(kf);
    
    % From Fisher 2007, data handling section, fexc must be an integer
    % multiple (K_integer) of the fwind, in other words, the exc_period must be an
    % integer multiple of the wind_period. Wind_period (the time for one window)
    % must be equal to Nfft multipled by the time for one sampling period.
    % Together, this means:
    %
    %           fexc = K_integer * fsamp / Nfft
    %                        i.e.
    %           Nfft = K_integer * fsamp / fexc
    % 
    % Now, in theory as long as K_integer >= 1, this will center Nfft. 
    % Therefore, for now, arbitrarily pick K_integer to get the best FRF results
    % while satisfying the above conditions. 
    
    checkperiod = fsamp/fexc; % equivalent to excitation period/sampling period
    K_integer_max = floor(lengthdata/checkperiod) % gives out the maximum K_integer for maximum Nfft size <= lengthdata  
    %K_integer = K_integer_pre - floor(0.5*K_integer_pre)
    %K_integer = K_integer*10
    %K_integer = 20;%20 for fexc(2) 268/0.95
    %Remember here that big Nfft = small bin size. Bin size too small will
    %cause response to spread to adjacent bins...Bin too large will give low
    % coherence because of averaging over frequencies to far away from
    % fexc. Need to manually find what K_integer is best by looking at FRF
    
    K_integer = K_integerOPT(kf)
    
    %if K_integer > k_integer_max
     %   return
    %end
   
  
    Nfft = K_integer*floor(fsamp/fexc) % both Nfft and Noverlap must be integers
    Noverlap = floor(Nfft/2); % typically, we use 50% overlap, i.e Nfft/2
    
    [P,F]=spectrum(xtrap,vpd,Nfft,Noverlap,hanning(Nfft),fsamp,0.95);
    %[P,F]=spectrum(xtrap,vpd,mk,Noverlap,hanning(Nfft),fsamp,0.99);
    
    % P has 5 colums:
    %       P(:,1) = input = Pxx = Power spectral density of input
    %       P(:,2) = output = Pyy = Power spectral density of output
    %       P(:,3) = ?????? = Spectra density of signal 1 and 2
    %       P(:,4) = TF estimate = Pxy 
    %       P(:,5) = coherence 
    if_range=find(abs(F-fexc)<0.001*fexc); % to find a small range of freq. around excitation freq.
    if_des=if_range(floor(median(find(P(if_range,5)==max(P(if_range,5)))))); % picking the index that gives the maximum coherence in the vicinity of the excitation freq.
    %if_des=if_range(floor(median(find(P(if_range,1)==max(P(if_range,1))))));
    %%picking the index that gives the max Pxx in the vicinity of the excitation freq.
    
    fexc= F(if_des);
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
PHfrf =abs(PHfrf - round(PHfrf./pi).*pi)*180/pi;% ???

PYm=mean(PYY,2);
Fi=FP(1:(end-1));
Pyi=cumtrapz(FP,PYm);
  
jkp=[];
for kkp=1:(numel(FRF_f)/2);
    jfrf=2*(kkp-1)+[1,2];
    jkkp=find(Fexc>=FRF_f(jfrf(1)) & Fexc<=FRF_f(jfrf(2))&Fexc>0);
    jkp=[jkp,jkkp];
end

Fexc=Fexc(jkp);
Mfrf=Mfrf(jkp);
PHfrf=PHfrf(jkp);
C=C(jkp);
H=H(jkp);
Pyy=Pyy(jkp);
Pyy2=Pyy2(jkp);


%%%%% plotting everything 
figure, subplot(211), semilogx(Fexc,Mfrf,'o','linewidth',2)
ylabel('Mag, dB'), set(gca,'Xlim',[0.01 5e3]);
subplot(212), semilogx(Fexc,PHfrf,'o','linewidth',2)
ylabel('Phase, deg'), xlabel('f, Hz'), set(gca,'Xlim',[0.01 5e3]);
hfrf=gcf;

figure, semilogx(Fexc,C,'o','linewidth',2)
ylabel('Coherence'), xlabel('f, Hz')

figure(hfrf), subplot(211),

figure, loglog(FP,PYm,'b','linewidth',2), hps=gcf;
xlabel('f (Hz)'), ylabel('Power Spectrum ($V^2*s$)')



cd(currentFolder)
savename=horzcat('name',num2str(chxpd));
save(savename)
   


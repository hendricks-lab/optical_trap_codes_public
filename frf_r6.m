function [mag,phase]=frf_r6(u,f,Rbead,Zbead,kT,falias)%,f0);

%u=[Minf, fp, Dfit, phi]
j=sqrt(-1);
Minf=1/u(1);
gammar=u(2);
alpha=u(3);
ktrap=u(4);
kcyt0=u(5);

kcyt1=u(6);
%keq0=u(6);
%w0=keq0/(gammar*9.42e-6);
%kcyt1=(u(6)-ktrap-kcyt0)/(w0^alpha);
% wm0=(gammar*9.42e-6)/m;
% kcyt1=(keq0-ktrap-kcyt0)/(wm0^alpha);

m=u(7)*1e-21;
nu=u(8)*1e12;

gamma_frf=gammar*9.42e-6.*ones(size(f));%gamma_r5(u,f,Rbead,Zbead,kT);

wfrf=2*pi.*f;

% kcyt=kcyt0+kcyt1.*wfrf.^alpha;
% keq=ktrap+kcyt;
% num_real=-m.*wfrf.^2-imag(gamma_frf).*wfrf + kcyt;
% num_imag=real(gamma_frf).*wfrf;
% den_real=keq-m.*wfrf.^2-imag(gamma_frf).*wfrf;
% den_imag=real(gamma_frf).*wfrf;
% Mth=Minf.*sqrt(num_real.^2+num_imag.^2)./...
%     sqrt(den_real.^2+den_imag.^2);
% PHth=atan2(num_imag,num_real)-...
%     atan2(den_imag,den_real);
% Mpd=(2*pi.*falias)./sqrt(wfrf.^2+(2*pi.*falias).^2);

Num=kcyt0-m.*wfrf.^2+gamma_frf.*wfrf.*j+(kcyt1.*(j.*wfrf).^(alpha))./gamma(alpha);
Den=ktrap+kcyt0-m.*wfrf.^2+gamma_frf.*wfrf.*j+(kcyt1.*(j.*wfrf).^(alpha))./gamma(alpha);

Mth=Minf.*sqrt(real(Num).^2+imag(Num).^2)./...
     sqrt(real(Den).^2+imag(Den).^2);
 PHth=atan2(imag(Num),real(Num))-...
     atan2(imag(Den),real(Den));
Mpd=(2*pi.*falias)./sqrt(wfrf.^2+(2*pi.*falias).^2);

PHpd=-atan2(wfrf,2*pi*falias);
mag=Mth.*Mpd;
phase=(180/pi).*(PHth+PHpd);

end
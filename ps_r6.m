function PS=ps_r6(u,f,fsamp,Rbead,Zbead,kT,falias)%,f0);

j=sqrt(-1);
Minf=1/u(1);
gammar=u(2);
alpha=u(3);
ktrap=u(4);
kcyt0=u(5);

kcyt1=u(6);
%keq0=u(6);
%w0=keq0/(gammar*9.42e-6);
%kcyt1=(keq0-ktrap-kcyt0)/(w0^alpha);
% wm0=(gammar*9.42e-6)/m;
% kcyt1=(keq0-ktrap-kcyt0)/(wm0^alpha);

m=u(7)*1e-21;
nu=u(8)*1e12;



gamma_fps=gammar*9.42e-6.*ones(size(f));%gamma_r5(u,f,Rbead,Zbead,kT);

% gamma0=kT/(Dfit*(1/Minf)^2);
% fm=gamma0/(2*pi*m);
% gamma_fps=gamma_r5(u,f,Rbead,Zbead,kT);

w=2*pi.*f;
%kcyt=kcyt0+kcyt1.*w.^alpha;

%gps_re=real(gamma_fps./gamma0);
%gps_im=imag(gamma_fps./gamma0);
% Pth_t1s=((Dfit.*gps_re)./((pi^2).*((feq+f.*gps_im-(f.^phi)./fm).^2 + (f.*gps_re).^2))).*...
%     (falias^2./(falias^2 + f.^2));

% Pth_t1s=(Minf.^2).*((4*kT.*gamma_fps)./(ktrap+kcyt0+kcyt1).^2) .* 1./(1+(gamma_fps.*w./(ktrap+kcyt)).^2).*...
%     (falias^2./(falias^2 + f.^2));

% finv=fsamp-f;
% winv=2*pi.*finv;
% Pth_t2s=(Minf.^2).*((4*kT.*gamma_fps)./(ktrap+kcyt0+kcyt1).^2) .* 1./(1+(gamma_fps.*winv./(ktrap+kcyt)).^2).*...
%     (falias^2./(falias^2 + finv.^2));
% 
% Pth_t1s=(Minf.^2).*((4*kT.*gamma_fps)./(ktrap+kcyt).^2) .* 1./(1+(gamma_fps.*w./(ktrap+kcyt)).^2).*...
%     (falias^2./(falias^2 + f.^2));

% Pth_t1s=(Minf.^2).*((4*kT.*gamma_fps)./(keq0).^2) .* 1./(1+(gamma_fps.*w./(ktrap+kcyt)).^2).*...
%     (falias^2./(falias^2 + f.^2));

% Pth_t1s=(Minf.^2).*((4*kT.*gamma_fps)./(keq0)) .* (ktrap+kcyt)./((ktrap+kcyt).^2+(gamma_fps.*w).^2).*...
%     (falias^2./(falias^2 + f.^2));


g=ktrap + kcyt0 + j.*w.*gamma_fps + (kcyt1.*(j.*w).^(alpha))./gamma(alpha); %g=6*pi*R*G

%Pth_t1s=(Minf.^2).*4.*kT.*gamma_fps./(real(g).^2+imag(g).^2);

Pth_t1s=(Minf.^2).*2.*kT.*imag(g)./(w.*((real(g).^2+imag(g).^2))).*...
     (falias^2./(falias^2 + f.^2)); %see Lau et al. (2003) PRL, eqn. 2

PS=real(Pth_t1s);%+real(Pth_t2s);
    
end
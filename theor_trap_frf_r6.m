function Mthr=theor_trap_frf_r6(u,fdata,ydata1,phdata1,W,ydata2,fsamp,Rbead,Zbead,kT,falias)%,f0);

i=sqrt(-1);

n1=numel(ydata1);
n2=numel(phdata1);
ffrf=fdata(1:n1);
fph=fdata([n1+1:n2]);
fps=fdata((n1+n2+1):end);
Wf=W*numel(fps)/(n1+n2);
%Wf=Wf.*exp(-ffrf./100);%Wf-ffrf./Wf;%
%Wf=Wf.*exp(-ffrf./50);

[Mth,PHth]=frf_r6(u,ffrf,Rbead,Zbead,kT,falias);
Mthr1=Wf.*(Mth-ydata1)./ydata1;

if isempty(phdata1)==1, %zph=0
    PHthr1=[];
else
   PHthr1=(Wf).*(PHth-phdata1)./phdata1; 
end
     
Pth_t=ps_r6(u,fps,fsamp,Rbead,Zbead,kT,falias);%,
IPth=(cumtrapz(fps,Pth_t)-ydata2)./ydata2;
IPth(1)=0;

Mthr=[Mthr1,PHthr1,IPth];
end
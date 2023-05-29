%���ؽ���ģ�͹�ʽ
%ʹ����˫�ؽ���ģ�͵����ݣ���֤��������Ƿ���ȷ
%����������ʲô�ط�����ȷ����Ҫ����
%���������ܵ���ѹ��ͼ����������Ҫ���ĵ���
%����ģ��������
%   ����ģ��Ϊ���Σ���Ϊ�������򣬴���������Ϊѹ������ֻ���Ǵ����ѷ죩��
%   ���м�����ѷ졢���ʡ��ܶ����򣨲��д���ģ�ͣ��ѷ�Ϊ�����ռ䣩��
%   ������ǻ��ʺ��ܶ����ڵ�˫�ؽ��������ܶ�������ʣ�����Ϊ�����տռ䣩
clear;clc;close all
% p0=13.9*1e6;              pw=1.8*1e6; 
p0=13.9*1e6;              pw=1.8*1e6;            mu=0.0008;  
B=1.19;                 q=10/24/3600;
N=8;                      num=1000;               S=2;%2
thetaf=0.05;                     df=2;
%����
xf=100;            h=15.44;%=z
wf=0.0746;          wfD=wf./xf;
yf=100;            yfD=yf./xf;
bf=0.0746;          bfD=bf./xf;
ye=200;            yeD=ye./xf;
%��϶��
phi3v0=0.08;               phi3m0=0.08;
phi2f0=0.0045;               phi2m0=0.08;               
phi2v0=0.08;
phi10=0.0045;
%��͸��
k3m=0.23./1e15;                k3v=0.23./1e15;%(ָ�����ڿ׶�����Χ�ѷ�֮���С�ѷ����͸��)
k2f0=500./1e15;             k2m0=0.23./1e15;          
k2v0=0.23./1e15;
k1=10000./1e15;
%��������
lamda2v=0.1;
omega2v=0.45;
lamda2m=0.1;
omega2m=0.45;
omega2f=0.1;
lamda3v=0.1;
omega3v=0.45;
%ѹ��ϵ��
Ct3m=0.000466./1e6;             Ct3v=0.000466./1e6;          
Cl3=0.000466./1e6;             Cl2f=0.000466./1e6./2.0;
Ct2v=0.000466./1e6;             Ct2m=0.000466./1e6;    
Ct2f=0.000466./1e6;
Ct1=0.000466./1e6;

G=0.011.*1e6;                         GD=G.*Cl2f.*xf;
gD=G.*k2f0.*h.*xf./(q.*mu.*B);                             
phiCt2=phi2m0.*Ct2m+phi2f0.*Ct2f+phi2v0.*Ct2v;
% C=10/1e6;     CD=C./(h.*xf.*xf.*phiCt2);
CD=0.08;

tini(:,1)=logspace(2,9,num);%��λ��s
T(:,1)=tini.*k2f0/(phiCt2*mu*xf^2); %����λ�ʱ��
i=1:1:N;
s=log(2)./T*i;

%��ʽ�����˸���
for i=1:N
   k=fix((i+1)./2.0):min(i,N./2.0);
   s2=sum(k.^(N./2.0).*factorial(2.*k)./...
       (factorial(N./2.0-k).*factorial(k).*...
       factorial(k-1).*factorial(i-k).*factorial(2.*k-i)));
   vs(i,1)=(-1).^(N./2.0+i).*s2;
end

%���ʺͷ춴˫�ؽ���ģ��
a31=GD.*Cl3./Cl2f;
omega3m=phi3m0.*Ct3m./phiCt2;
b31=lamda3v.*k2f0./k3m.*omega3v.*s./(omega3v.*s+lamda3v)+...
    omega3m.*k2f0.*s./k3m;
r31=(a31+(a31.^2+4.*b31).^0.5)./2.0;
r32=(a31-(a31.^2+4.*b31).^0.5)./2.0;
a32=r31./r32.*exp((r31-r32).*yeD);
F32=(r31.*exp(r31.*yfD)-a32.*r32.*exp(r32.*yfD))./...
    (exp(r31.*yfD)-a32.*exp(r32.*yfD));

%���ʣ��춴���ѷ����ؽ�������
% lamda2v=alphav2.*k2v0.*xf.*xf./k2f0;
% omega2v=phi2v0.*Ct2v./phiCt2;
% lamda2m=alpham2.*k2m0.*xf.*xf./k2f0;
% omega2m=phi2m0.*Ct2m./phiCt2;
% omega2f=phi2f0.*Ct2f./phiCt2;
b21=(lamda2v.*omega2v.*s./(omega2v.*s+lamda2v)+...
    lamda2m.*omega2m.*s./(omega2m.*s+lamda2m)+...
    omega2f.*s).*(1.0./bfD).^thetaf;
m21=k3m./k2f0.*(yfD./bfD).^(thetaf+2-df);
n21=b21.^0.5.*yfD.^(thetaf./2.0).*...
    besseli((1-df)./(2+thetaf),2.*b21.^0.5./(2+thetaf).*yfD.^(1+thetaf./2.0));
n22=b21.^0.5.*yfD.^(thetaf./2.0).*...
    besselk((1-df)./(2+thetaf),2.*b21.^0.5./(2+thetaf).*yfD.^(1+thetaf./2.0));
n23=m21.*F32.*...
    besseli((3+thetaf-df)./(2+thetaf),2.*b21.^0.5./(2+thetaf).*yfD.^(1+thetaf./2.0));
n24=m21.*F32.*...
    besselk((3+thetaf-df)./(2+thetaf),2.*b21.^0.5./(2+thetaf).*yfD.^(1+thetaf./2.0));
n25=wfD.^((3+thetaf-df)./2.0).*...
    (besseli((3+thetaf-df)./(2+thetaf),2.*b21.^0.5./(2+thetaf).*wfD.^(1+thetaf./2.0))+...
    (n21-n23)./(n22+n24).*...
    besselk((3+thetaf-df)./(2+thetaf),2.*b21.^0.5./(2+thetaf).*wfD.^(1+thetaf./2.0)));
a21=b21.^0.5.*wfD.^((3+2.*thetaf-df)./2.0).*...
    besseli((1-df)./(2+thetaf),2.*b21.^0.5./(2+thetaf).*wfD.^(1+thetaf./2.0));
a22=b21.^0.5.*wfD.^((3+2.*thetaf-df)./2.0).*...
    besselk((1-df)./(2+thetaf),2.*b21.^0.5./(2+thetaf).*wfD.^(1+thetaf./2.0));
F21=a21./n25+(n23-n21)./(n22+n24).*a22./n25;

%ѹ������HF��
a11=phi10.*Ct1.*k2f0.*s./phiCt2./k1-k2f0./k1./wfD.*(wfD./bfD).^(df-thetaf-2).*F21;
r11=a11.^0.5;
r12=-a11.^0.5;
C0=-k2f0./2./k1./wfD;
C11=C0./s./(r11-r11.*exp(r11-r12));
C12=-C11.*r11./r12.*exp(r11-r12);

xD=0.0;

p1Ds=C11.*exp(r11.*xD)+C12.*exp(r12.*xD);
%���ؽ���ģ��
 
pwDs=(s.*p1Ds+S)./(s.*(1+s.*CD.*(s.*p1Ds+S)));

pwD=log(2)./T.*(pwDs*vs);

% ��ѹ�����£���ʽ�ռ����;�����
qDs=1./(s.^2.*pwDs);%

% ��ʵ�ռ��£��;�����
qD=log(2)./T.*(qDs*vs);%

qtotal=qD*k2f0*h*(p0-pw)/mu/B*24*3600;     %   ��λ��m3/Day
%����������ͨ���ȷ��ݵõ�ʱ������β���������λ�������ͨ��������ٷ��ݣ���ֵӦ����Ψһ�İɣ�����
dt=T-circshift(T,[1,0]);
dt(1,1)=T(1,1);
Qtotal=zeros(num,1);   
for ii=2:num
  if isnan(qtotal(ii,1))
    Qtotal(ii)=Qtotal(ii-1)+0;
  elseif qtotal(ii-1)>=0
      Qtotal(ii)=Qtotal(ii-1)+(qtotal(ii-1)+qtotal(ii))/2*dt(ii); 
  else 
     Qtotal(ii)=Qtotal(ii-1)+0;
  end
end
 
t=tini/24/3600;%��λ����

figure(1)
%  col={'g','r','m','k','b','c','r','b','g','k','m','y','w'};%������ɫ�ľ���
%     xiankuan={1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5};
loglog(T,pwD)%ѹ����������
hold on
[mr mc]=size(pwD);
for i=1:1:mr-1
    dpwD(i,1)=(pwD(i+1,1)-pwD(i,1))./(log(T(i+1,1))-log(T(i,1)));%%Ϊʲô�����õ���log����-log����
    ddt(i,1)=T(i,1);
end
    loglog(ddt,dpwD,'--');%ѹ���������ߵĻ���
grid on
figure(2)
plot(t,qD)%�ղ�����������
hold on
[mr mc]=size(qD);
for i=1:1:mr-1
    dq(i,1)=-(qD(i+1,1)-qD(i,1))./(log(t(i+1,1))-log(t(i,1)));%%Ϊʲô�����õ���log����-log����
    ddt(i,1)=t(i,1);
end
%loglog(ddt,dq,col{kol},'linestyle','--','linewidth',1);%�ղ����������ߵĻ���
    grid on

% figure(3)
% [out_t,out_Qtotal]=zdrawQ(Qtotal,t,3,3);
% hold on
% % print(gcf,'-dpng','-r1200',...
% %      'C:\Users\HP\Desktop\zuhui./ģ���ܲ�ͼ(����)')
% 
% % print(gcf,'-dpng','-r1200',...
% %      'C:\Users\HP\Desktop\huedetu./�����ѹ��ͼ')

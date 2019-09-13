%%%%%%%%%%%% The basic code can achieve SSA, which could be used for reducing noise,
%% and principal component extraction. In addition, the time series prediction
%% can be achieved on the basis of the code. The code is very useful.
%% The author is Wanqiu Li, 2019,09,13

%%%%%%%%%%%%%%%%%%%%%%%�����ı�����.txt%%%%%%%%%%%%%%%%%%%%%%%
%filename1=['GSLXTcc.txt'];
%fid1=fopen(filename1,'r');
%[gslx,COUNT1]=fscanf(fid1,'%f',inf);
%fprintf(fid1,'%f',gslx);
%fclose(fid1);
%height(1:2152)=gslx(6457:8608);

height=(1:2152)';
aa=load('GSLXTcc.txt');
yya=aa(:);
height=yya(6457:8608);  %%%  cm unit
x1=height';
L=1076;
%%%%%%%%%%���е���CORSվ���������%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % pi=3.1415926;
% % % % % % % % % % % % % t=(1:2000);
% % % % % % % % % % % % % H0=10*sin((2*pi/365)*t);
% % % % % % % % % % % % % %p=plot(t,H0,'r');
% % % % % % % % % % % % % %%%%%%������˹������%%%%%%%%%%%%
% % % % % % % % % % % % % noise=mvnrnd(0,1,2000)';
% % % % % % % % % % % % % H=H0+noise;
% % % % % % % % % % % % % %pp=plot(t,H,t,noise);
% % % % % % % % % % % % % %%%%%%%%%%%%�Բ�����ʱ�����н���SSA%%%%%%%%%%%%%
% % % % % % % % % % % % % %function [d,r,vr]=ssa(x1,L)
% % % % % % % % % % % % % x1=H;
% % % % % % % % % % % % % L=1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step1 : Build trayectory matrix

   N=length(x1); 
   if L>N/2;L=N-L;end
	K=N-L+1; 
   X=zeros(L,K);  
	for i=1:K
	  X(1:L,i)=x1(i:L+i-1); 
	end
    
% Step 2: SVD����ֵ�ֽ�

   S=X*X'; %
	[U,autoval]=eig(S);%eig���ؾ��������ֵ������������U������������autoval������ֵ
	[d,i]=sort(-diag(autoval));
    %x= diag(v,k)������v��Ԫ����Ϊ����X�ĵ�k���Խ���Ԫ�أ���k=0ʱ��vΪX�����Խ��ߣ���k>0ʱ��vΪ�Ϸ���k���Խ��ߣ���k<0ʱ��vΪ�·���k���Խ��ߡ�
    %sort(A��dim)��ʾ��A�еĸ���Ԫ���������У�sort(A)Ĭ�ϵ�������d������õ�������i������d��ÿһ���Ӧ��A���������
   d=-d;
   U=U(:,i);
   sev=sum(d); %����ֵ�ĺ�
	plot((d./sev)*100),hold on,plot((d./sev)*100,'rx');
	title('Singular Spectrum');xlabel('Eigenvalue Number');ylabel('Eigenvalue (% Norm of trajectory matrix retained)')  %�켣����ķ���
   V=(X')*U; 
   rc=U*V';
    %����ֵ�ֽ�S=U*D*V
% Step 3: Grouping

  % I=input('Choose the agrupation of components to reconstruct the series in the form I=[i1,i2:ik,...,iL] ')
   I=(20:20);
% % % % % % %    I=(1:40);
   Vt=V';
   rca=U(:,I)*Vt(I,:);

% Step 4: Reconstruction

   y=zeros(N,1);  
   Lp=min(L,K);
   Kp=max(L,K);

   for k=0:Lp-2
     for m=1:k+1;
      y(k+1)=y(k+1)+(1/(k+1))*rca(m,k-m+2);
     end
   end

   for k=Lp-1:Kp-1
     for m=1:Lp;
      y(k+1)=y(k+1)+(1/(Lp))*rca(m,k-m+2);
     end
   end

   for k=Kp:N
      for m=k-Kp+2:N-Kp+1;
       y(k+1)=y(k+1)+(1/(N-k))*rca(m,k-m+2);
     end
   end

   figure;subplot(2,1,1);hold on;xlabel('Data poit');ylabel('Original and reconstructed series')
   plot(x1,'b');grid on;plot(y,'r')
% plot(t,y)
   r=x1'-y;
   subplot(2,1,2);plot(r,'g');xlabel('Data poit');ylabel('Residual series');grid on
   vr=(sum(d(I))/sev)*100;
%%%%��RC�ɷֱַ��������
   %plot(y,'b');
%%%%%%%%%%%%%%������������,����ͼ%%%%%%%%%%%%%%%
%   xy1=y; 
%   xy2=y; 
%   xy3=y; 
%   xy4=y; 
%   xy5=y; 
%   xy6=y; 
%   xy7=y;
%   xy8=y;
%   xy9=y;
%   xy10=y;
%   xy11=y;
%   xy12=y;
%   xy13=y;
%   xy14=y;
%   xy15=y;
%   xy16=y;
%   xy17=y;
%   xy18=y;
%   xy19=y;
   xy20=y;

rc=[xy1,xy2,xy3,xy4,xy5,xy6,xy7,xy8,xy9,xy10,xy11,xy12,xy13,xy14,xy15,xy16,xy17,xy18,xy19,xy20];   
coef=zeros(20,20);
rr=zeros(20,20);
xishu=zeros(2,2);
for jj=1:20
     for kk=1:20
         if jj==kk
            rr(jj,kk)=1;
         else
            xishu=corrcoef(rc(:,jj),rc(:,kk));  
            rr(jj,kk)=xishu(1,2);
         end      
     end   
end
  
[aa,bb]=meshgrid(0.5:1:19.5,0.5:1:19.5); 
zz=rr; 
surf(aa,bb,zz)    %%% meshɫ���ǲ������ģ�surfɫ���Ǿ�����ֵ��
view(2)
%contourf(aa,bb,zz)
colorbar
title('���ϵ��')
%%%%%%%%%%%%%%%%%�ϲ��ɷ֣���ͼ��Ƶ�׷���%%%%%%%%%%%%%%%%%%%%%%%%%
%��һ�л�ͼ
figure;subplot(8,2,1);hold on;ylabel('RC1');plot(xy1,'b');title('(a) ��վʱ��RC�ϲ�');
n=length(xy1');
t=1:n;
N=2^12;
f=linspace(0,1,N+1);
yy=[xy1',zeros(1,N-n)];
Y=fft(yy);
P=Y.*conj(Y)/N;
D=P(1:N/2+1);
subplot(8,2,2);plot(f(1:N/2+1),D,'r');ylabel('���ܶ�');title('(b) RC���ɷ�Ƶ��ͼ');
%�ڶ��л�ͼ
subplot(8,2,3);plot(xy2+xy3,'b');ylabel('RC2+RC3');
n=length((xy2+xy3)');
t=1:n;
N=2^12;
f=linspace(0,1,N+1);
yy=[(xy2+xy3)',zeros(1,N-n)];
Y=fft(yy);
P=Y.*conj(Y)/N;
D=P(1:N/2+1);
subplot(8,2,4);plot(f(1:N/2+1),D,'r');
%�����л�ͼ
subplot(8,2,5);plot(xy4+xy5,'b');ylabel('RC4+RC5');  
n=length((xy4+xy5)');
t=1:n;
N=2^12;
f=linspace(0,1,N+1);
yy=[(xy4+xy5)',zeros(1,N-n)];
Y=fft(yy);
P=Y.*conj(Y)/N;
D=P(1:N/2+1);
subplot(8,2,6);plot(f(1:N/2+1),D,'r');
%�����л�ͼ   
subplot(8,2,7);plot(xy6+xy7,'b');ylabel('RC6+RC7');
n=length((xy6+xy7)');
t=1:n;
N=2^12;
f=linspace(0,1,N+1);
yy=[(xy6+xy7)',zeros(1,N-n)];
Y=fft(yy);
P=Y.*conj(Y)/N;
D=P(1:N/2+1);
subplot(8,2,8);plot(f(1:N/2+1),D,'r');
%�����л�ͼ
subplot(8,2,9);plot(xy8+xy9,'b');ylabel('RC8+RC9');
n=length((xy8+xy9)');
t=1:n;
N=2^12;
f=linspace(0,1,N+1);
yy=[(xy8+xy9)',zeros(1,N-n)];
Y=fft(yy);
P=Y.*conj(Y)/N;
D=P(1:N/2+1);
subplot(8,2,10);plot(f(1:N/2+1),D,'r');
%�����л�ͼ
subplot(8,2,11);plot(xy10+xy11,'b');ylabel('RC10+RC11');
n=length((xy10+xy11)');
t=1:n;
N=2^12;
f=linspace(0,1,N+1);
yy=[(xy10+xy11)',zeros(1,N-n)];
Y=fft(yy);
P=Y.*conj(Y)/N;
D=P(1:N/2+1);
subplot(8,2,12);plot(f(1:N/2+1),D,'r');
%�����л�ͼ
subplot(8,2,13);plot(xy12+xy13,'b');ylabel('RC12+RC13');
n=length((xy12+xy13)');
t=1:n;
N=2^12;
f=linspace(0,1,N+1);
yy=[(xy12+xy13)',zeros(1,N-n)];
Y=fft(yy);
P=Y.*conj(Y)/N;
D=P(1:N/2+1);
subplot(8,2,14);plot(f(1:N/2+1),D,'r');
%�ڰ��л�ͼ
subplot(8,2,15);plot(xy14+xy15,'b');xlabel('ʱ��/��');ylabel('RC14+RC15');
n=length((xy14+xy15)');
t=1:n;
N=2^12;
f=linspace(0,1,N+1);
yy=[(xy14+xy15)',zeros(1,N-n)];
Y=fft(yy);
P=Y.*conj(Y)/N;
D=P(1:N/2+1);
subplot(8,2,16);plot(f(1:N/2+1),D,'r');xlabel('Ƶ��/Hz');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%ͳ����
   %cgmean=mean(y)
   chonggou=std(y)
   cancha=std(r)
   yuanshi=std(x1)
   junzhi=mean(r)
  % bias=std(y-x1')
%���ƹ��������ߡ���Ƶ�ʣ��ǵ�λʱ������������Ա仯�Ĵ�����
%�����������˶�Ƶ���̶ȵ��������÷���f��ͱ�ʾ����λΪ���֮һ������Ϊs
   fs=1/365;
   gg=fft(r,N);
   mag=abs(gg); 
   f=t*fs/N;
   plot(f,mag);
% junzhi1=mean(H0)
% junzhi2=mean(y)
  
%%��������ֵ�ֲ�
figure;
d=d'
plot(d,'ro')

%%%%������״ͼ
%bias=std(y-x1)
yyy=[63.1,41.4,0.76,0.61,0.75,0.61,0.55,0.75];
bar(yyy)

%%���Ƹ�˹�ֲ�ֱ��ͼ
j=0;
for i=1:2000
    if((r(i)<3.5)&&(r(i)>3.0)) 
       j=j+1; 
    end 
end

fanwei=[-3.25,-2.75,-2.25,-1.75,-1.25,-0.75,-0.25,0.25,0.75,1.25,1.75,2.25,2.75,3.25];
count=[4,12,30,96,210,311,356,388,284,172,81,41,9,5];
bar(fanwei,count);
%hold on,plot(fanwei,count,'r');




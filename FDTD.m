%sweeping freq at (0.3,0.3)
clear all
close all
%% single oscillator model
thickness=0.01;
width=1.5;
length=1;
density=392;%density of Sitka Spruce;
volume=thickness*width*length;
m=volume*density;%mass
E_long=12*10^9; %looked up on lecture slides, youngs modulus in Pa
%% ISOTROPIC BOARD
%definitions of variables:
poisson_variable=0.3; %given in the paper
D=thickness^3*E_long/(12*(1-poisson_variable^2));
p=density;
time_diff = 5e-06;%in sec
count=1
length=1;
width=1.4;
x_diff = 0.02;%2cm
y_diff = 0.02;%2cm
x_int=1/x_diff;

force_length=0.3*x_int; %at 0.3m
force_height=0.3*x_int;
R = 1; %given in assignment
B = R*time_diff/(2*thickness*density);
a1 = 2 / (1+B);
a2 = (-1+B) / (1+B);
a3 = -D*(time_diff)^2/(p*thickness*(x_diff)^4);
a4 = -D*(time_diff)^2/(p*thickness*x_diff^4);
a5 = -2*D*(time_diff^2)/(p*thickness*x_diff^4);
a6 = time_diff^2/(p*thickness*x_diff^2);
z(:,:,1)=0; %initiliaze
i=3;
j=3;
n=3;
ll=10:10:100;
lll=100:100:1000;
llll=1000:1000:10000;
freq_array=[ll,lll,llll];
for f=freq_array %example freq
%1/f*100
time_total_sec=0.1;%we should have at least 50 periods1.5;
z=zeros(length/x_diff,width/y_diff,round(time_total_sec*1/time_diff)); %in (cm,cm,ms)
F=zeros(length/x_diff,width/y_diff,round(time_total_sec*1/time_diff));
for t=time_diff:time_diff:time_total_sec
F(force_length,force_height,round(t*1/time_diff))=1*sin(2*pi*f*t);
end
t=(time_diff:time_diff:time_total_sec)';
temp=squeeze(F(force_length,force_height,:));
% plot(t,temp)
% title('Force at (0.3m,0.3m)')
% xlabel('Time (s)')
% ylabel('Magnitude')
% set(gca,'FontSize',16)
%%
for n=3:time_total_sec*1/time_diff-3
for j=3:width*x_int-3
for i=3:length*x_int-3
z(i,j,n+1) = a1*z(i,j,n) + a2*z(i,j,n-1) + a3*(z(i+2,j,n)-4*z(i+1,j,n)+6*z(i,j,n)-4*z(i-1,j,n)+z(i-2,j,n))...
+a4*(z(i,j+2,n)-4*z(i,j+1,n)+6*z(i,j,n)-4*z(i,j-1,n)+z(i,j-2,n))+a5*(z(i+1,j+1,n)...
+z(i+1,j-1,n)+z(i-1,j+1,n)+z(i-1,j-1,n)-2*z(i+1,j,n)-2*z(i-1,j,n)-2*z(i,j+1,n)...
-2*z(i,j-1,n)+4*z(i,j,n))+a6*F(i,j,n);
end
end
end
%%
% figure
% h=surf(squeeze(z(:,:,500)))
% set(h,'LineStyle','none')
% title('Displacement graph of sound board, f=500Hz')
% set(gca,'FontSize',16)
% xlabel('Length')
% ylabel('Height')
%%
% figure
temp=squeeze(z(force_length,force_height,:));
%z(force_length,force_height,end-3)
% plot(t,temp)
% title('Displacement at the force location')
% xlabel('Samples')

%ylim([-20000 20000])
%%
% hold on
t=(time_diff:time_diff:time_total_sec)';
dx = temp(2:end) - temp(1:end-1);
v = dx./time_diff;
% plot(v)
F_array=squeeze(F(force_length,force_height,:));
F_array=F_array(1:end-1);
% plot(t(1:end-1),abs(F_array)./abs(v))
% %ylim([0,2*10^5])
% %xlim([0.7,3])
% title('Impedance for at the force location, f=5Hz')
% xlabel('Time (s)')
% ylabel('|Z| (kg/s)')
%
% set(gca,'FontSize',16)
%%
% figure
% subplot(3,1,1)
% plot(t(1:end-1),v,'linewidth',2)
% title('Velocity(up) and Force(down) Amplitude Comparison')
% set(gca,'FontSize',15)
% ylabel('Amplitude')
% xlabel('Time (s)')
%xlim([1,996])
% subplot(3,1,2)
% plot(t(1:end-1),F_array,'r','linewidth',2)
% ylabel('Amplitude')
% xlabel('Time (s)')
% set(gca,'FontSize',16)
% subplot(3,1,3)
v_sample = v((end/8):(end/2)); %picked random stuff from the middle
% plot(v_sample)
% ylabel('Amplitude')
% xlabel('Samples')
% set(gca,'FontSize',16)
Z_point(count)=1/(abs((max(v_sample)-min(v_sample)))/2); %take the peak to peak amplitude
% a better way would be to filter the response, average the peaks and find
%a more reliable value but yeah.
count=count+1
end
%%
close all
figure
plot(freq_array,Z_point,'r','linewidth',2)
xlabel('Frequency in Hz')
ylabel('|Z| (kg/s)')
title('Driving point impedance for different frequencies')
%ylim([0 10000])
set(gca,'FontSize',18)
hold on
stem(freq_array,Z_point)
figure
loglog(freq_array,Z_point,'--om','linewidth',2)
hold on
loglog(freq_array,Z_point,'m','linewidth',2)
set(gca,'FontSize',18)
xlabel('Frequency in Hz')
ylabel('|Z| (kg/s)')

title('Log graph, Driving point impedance for different frequencies')
ylim([10,10000])
 

%PART B)
%sweepind location fixed freq
clear all
close all
%% single oscillator model
thickness=0.01;
density=392;%density of Sitka Spruce;
length=1;
width=1.4;
volume=thickness*width*length;
m=volume*density;%mass
E_long=12*10^9; %looked up on lecture slides, youngs modulus in Pa
%% ISOTROPIC BOARD
%definitions of variables:
poisson_variable=0.3; %given in the paper
D=thickness^3*E_long/(12*(1-poisson_variable^2));
p=density;
time_diff = 5e-06;%in sec
count=1
length=1;
width=1.4;
x_diff = 0.02;%2cm
y_diff = 0.02;%2cm
x_int=1/x_diff;
y_int=1/y_diff;
R = 1; %given in assignment
B = R*time_diff/(2*thickness*density);
a1 = 2 / (1+B);
a2 = (-1+B) / (1+B);
a3 = -D*(time_diff)^2/(p*thickness*(x_diff)^4);
a4 = -D*(time_diff)^2/(p*thickness*x_diff^4);
a5 = -2*D*(time_diff^2)/(p*thickness*x_diff^4);
a6 = time_diff^2/(p*thickness*x_diff^2);
z(:,:,1)=0; %initiliaze
i=3;
j=3;
n=3;
ll=10:10:100;
lll=100:100:1000;
llll=1000:1000:10000;
freq_array=[ll,lll,llll];
f=349.23;%example freq
%1/f*100
time_total_sec=1/f*6;%we should have at least 50 periods1.5;
z=zeros(length/x_diff,width/y_diff,round(time_total_sec*1/time_diff)); %in (cm,cm,ms)
zero_array=zeros(length/x_diff,width/y_diff,round(time_total_sec*1/time_diff));
F=zero_array;
count_length=1;
count_width=1;

for force_length=(x_diff*4)*x_int : x_diff*x_int : length*x_int-(x_diff*4)*x_int %at 0.3m
count_width=1;
for force_width=(y_diff*4)*y_int : y_diff*y_int : width*y_int-(y_diff*4)*y_int;
F=zero_array;
for t=time_diff:time_diff:time_total_sec
F(force_length,force_width,round(t*1/time_diff))=1*sin(2*pi*f*t);
end
if 1==0
t=(time_diff:time_diff:time_total_sec)';
temp=squeeze(F(force_length,force_width,:));
plot(t,temp)
title('Force at (0.3m,0.3m)')
xlabel('Time (s)')
ylabel('Magnitude')
set(gca,'FontSize',16)
end
%%
for n=3:time_total_sec*1/time_diff-3
for j=3:width*x_int-3
for i=3:length*x_int-3
z(i,j,n+1) = a1*z(i,j,n) + a2*z(i,j,n-1) + a3*(z(i+2,j,n)-4*z(i+1,j,n)+6*z(i,j,n)-4*z(i-1,j,n)+z(i-2,j,n))...
+a4*(z(i,j+2,n)-4*z(i,j+1,n)+6*z(i,j,n)-4*z(i,j-1,n)+z(i,j-2,n))+a5*(z(i+1,j+1,n)...
+z(i+1,j-1,n)+z(i-1,j+1,n)+z(i-1,j-1,n)-2*z(i+1,j,n)-2*z(i-1,j,n)-2*z(i,j+1,n)...
-2*z(i,j-1,n)+4*z(i,j,n))+a6*F(i,j,n);
end
end
end
%%
% figure
% h=surf(squeeze(z(:,:,500)))
% set(h,'LineStyle','none')
% title('Displacement graph of sound board, f=500Hz')
% set(gca,'FontSize',16)
% xlabel('Length')
% ylabel('Height')
%%
% figure
temp=squeeze(z(force_length,force_width,:));
%z(force_length,force_height,end-3)
% plot(t,temp)
% title('Displacement at the force location')
% xlabel('Samples')
%ylim([-20000 20000])
%%
% hold on
t=(time_diff:time_diff:time_total_sec)';
dx = temp(2:end) - temp(1:end-1);
v = dx./time_diff;
% plot(v)
F_array=squeeze(F(force_length,force_width,:));
F_array=F_array(1:end-1);
% plot(t(1:end-1),abs(F_array)./abs(v))
% %ylim([0,2*10^5])
% %xlim([0.7,3])
% title('Impedance for at the force location, f=5Hz')
% xlabel('Time (s)')
% ylabel('|Z| (kg/s)')
%
% set(gca,'FontSize',16)
%%
% figure

% subplot(3,1,1)
% plot(t(1:end-1),v,'linewidth',2)
% title('Velocity(up) and Force(down) Amplitude Comparison')
% set(gca,'FontSize',15)
% ylabel('Amplitude')
% xlabel('Time (s)')
%xlim([1,996])
% subplot(3,1,2)
% plot(t(1:end-1),F_array,'r','linewidth',2)
% ylabel('Amplitude')
% xlabel('Time (s)')
% set(gca,'FontSize',16)
% subplot(3,1,3)
v = v(1:(end-100)); %picked random stuff from the middle
% plot(v_sample)
% ylabel('Amplitude')
% xlabel('Samples')
% set(gca,'FontSize',16)
Z_point(count_length,count_width)=1/(abs((max(v)-min(v)))/2); %take the peak to peak amplitude
% a better way would be to filter the response, average the peaks and find
%a more reliable value but yeah.
count_width=count_width+1
end
count_length=count_length+1
end
%%
%plot(Z_point(:,5))
%%
figure
h=surf(Z_point)
set(h,'LineStyle','none')
title('Impedance graph of soundboard for different force points, f=349.23Hz')
set(gca,'FontSize',16)
xlabel('Length')
ylabel('Width')
zlabel('|Z| (kg/s)')
zlim([0 500])
colorbar
lim = caxis
caxis([0 500])
%% Part C and D
%% part c
epsilon=0.01;
freq=[349.23,440,523.25,659.25,783.99];
Z_freq=[105.2,136.4,253.902,154.9,172.4];
wo=2*pi*Z_freq;
Yb=1./Z_freq;
pS=0.0108;%kg/m^3
F0=1;
time_step=0.01;
t=0:time_step:2;
T=1200;%from the slides
for i=1:numel(Yb)
Z0(i)=sqrt(pS*T);
beta=Yb(i)/pi*1i*Z0(i);%formulas from the book

mu=(epsilon^2+beta^2)^(1/2);
a(1,i)=beta+epsilon+mu;%ora(2,i)=beta+epsilon-mu;%orVb =
2.*pi.*F0.*beta/mu/Z0(i).*exp(1j.*(epsilon+beta).*wo(i).*t).*(mu.*cos(mu.*wo(i).*t)+1j.*epsilon.*sin(mu.*wo(i).*t)).*e
xp(1j.*wo(i).*t);
plot(t,mag2db(abs(Vb)),'LineWidth',2)
decay_time=find(abs(Vb) < 0.002);
decay_time_array(i)=decay_time(1)*time_step;
hold on
end
real(a)
a
legend('349.23','440','523.25','659.25','783.99')
title('Different Strings Bridge Velocity vs time')
xlabel('Time in s')
ylabel('Magnitude in db')
set(gca,'FontSize',16)
%a^2 - 2*(beta+epsilon)*a + 2*epsilon*beta = 0;
%beta=(a^2+2*epsilon*a)/(-2*a+2*epsilon);
%% part 
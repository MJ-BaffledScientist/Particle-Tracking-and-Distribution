%clear all;

%info for old phase 4, young phase 4 and 5
%frame interval of 32ms and 9.75pixels = 1micrometer
%info for morning pALL
%frame interval of 30ms and 453pixels = 70micrometer
%info for afternoon pALL and original
%frame interval of 30ms and 410pixels=40micrometer

load('HeLaM_Lysobrite_CombinedExp1-3_Kif5c')
%load("Complete_longphase.mat");
%set the scale of distance per pixel
%lengthScale = 40e-6/410;
lengthScale = 70e-6/270;
lengthScale = lengthScale*1e6;
%set time scale of seconds per frame
timeScale = 30e-3;
% timeScale = 98.5e-3;
framerate = 1/timeScale;
%parameter for boundary in scale
boundaries = [0.0:1:15];

[s, fpp, count] = calcFPP_parallel(tr,framerate,lengthScale,boundaries,1);

saveas(gcf,['HeLaM_Lysobrite_FPP_Kif5c'],'epsc')
saveas(gcf,['HeLaM_Lysobrite_FPP_Kif5c'],'png')
savefig([['HeLaM_Lysobrite_FPP_Kif5c']])
    
% %create function that plots average velocity with length travelled
% [Laglen,Lagvel,avglen,avgvel,len, vel, time] = Velocity_Length(tr,timeScale,lengthScale);
% 
% figure
% subplot(3,2,1)
% plot(avglen,avgvel,'.')
% xlabel('Average Length')
% ylabel('Average Velocity')
% 
% %figure
% subplot(3,2,2)
% plot(time,avglen,'.')
% xlabel('Time')
% ylabel('Average Length')
% 
% %figure
% subplot(3,2,3)
% plot(time,avgvel,'.')
% xlabel('Time')
% ylabel('Average Velocity')
% 
% %figure
% subplot(3,2,4)
% plot(len,vel,'.')
% xlabel('Length')
% ylabel('Velocity')
% 
% %figure
% subplot(3,2,5)
% plot(len,time,'.')
% xlabel('Length')
% ylabel('Time')
% 
% %figure
% subplot(3,2,6)
% plot(time,vel,'.')
% xlabel('Time')
% ylabel('Velocity')
% 
% figure
% hold on
% for i = 1:length(Laglen)
%     plot(Lagvel{i},Laglen{i},'-')
% end
% ylabel('Length')
% xlabel('Velocity')
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% 
% %probability vs length travelled




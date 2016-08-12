% Script to load up ASP tile parameters

% Load ASP paramaters
%par = LoadASPpar();
% weights
wlow= 1;
wmid= 1;%0.2;%-0.2;
whigh= 1;%0.1;%-0.1;
wDC= 1; %0.005;%0.05;



w = [wlow  wlow  wlow  wlow  wlow  wlow  wlow  wlow  wlow  wlow  wlow  wlow  wlow  wlow  wlow  wlow  ... 
     wmid  wmid  wmid  wmid  wmid  wmid  wmid  wmid  wmid  wmid  wmid  wmid  wmid  wmid  wmid  wmid  ...
     whigh whigh whigh whigh whigh whigh whigh whigh whigh whigh whigh whigh whigh whigh whigh whigh ...
     wDC];

b_low=12;
b_mid=18;
b_high=24;

par(1,:)=  [b_low,0,180]; 
par(2,:)=  [b_low,0,0];
par(3,:)=  [b_low,0,270];  
par(4,:)=  [b_low,0,90];

par(5,:)=  [b_low,90,270];  
par(6,:)=  [b_low,90,90];  
par(7,:)=  [b_low,90,180];
par(8,:)=  [b_low,90,0];  

par(9,:)=  [b_low,45,270];
par(10,:)= [b_low,45,90];
par(11,:)= [b_low,45,180];   
par(12,:)= [b_low,45,0];  

par(13,:)= [b_low,-45,270];  
par(14,:)= [b_low,-45,90];
par(15,:)= [b_low,-45,180]; 
par(16,:)= [b_low,-45,0];

par(17,:)= [b_mid,22.5,180];                         
par(18,:)= [b_mid,22.5,0]; 
par(19,:)= [b_mid,22.5,270];                         
par(20,:)= [b_mid,22.5,90];

par(21,:)= [b_mid,-67.5,270];    
par(22,:)= [b_mid,-67.5,90];  
par(23,:)= [b_mid,-67.5,180];    
par(24,:)= [b_mid,-67.5,0];  

par(25,:)= [b_mid,-22.5,180];                       
par(26,:)= [b_mid,-22.5,0]; 
par(27,:)= [b_mid,-22.5,270];                       
par(28,:)= [b_mid,-22.5,90]; 

par(29,:)= [b_mid,67.5,180];                 
par(30,:)= [b_mid,67.5,0];
par(31,:)= [b_mid,67.5,270];                 
par(32,:)= [b_mid,67.5,90];

par(33,:)= [b_high,0,270];
par(34,:)= [b_high,0,90];
par(35,:)= [b_high,0,180];
par(36,:)= [b_high,0,0];

par(37,:)= [b_high,90,180];
par(38,:)= [b_high,90,0];
par(39,:)= [b_high,90,270];
par(40,:)= [b_high,90,90];

par(41,:)= [b_high,45,180];
par(42,:)= [b_high,45,0];
par(43,:)= [b_high,45,270];
par(44,:)= [b_high,45,90];

par(45,:)= [b_high,-45,180];
par(46,:)= [b_high,-45,0];
par(47,:)= [b_high,-45,270];
par(48,:)= [b_high,-45,90];

par(49,:)=[0,0,0];
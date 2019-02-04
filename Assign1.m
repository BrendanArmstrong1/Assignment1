clear
clearvars -GLOBAL
clc
clf
close all

boxwidth = 200e-9; %Width of arena
boxlength =  100e-9; %lengthidth of arena 
k = 1.38064852e-23; %Boltzmann constant
mo = 9.109e-31; %rest mass
m = .26*mo; %mass of electron
T = 500; %temperature in kelvin
sz = 5; %point size for electrons
Numberelectrons = 20;
tau = 0.2e-12;
vth = sqrt(k*T/m); %thermal velocity
dt = 1e-15; %unit time step size
Pscat = 1-exp(-dt/tau);

fracshow = 1; %this section is show showing certain electrons to not crowd the page
Select = round(linspace(1,Numberelectrons,Numberelectrons*fracshow)); % picks a fraction of electrons to show on plot
Show = zeros(1,Numberelectrons);


%define boxes

wallwidth = boxwidth*.2;
wallhight = boxlength/2*.85;


figure(1)
hold on
axis([0 boxwidth 0 boxlength]);
for i=0:5e-10:wallhight
    scatter(boxwidth/2-wallwidth/2,i,15,[0 0 0],'filled')
    scatter(boxwidth/2+wallwidth/2,i,15,[0 0 0],'filled')
    scatter(boxwidth/2-wallwidth/2,boxlength-i,15,[0 0 0],'filled')
    scatter(boxwidth/2+wallwidth/2,boxlength-i,15,[0 0 0],'filled')
end

for i=boxwidth/2-wallwidth/2:5e-10:boxwidth/2+wallwidth/2
    scatter(i,wallhight,15,[0 0 0],'filled')
    scatter(i,boxlength-wallhight,15,[0 0 0],'filled')
end

%---------------- Initialize electron parameters ---------------
for i=1:Numberelectrons
    electrons.R(i) = rand(1,1);
    electrons.G(i) = rand(1,1);
    electrons.B(i) = rand(1,1);
    
    electrons.yposition(i) = rand*boxlength; %random y position from 0 to boxlength
    
    electrons.xposition(i) = rand*boxwidth; %random x position from 0 to boxwidth
    while (electrons.xposition(i) <= boxwidth/2+wallwidth/2 && electrons.xposition(i) >= boxwidth/2-wallwidth/2 && electrons.yposition(i) >= boxlength-wallhight && electrons.yposition(i) <= wallhight)
        electrons.xposition(i) = rand*boxwidth;
        electrons.yposition(i) = rand*boxlength;
    end
    
    
    
   
    electrons.xvel(i) = vth*randn(1,1); %thermal velocity in x direction
    electrons.yvel(i) = vth*randn(1,1); %thermal velocity in y direction
    electrons.xp(i) = electrons.xposition(i)-dt*electrons.xvel(i); %last position of electron required for integration, initialized to zero
    electrons.xpp(i) = electrons.xposition(i) - 2*dt*electrons.xvel(i); % two positions earlier
    electrons.yp(i) = electrons.yposition(i)-dt*electrons.yvel(i); 
    electrons.ypp(i) = electrons.yposition(i) - 2*dt*electrons.yvel(i);
end
% --------------------------------------------------------------
t=0;
for i=1:length(Select)
    Show(Select(i)) = 1;
end

while t<1
    figure(1);
    scatter(electrons.xposition.*Show,electrons.yposition.*Show,sz,[electrons.R' electrons.G' electrons.B']);
    t = t+dt; % increase time step
    
    
    
    
   
    
    leftcheck = -electrons.xpp + 2*electrons.xp <= 0; % find elements of xposition that are less than zero
    rightcheck = -electrons.xpp + 2*electrons.xp >= boxwidth; % find elements larger than box width
    
    if sum(leftcheck) > 0
        electrons.xposition(leftcheck == 1) = electrons.xposition(leftcheck == 1) + boxwidth;
        electrons.xp(leftcheck == 1) = electrons.xp(leftcheck == 1) + boxwidth;
        electrons.xpp(leftcheck == 1) = electrons.xpp(leftcheck == 1) + boxwidth;
        leftcheck(leftcheck == 1) = 0;
    end
    if sum(rightcheck) > 0
        electrons.xposition(rightcheck == 1) = electrons.xposition(rightcheck == 1) - boxwidth;
        electrons.xp(rightcheck == 1) = electrons.xp(rightcheck == 1) - boxwidth;
        electrons.xpp(rightcheck == 1) = electrons.xpp(rightcheck == 1) - boxwidth;
        rightcheck(rightcheck == 1) = 0;
    end
    
    topcheck = -electrons.ypp + 2*electrons.yp >=boxlength;
    bottomcheck = -electrons.ypp + 2*electrons.yp <= 0;
    
    if sum(topcheck) > 0
        electrons.yp(topcheck == 1) = electrons.yp(topcheck == 1) + 2*(boxlength - electrons.yp(topcheck == 1));
        electrons.ypp(topcheck == 1) = electrons.ypp(topcheck == 1) + 2*(boxlength - electrons.ypp(topcheck == 1));
    end
    if sum(bottomcheck) > 0
        electrons.yp(bottomcheck == 1) = electrons.yp(bottomcheck == 1) - 2*electrons.yp(bottomcheck == 1);
        electrons.ypp(bottomcheck == 1) = electrons.ypp(bottomcheck == 1) - 2*electrons.ypp(bottomcheck == 1);
    end
    averagevel2 = sum(electrons.xvel.^2+electrons.yvel.^2)/Numberelectrons;
    clc
    temp = .5*m*averagevel2*2/2/k
    
    
    scat = rand(1,Numberelectrons) <= Pscat;
    electrons.xvel(scat) = vth*randn(1,length(electrons.xvel(scat)));
    electrons.yvel(scat) = vth*randn(1,length(electrons.yvel(scat)));
    electrons.xp(scat) = electrons.xposition(scat)-dt*electrons.xvel(scat); 
    electrons.xpp(scat) = electrons.xposition(scat) - 2*dt*electrons.xvel(scat);
    electrons.yp(scat) = electrons.yposition(scat)-dt*electrons.yvel(scat); 
    electrons.ypp(scat) = electrons.yposition(scat) - 2*dt*electrons.yvel(scat);
    
    
    electrons.xposition = -electrons.xpp + 2*electrons.xp; 
    electrons.xpp = electrons.xp;
    electrons.xp = electrons.xposition;
    
    electrons.yposition = -electrons.ypp + 2*electrons.yp;    
    electrons.ypp = electrons.yp;
    electrons.yp = electrons.yposition;
    
    
    
    
    pause(.001);
end
clear all, close all;
clc;

%% generate initial conditions
%% trial 1 - uniform distribution with hard circular boundaries
L = 1;% Galaxy size
N = 5000; % number of stars
randVec = zeros(1,L);
k=1;
while k<=N % use a while loop to ensure that we get L values even some random numbers might lay outside the unit circle
    randVec(1,k) = (rand(1,1)-0.5)*2;
    randVec(2,k) = (rand(1,1)-0.5)*2;
    if norm(randVec(:,k)) <= 1 % check if the generated random number lays within the unit circle
        k = k+1;
    end
end

% mapping random numbers onto square
x = L*randVec(1,:);
y = L*randVec(2,:);
r = sqrt(x.^2 + y.^2);
phi = atan2(y,x);


%% plotting
close all
figure(1)
plot(r.*cos(phi),r.*sin(phi),'.w','Markersize',10)
% plot(x,y,'.w','Markersize',10)
set(gca,'color','k')
xlim([-L,L])
ylim([-L,L])
axis equal


%% velocity distribution function
close all;
d = 0:1/1000:1;
calcOmega = @(x) 1./(x*10+1) .* ( 1 - 1./((2*x*10).^2+1) );
plot(d,calcOmega(d))


%% every particle gets a velocity according to the velocity distribution function
for k=1:N
    omega(k) = calcOmega(r(k));
end


%% simulating a few time steps
N_steps = 200;
Ts = 0.1;

close all
figure('color','k','position',[958    42   962   954])
plotHandle = plot(r.*cos(phi),r.*sin(phi),'.w','Markersize',1);
% plot(x,y,'.w','Markersize',10)
set(gca,'color','k')
xlim([-L,L])
ylim([-L,L])
axis equal
hold all
for k=1:N_steps
	phi = phi + Ts*omega; % euler forward
    set(plotHandle,'XData',r.*cos(phi),'YData',r.*sin(phi))
    drawnow
    pause(0.05)
end

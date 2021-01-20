clear all, close all;
clc;

%% parameters
L = 1;% Galaxy size
N = 2000; % number of stars
starDistr = 2;  % start distribution:  1=uniform, 2=gaussian

centerDensity = 1; % only for starDistr=2

%% generate initial conditions
if starDistr == 1
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
    
elseif starDistr == 2
    randVec(1,:) = randn(1,N)/3;
    randVec(2,:) = rand(1,N)*2*pi;

    r = randVec(1,:)/centerDensity*L;
    phi = randVec(2,:);
    subplot(211); hist(r,100);
    subplot(212); hist(phi,100)
    
else
    error('Wrong type of star distribution')
end



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
omega = zeros(1,N);
for k=1:N
    omega(k) = calcOmega(abs(r(k)));
end


%% simulating a few time steps
N_steps = 200;
Ts = 0.1;
L_plot = max(r);

close all
figure('color','k','position',[958    42   962   954])
plotHandle = plot(r.*cos(phi),r.*sin(phi),'.w','Markersize',1);
% plot(x,y,'.w','Markersize',10)
set(gca,'color','k')
xlim([-L_plot,L_plot])
ylim([-L_plot,L_plot])
axis equal
hold all
for k=1:N_steps
	phi = phi + Ts*omega; % euler forward
    set(plotHandle,'XData',r.*cos(phi),'YData',r.*sin(phi))
    xlim([-L_plot,L_plot])
    ylim([-L_plot,L_plot])
    drawnow
    pause(0.05)
end

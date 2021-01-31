clear all, close all;
clc;

%% parameters
L = 1;% Galaxy size
N = 5000; % number of stars
starDistr = 2;  % start distribution:  1=uniform, 2=gaussian
speedFunction = 2; % 1=Kelper type speed function, 2=real star speed function
N_animSteps = 100;
Ts = 0.06;
centerDensity = 1; % only for starDistr=2
armDensity = 0.4; % [0,1]
armProm = 10; % arm prominence
N_arms = 3;

exportType = 1; % export to 1=GIF, 2=mpeg4

%% parameters for export to MPEG4 or GIF 

if exportType==1
    gifname = 'galaxyAnimation.gif';
    createGIFimage = true;
    idx_start = 1;
elseif exportType==2
    
end

%% generate initial conditions -> distribute starts

if starDistr == 1           % stars distributed uniformally over radius and angle
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
    
elseif starDistr == 2           % starts distributed in a sinusoidal fashion to generate arms
    % Idea: 
    %   1) generate a uniform distribution in radial direction -> randVec(1,:)
    %   2) generate a sinusoidal distribution in angular direction -> randVec(2,:)
    randVec(1,:) = abs(randn(1,N)/3);   % uniform distribution in radial direction
    help = 2*pi*rand(1,round(N*(1-armDensity))); % start with a uniform distribution in angular direction
    for k=1:N_arms-1
        % depending how many arms there are, generate a map from [0,2*pi) onto [-1,1]  
        help = [help, (1/armProm*asin(2*(rand(1,round(N*armDensity/N_arms))-0.5))+ pi/2 + 2*pi/N_arms*(k-1))];
    end
    help = [help, (1/armProm*asin(2*(rand(1,round(N*armDensity/N_arms))-0.5))+pi/2+2*pi/N_arms*(N_arms-1))];
    try
        randVec(2,:) = [help, (1/armProm*asin(2*(rand(1,round(N*armDensity/N_arms))-0.5))+pi/2+2*pi/N_arms*(N_arms-1))];
    catch
        help =         [help, (1/armProm*asin(2*(rand(1,round(N*armDensity/N_arms))-0.5))+pi/2+2*pi/N_arms*(N_arms-1))];
        if numel(help) > N
            help(N+1:end) = [];
            randVec(2,:) = help;
        elseif numel(help) < N
            help(end+1:end+(N-numel(help))) = 2*pi*rand(1,N-numel(help));
            randVec(2,:) = help;
        end
    end

    r = randVec(1,:)/centerDensity*L;
    phi = randVec(2,:);
%     subplot(211); hist(r,100);
%     subplot(212); hist(phi,100); xlim([0,2*pi])
    
else
    error('Wrong type of star distribution')
end



%% plotting
close all
figure(1)
plot(r.*cos(phi),r.*sin(phi),'.w','Markersize',1)   % using polar coordinates
% plot(x,y,'.w','Markersize',10)                    % using cartesian coordinates
set(gca,'color','k')
xlim([-L,L])
ylim([-L,L])
axis equal


%% velocity distribution function
close all;
d = 0:1/1000:1;
if speedFunction == 1   % Kelper type speed function
    calcOmega = @(x) 1/0.4 * 1./(x*10+1) .* ( 1 - 1./((2*x*10).^2+1) );
else                    % real star speed function
    calcOmega = @(x) 0.25*(log10(x+0.01)+4) - 0.3;
end
plot(d,calcOmega(d)), grid on;


%% every particle gets a velocity according to the velocity distribution function
omega = zeros(1,N);
for k=1:N
    omega(k) = calcOmega(abs(r(k)));
end


%% simulating a few time steps
L_plot = max(r);

close all
figure('color','k','position',[500          42        640         480])
plotHandle = plot(r.*cos(phi),r.*sin(phi),'.w','Markersize',1);
% plot(x,y,'.w','Markersize',10)
set(gca,'color','k')
xlim([-L_plot,L_plot])
ylim([-L_plot,L_plot])
axis equal
hold all
for k=1:N_animSteps
	phi = phi + Ts*omega; % euler forward
    set(plotHandle,'XData',r.*cos(phi),'YData',r.*sin(phi))
    xlim([-L_plot,L_plot])
    ylim([-L_plot,L_plot])
    drawnow
    pause(0.05)
    
    exportImagesToGif( getframe(gcf), 'collect' );
%     if exportType==1
%         f = getframe(gcf);
%         if k==idx_start
%             [im,map] = rgb2ind(f.cdata,256,'nodither');
%         else
%             im(:,:,1,k-idx_start+1) = rgb2ind(f.cdata,map,'nodither');
%         end
%     end
end

[ gifFile ] = exportImagesToGif( getframe(gcf), 'finalize', 'Ts', Ts, 'loopCnt', inf, 'gifName', gifname );
% if exportType==1
%     %% export simulation to an animated gif
%     % save to fle
%     imwrite(im,map,gifname,'DelayTime',Ts,'LoopCount',inf) %g443800
%     alreadysaved = 1;
%     fprintf('done.\n');
% end

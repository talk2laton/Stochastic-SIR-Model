function SIRModelSimulation(Np, Ni, Ds, Prinf, Tr, F, videoname)
% This file simulates SIR model for spread of infectious diseases
close all;
%     Np     = the population of the community. (Default = 100, so you easily  
%              read the percentage of the community infectious/infected)
%     Ni     = the number of index case (initial infecious individual, Default = 1).
%     Ds     = is the minimum safe distance below which the probability of being 
%              infected is Prinf (Measure of how contagious the disease is: 
%              Lower means more contagious. Default = 0.05).
%     Prinf  = probability that you will be infected if you are closer than Ds 
%              to an infectious iindividial (Measure of how precautious people 
%              are: washing hands, not touching mucus membrane. Default = 0.2).
%     Tr     = Time it takes for an individual to be removed from the population 
%              (Death, Quarantined, Recovered and become immuned. Default = 1).
%     F      = Relative repulsive force strength (measure of social distancing Default = 1)
% videoname  = Simulation video file name
if(nargin < 1)
    Np = 100; Ni = 1; Ds = 0.05; Prinf = 0.2; Tr = 1; F = 1; 
elseif (nargin < 2) 
    Ni = 1; Ds = 0.05; Prinf = 0.2; Tr = 1; F = 1;
elseif (nargin < 3) 
    Ds = 0.05; Prinf = 0.2; Tr = 1; F = 1; 
elseif (nargin < 4) 
    Prinf = 0.2; Tr = 1; F = 1;
elseif (nargin < 5) 
    Tr = 1; F = 1; 
elseif (nargin < 6) 
    F = 1; 
end
writevideo = exist('videoname','var');
Colors = 'grb'; dt = 1/30; S = (1:Np)'; status = ones(Np,1); t = 0;
I = randi(Np,Ni,1); T = zeros(size(I)); S(I) = []; status(I) = 2; 
% Np+1 particle is dummy removed particle placed outside the domain, just to force removed to appear in the legend 
Position = [rand(Np,2);-1,-1]; status(Np+1) = 3; R = Np+1; 
% To make sure that the velocity of each particle = 1 makesure dummy doesn't move;
Velocity = rand(Np+1,2)-0.5; Vnorm = vecnorm(Velocity,2,2); Velocity = Velocity./[Vnorm, Vnorm]; Velocity(end,:) = 0;
figure('Position',[200 500 1300 400],'Color', [1.00 1.00 1.00]); 
ax1 = axes('Position',[0.05 0.1 0.4 0.8]);  ax2 = axes('Position',[0.55,0.1,0.4,0.80]); 
axes(ax1); axis([0,10*Tr,0,Np]); hold on; box on; ax1.TickLabelInterpreter = 'Latex';
infectious = numel(I); susceptible = numel(S); removed = numel(R)-1;
S1 = fill([0, t, t, 0, 0],[0, 0, Np, Np, 0],Colors(1));
I1 = fill([t;flipud(t)],[zeros(size(t));flipud(infectious)],Colors(2));
R1 = fill([t;flipud(t)],[zeros(size(t));flipud(removed)],Colors(3)); 
lgd = legend([S1,I1, R1],{'Susceptible', 'Infectious', 'Removed'},'FontSize',12, 'interpreter','latex', 'location', 'northeastoutside'); 
title(lgd,'Legend','FontSize',12, 'interpreter','latex');
axes(ax2); axis([0,1,0,1]); hold on; ax2.XTick = []; ax2.YTick = []; box on; 
S2 = plot(Position(status == 1,1), Position(status == 1,2), 'o', 'MarkerEdgeColor','k', 'MarkerFaceColor', Colors(1));
I2    = plot(Position(status == 2,1), Position(status == 2,2), 'o', 'MarkerEdgeColor','k', 'MarkerFaceColor', Colors(2));
R2     = plot(Position(status == 3,1), Position(status == 3,2), 'o', 'MarkerEdgeColor','k', 'MarkerFaceColor', Colors(3));
lgd = legend([S2, I2, R2],{'Susceptible', 'Infectious', 'Removed'}, 'FontSize',12, 'interpreter','latex', 'location', 'northeastoutside'); 
title(lgd,'Legend','FontSize',12, 'interpreter','latex');
drawnow; 
if(writevideo)
    vidObj = VideoWriter(videoname); set(vidObj, 'FrameRate',10); open(vidObj); image = getframe(gcf); writeVideo(vidObj,image);
end
while(t < 10*Tr)
    %% Accounting for Repulsion (Social distancing)
    [X1,X2] = meshgrid(Position(:,1));  [Y1,Y2] = meshgrid(Position(:,2));
    Xc = X1 - X2; Yc = Y1 - Y2; D = ((X1 - X2).^2 + (Y1 - Y2).^2).^(-1.5); 
    %Taking out the forces from the dummy element and interraction of a particle with itself
    D(1:Np+2:end) = 0; D(:, end);                         
    Acceleration = -F*1e-2*[sum(Xc.*D,2),sum(Yc.*D,2)]; Acceleration(end,:) = 0; 
    %% Updating Position and Velocity
    Position = Position + dt*Velocity; T = T + dt; 
    Velocity = Velocity + Acceleration*dt; Vnorm = vecnorm(Velocity,2,2);
    Velocity(Position<0) = -Velocity(Position<0); Position(Position<0) = 0;
    Velocity(Position>1) = -Velocity(Position>1); Position(Position>1) = 1;
    Velocity(Vnorm > 1,:) = Velocity(Vnorm > 1,:)./[Vnorm(Vnorm > 1), Vnorm(Vnorm > 1)];
    Position(end,:) = -1; Velocity(end,:) = 0; Velocity(isnan(Velocity)) = 0;
    %% Distance between infectious and susceptible person
    [Xs,Xi] = meshgrid(Position(S,1),Position(I,1));  [Ys,Yi] = meshgrid(Position(S,2),Position(I,2));
    d = sqrt((Xs - Xi).^2 + (Ys - Yi).^2); D = min([d;d]);  
    dn = sum(D < Ds); 
    %% Probability (How people takes preventive measures like handwashing, not touching any mucus membrane)
    P = ones(size(D)); P(D < Ds) = rand(1,dn);
    %% Updating S and I;
    I = [I;S(D < Ds & P < Prinf)]; T = [T;zeros(size(S(D < Ds & P < Prinf)))]; 
    status(S(D < Ds & P < Prinf)) = 2;  S(D < Ds & P < Prinf) = []; 
    %% Updating R and I;
    R = [R;I(T > Tr)]; status(I(T > Tr)) = 3; I(T > Tr) = []; T(T > Tr) = [];
    %% Updating graphics
    infectious = [infectious;numel(I)]; susceptible = [susceptible;numel(S)]; removed = [removed;numel(R)-1]; t = [t;t(end)+dt];
    S1.XData = [0, t(end), t(end), 0, 0]; S1.YData = [0, 0, Np, Np, 0];
    I1.XData = [t;flipud(t)]; I1.YData = [zeros(size(t));flipud(infectious)];
    R1.XData = [t;flipud(t)]; R1.YData = [Np - zeros(size(t));Np - flipud(removed)]; 
    S2.XData = Position(status == 1,1); S2.YData = Position(status == 1,2);
    I2.XData = Position(status == 2,1); I2.YData = Position(status == 2,2);
    R2.XData = Position(status == 3,1); R2.YData = Position(status == 3,2); drawnow; 
    if(writevideo)
        image = getframe(gcf); writeVideo(vidObj,image);
    end
end
if(writevideo)
    close(vidObj);
end
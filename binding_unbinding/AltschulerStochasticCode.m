% Sstochastic Altschuler model for polarization along cell membrane
%
% Diffusion of X(t) particles on the cell membrane
%
% Three types of chemical reactions can happen on the membrane:
% (1) association: a particle will randomly appear on the membrane 
% (2) disassociation: a particle will randomly dissapear from the membrane
% (3) association thru recruitment: a particle will recruit another
% particle to the membrane
% 
% Last updated: 4/29/2019

clear
close all;

N       = 100;                  % total number of molecules in the cell (conserved)
kon     = 0.1;                  % units: 1/min
kfb     = 1.0;                  % units: 1/min
koff    = 0.9;                  % units: 1/min
D       = 0.012;                % units: um^2/min

Tend = 6.0;                     % end simulation time
p    = 0.5;                     % probability of hoping left or right
dt   = 0.001;                   % temporal discretization (units: min)
Nt   = Tend/dt;                 % total time steps
dx   = sqrt(2*D*dt);            % spatial discretization (units: um)
L    = 10.0;                    % length of membrane box (units: um)

MAX_OUTPUT_LENGTH = 10000;
pos = zeros(N,Nt);                      % array of positions at all time levels
n   = zeros(N,Nt);                      % state of the particle (0 inactive, 1 active)
T = zeros(MAX_OUTPUT_LENGTH,1);         % times of chemical reactions
X = zeros(MAX_OUTPUT_LENGTH,1);         % number of molecules on the membrane    

% initial conditions (T=0)
rxn_count       = 1;
X(1)            = 0.1*N;                % # of particles on membrane
T(1)            = 0.0;                   
pos(1:X(1),1)   = L*rand(X(1),1);       % random locations for mem-bound particles
n(1:X(1),1)     = 1;                    % activate mem-bound particles

% Loop in time
%

% (1) Find number of bound molecules X(t) via a Markov process and
% corresponding chemical reaction times T(t)
%
while T(rxn_count)<Tend
    
  % Sample earliest time-to-fire (tau)
  r = rand(2,1);
  nx = X(rxn_count);
  a0 = koff*nx + (kon+kfb*nx/N)*(N-nx);
  tau = -log(r(1))/a0;
  
  % Update the number of molecules on the membrane
  X(rxn_count+1,1) = nx + (r(2)<((kon+kfb*nx/N)*(N-nx)/a0))*1.0 + (r(2)>=((kon+kfb*nx/N)*(N-nx)/a0))*(-1.0);
  T(rxn_count+1)   = T(rxn_count) + tau;
  rxn_count = rxn_count + 1;
end
totalrxns = rxn_count;

figure(1);
scatter(T(1:totalrxns),X(1:totalrxns));

% (2) Perform Brownian motion for all bound molecules between consecutive 
% chemical reactions
%
% X = [10 11 12 11 10];
% T = [0 1 5 8 9];
% TT = [T,Tend];
% XX = [X X(end)];
t = 1;
tic
%for i=1:5
for i=1:totalrxns
    K = X(i);
    
    % Between reactions, perform Brownian motion with periodic BC
    while(t*dt<T(i+1) && t<Nt)
        r = rand(K,1);    % coin flip
        n(1:K,t+1) = 1;
        pos(1:K,t+1) = pos(1:K,t) + dx*((r<p)*1.0 + (r>(1-p))*(-1.0));

        % periodic boundary conditions
        for j=1:K
            pos(j,t+1) = pos(j,t+1) + (-L).*(pos(j,t+1)>L) + (L).*(pos(j,t+1)<0.0);
        end
        
        t = t+1;
    end
    
    % Setup for next chemical rxn
    pon = kon/(kon+kfb*(N-K));
    heq = K/N;
    
    if(X(i+1)<K && i<totalrxns)                      % diassociation event (particle off)
        id = randi([1,K],1);
        oldcol = pos(id,1:end);
        othercols = pos([1:id-1,id+1:K],1:end);
        otherothercols = pos(K+1:end,1:end);
        newpos = [othercols;oldcol;otherothercols];
        pos = newpos;
        n(K,t) = 0;
    elseif(X(i+1)>K && i<totalrxns)                 % association event (on or recruitment)
        rr = rand(1,1); 
        id = randi([1,K],1);
        pos(K+1,t) = pos(K+1,t)+(rr<pon)*L*rand(1,1);   % on event
        pos(K+1,t) = pos(K+1,t)+(rr>=pon)*pos(id,t);    % recruitment event
        n(K+1,t) = 1;
    end
     
end
toc

% plot all particle trajectories
%
cc = [255 219 88]/256.*ones(Nt,1); % color string
time = linspace(0,Tend,Nt);
figure(2);
for j=1:max(X)
    hold on;
    %plot(linspace(0,Tend,Nt),pos(j,:));
    %scatter(linspace(0,Tend,Nt),pos(j,:),1,ccc);
    scatter(linspace(0,Tend,Nt),pos(j,:),1,cc);
    %plot(linspace(0,Tend,Nt),pos(j,:));
    box on;
    set(gca,'Color','k','fontsize',30);
    xlabel('Time (minutes)');
    ylabel('Membrane position');
    ylim([0 10])
%     plot(1*ones(100,1),linspace(0,10,100),'-r');
%     plot(5*ones(100,1),linspace(0,10,100),'-r');
%     plot(8*ones(100,1),linspace(0,10,100),'-r');
%     plot(9*ones(100,1),linspace(0,10,100),'-r');
end




% *** THIS PROGRAM DEALS W/ ONE CHEMICAL SPECIES WITH A BINDING AND UNBINDING REACTION, WITH NO FEEDBACK ***

% rates: r_bind = k_bind * (1-a)   
% 			 r_unbind = k_unbind * a

% Initialize the Basic Parameters: number of molecules of A (N), max # of reactions (T) , state vector (X), proportion bound over time (A), Reaction Coeffecients (k_bind, k_unbind) L is the length of the 1D representation of the membrane of a cell.
close all;

N = 1000; T = 400; k_bind = 0.5; k_unbind = 0.4; L = 1;  

% X(1, :) is a bool representing whether each particle is bound or not.
% X(2, :) is a float representing the position of each particle.
X = zeros(2,N); 
X(2,:) = rand(L, N);
% A(1, :) is a float representing the proportion of chemical A that is bound.
% A(2, :) is the time that this value is recorded.
A = zeros(2,T);
big_X = zeros(2,N,T);

function particle_index = find_bound_particle(X, N)
x = randi([1,N]);
switch X(1, x)
    case 1
        particle_index = x;
    case 0
        particle_index = find_bound_particle(X, N);
end
end

function particle_index = find_unbound_particle(X, N)
x = randi([1,N]);
switch X(1, x)
    case 0
        particle_index = x;
    case 1
       particle_index = find_unbound_particle(X, N);
end

end
% This odd index is due to using i-1 to gather the value of the last index in each of the arrays, and I didn't want to add a switch statement every time for i = 1. This is more computationally efficient
for i = 2:(T+1)

    % Setting the rates at the start of each loop
    r_bind = k_bind*(1-A(1,i-1));
    r_unbind = k_unbind*A(1,i-1);
    r_sum = r_bind + r_unbind;

    % Generating a random number for use in Gillespie Algo
    w = rand(1);

    tau = (1 / r_sum) * log(1 / w);

    % Picking the reaction using Gillespie Algo
    if (r_sum * w < r_bind)
        X(1,find_unbound_particle(X, N)) = 1;
    else
        X(1,find_bound_particle(X, N)) = 0;
    end

    % Updating proportion and time for A
    A(:, i) = [(sum(X(1, :)) / N), (A(2, i-1) + tau)];
    big_X(:,:,i) = X;
end

% The steady state  f the concentration for the ODE version of this is:
a_steady_state = (k_bind / (k_bind + k_unbind));


%Setting up the ODE to plot 
[t, a] =  ode45(@(t, a) k_bind * (1 - a) - k_unbind * a, [0, A(2,T)], A(1,1));

subplot(2,1,1);
%ODE plot
plot(t,a); hold on;

%Stochastic Plot
plot(A(2,:),A(1,:));

%Steady state comparison
yline(a_steady_state);

%Title and labels
title('Simple Binding-Unbinding Model w/ ODE Comparison');
xlabel('Time');
ylabel('Bound Proportion of Chemical A');

%Plotting a histogram of particle location
subplot(2,1,2);
hist(X(2,:),20);
title('Histogram of Particle Position at the Ending Time');
xlabel('Position of Particles');
ylabel('Number of Particles');


%% heatmap
pos_mat = squeeze(big_X(2,:,:));
pos_mat = discretize(pos_mat,20)./20; % discretize 1D domain

bool_mat = logical(squeeze(big_X(1,:,:)));

% create table
particle_vec = repmat([1:1:N]',T+1,1);
time_vec = reshape([1:1:(T+1)] .* ones(1,N)',[],1);
pos_vec = reshape(pos_mat, N*(T+1),[]); % reshape position
bool_vec = reshape(bool_mat, N*(T+1),[]); % reshape bind/unbind bool

tab = table(particle_vec,time_vec,pos_vec,bool_vec);

bind_tab = tab(bool_vec == true,:);
unbind_tab = tab(bool_vec == false & time_vec>1 , :); % avoid inital condition biasing overall

% heatmap
figure
subplot(2,1,1)
h = heatmap(bind_tab,"time_vec","pos_vec", Colormap=hot(100));
title('Binded Kymograph'); xlabel('Time'); ylabel('Position')

subplot(2,1,2)
h1 = heatmap(unbind_tab,"time_vec","pos_vec", Colormap=hot(100));
title('Unbinded Kymograph'); xlabel('Time'); ylabel('Position')



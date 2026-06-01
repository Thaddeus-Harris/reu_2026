% *** THIS PROGRAM DEALS W/ ONE CHEMICAL SPECIES WITH A BINDING AND UNBINDING REACTION, WITH NO FEEDBACK ***

% rates: r_bind = k_bind * (1-a)
% 			 r_unbind = k_unbind * a

N = 100; % The number of molecules
T = 8000; % Number of reactions
k_bind = 0.000001; % Binding rate constant
k_unbind = 0.9; % Unbinding rate constant
k_feedback = 1.0; % Feedback rate constant
L = 1; % Length of cell membrane
v_x = 0.0005; % Velocity for random walk
p = 0.5; % Probability of moving to the right in random walk (0.5 for no bias)

% X(1,:) is a bool representing whether each particle is bound or not.
X(1,:,:) = zeros(1,N,T);
X(1,1,1) = 1;
% X(2,:) is a float representing the position of each particle.
X(2,:,1) = (L * rand(1, N));

% A(1, :) is a float representing the proportion of chemical A that is bound.
% A(2, :) is the time that this value is recorded.
A = zeros(2,T);
tic
W_1 = rand(1,T+1);
W_2 = rand(1,T+1);
% This odd index is due to using i-1 to gather the value of the last index in each of the arrays, and I didn't want to add a switch statement every time for i = 1. This is more computationally efficient
for i = 2:(T+1)

	% set the position to the position of the previous time.
	X(:,:,i) = X(:,:,i-1);

	% Setting the rates at the start of each loop
	r_bind = k_bind*(N-sum(X(1,:,i)));
	r_unbind = k_unbind*sum(X(1,:,i));
	r_feedback = k_feedback*A(1,i-1)*(N - sum(X(1,:,i)));
	r_sum = r_bind + r_unbind + r_feedback;

	% Generating a random number for use in Gillespie Algo

	tau = (1 / r_sum) * log(1 / W_1(i));

	% Picking the reaction using Gillespie Algo
  if (r_sum * W_2(i) < r_bind)
		unbound_list = find(X(1,:,i) == 0);
		unbound_index = unbound_list(randi(length(unbound_list)));
		X(1,unbound_index,i) = 1;
		X(2,unbound_index,i) = L * rand(1);
	elseif (r_bind < r_sum * W_2(i) <= r_bind + r_unbind)
		bound_list = find(X(1,:,i) == 1);
		bound_index = bound_list(randi(length(bound_list)));
		X(1,bound_index,i) = 0;
	else 
		unbound_list = (X(1,:,i) == 0);
		unbound_index = randi(length(unbound_list));
		
		bound_list= (X(1,:,i) == 1);
		bound_index = randi(length(bound_list));

		X(1,unbound_index,i) = 1;
		X(2,unbound_index,i) = X(2,bound_index,i);
	end

	% Creates masks for bound & unbound values
	unbound = (X(1,:,i) == 0);
	bound = (X(1,:,i) == 0);

	X(2,unbound,i) = X(2,unbound,i-1);
	random_steps = ((2 * randi(2, 1, length(X(2,bound,i)))) - 3) * v_x * tau;
	X(2,bound,i) = mod(X(2,bound,i-1) + random_steps, L);

	% Updating proportion and time for A
	A(:, i) = [(sum(X(1, :,i)) / N), (A(2, i-1) + tau)];
end
toc

% The steady state  f the concentration for the ODE version of this is:
a_steady_state = 1 - (k_unbind/k_feedback);


%Setting up the ODE to plot
[t, a] =  ode45(@(t, a) k_bind * (1 - a) + k_feedback * (1 - a) * a  - k_unbind * a, [0, A(2,T)], A(1,1));

subplot(1,3,1);
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
subplot(1,3,2);
hist(X(2,:,T),20);

title('Histogram of Particle Position at the Ending Time');
xlabel('Position of Particles');
ylabel('Number of Particles');

subplot(1,3,3);

for n = 1:N
	scatter(A(2,:),X(2,n,:),"c",".");
	hold on;
end
title('Particle Trajectory');
ylabel('X Position');
xlabel = ('time');

% *** THIS PROGRAM DEALS W/ ONE CHEMICAL SPECIES WITH A BINDING AND UNBINDING REACTION, WITH NO FEEDBACK ***

% rates: r_bind = k_bind * (1-a)
% 			 r_unbind = k_unbind * a

% Initialize the Basic Parameters: number of molecules of A (N), max # of reactions (T) , state vector (X), proportion bound over time (A), Reaction Coeffecients (k_bind, k_unbind) L is the length of the 1D representation of the membrane of a cell.

N = 100; T = 4000; k_bind = 0.1; k_unbind = 0.5; k_feedback = 0.3; L = 10; v_x = 0.01;

% X(1,:) is a bool representing whether each particle is bound or not.
% X(2,:) is a float representing the position of each particle.
X(1,:,:) = zeros(1,N,T);
X(2,:,1) = (L * rand(1, N));
% A(1, :) is a float representing the proportion of chemical A that is bound.
% A(2, :) is the time that this value is recorded.
A = zeros(2,T);

% This odd index is due to using i-1 to gather the value of the last index in each of the arrays, and I didn't want to add a switch statement every time for i = 1. This is more computationally efficient
for i = 2:(T+1)

	% set the position to the position of the previous time.
	X(:,:,i) = X(:,:,i-1);

	% Setting the rates at the start of each loop
	r_bind = k_bind*(N-sum(X(1,:,i)));
	r_unbind = k_unbind*sum(X(1,:,i));
	r_feedback = k_feedback*A(1,i-1)*sum(X(1,:,i));
	r_sum = r_bind + r_unbind + r_feedback;

	% Generating a random number for use in Gillespie Algo
	w = rand(1);

	tau = (1 / r_sum) * log(1 / w);

	% Picking the reaction using Gillespie Algo
  if (r_sum * w < r_bind)
		unbound_list = find(X(1,:,i) == 0);
		unbound_index = unbound_list(randi(length(unbound_list)));
		X(1,unbound_index,i) = 1;
	elseif (r_sum * w < r_bind + r_unbind)
		bound_list = find(X(1,:,i) == 1);
		bound_index = bound_list(randi(length(bound_list)));
		X(1,bound_index,i) = 0;
		X(2,bound_index,i) = 10 * rand(1);
	else 
		unbound_list = find(X(1,:,i) == 0);
		unbound_index = unbound_list(randi(length(unbound_list)));
		bound_list= find(X(1,:,i) == 1);
		bound_index = bound_list(randi(length(bound_list)));
		X(1,unbound_index,i) = 1;
		X(2,unbound_index,i) = X(2,bound_index,i);
	end



		for n = 1:N
			switch X(1,n,i)
				case 0
          X(2,n,i) = X(2,n,i-1);
				case 1
		      X(2,n,i) = mod( (X(2,n,i-1) + (2 * randi(2) - 3) * v_x * tau), L);
			end
	end
	% Updating proportion and time for A
	A(:, i) = [(sum(X(1, :,i)) / N), (A(2, i-1) + tau)];
end

% The steady state  f the concentration for the ODE version of this is:
a_steady_state = (k_bind / (k_bind + k_unbind));


%Setting up the ODE to plot
[t, a] =  ode45(@(t, a) k_bind * (1 - a) - k_unbind * a, [0, A(2,T)], A(1,1));

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

% *** THIS PROGRAM DEALS W/ ONE CHEMICAL SPECIES WITH A BINDING AND UNBINDING REACTION, WITH NO FEEDBACK ***

% rates: r_bind = k_bind * (1-a)
% 			 r_unbind = k_unbind * a

N_x = 5000; % The number of molecules
N_x_starting_bound = N_x/10; % The number of molecules
N_y = 5000; % The number of molecules
N_y_starting_bound = N_y/10; % The number of molecules
T = 5000; % Number of reactions
k_bind = 0.00000005; % Binding rate constant
k_unbind = 0.9; % Unbinding rate constant
k_feedback = 1.0; % Feedback rate constant
L = 1; % Length of cell membrane
a_r = 0.00001; % Acceleration as experience by x particles in random walk
v_r = 0.0025; % Velocity for random walk
p = 0.5; % Probability of moving to the right in random walk (0.5 for no bias)

% X(1,:,t) is a bool representing whether each particle is bound or not.
X(1,:,:) = zeros(1,N_x,T);
X(1,1:N_x_starting_bound,1) = 1;
Y(1,:,:) = zeros(1,N_y,T);
Y(1,1:(N_y_starting_bound),1) = 1;
Y(1,1:N_y_starting_bound,1) = 1;
% X(2,:,t) is a float representing the x position of each particle.
% X(3,:,t) is a float representing the y position of each particle.
X(2:3,1:N_x_starting_bound,1) = (L * rand(2,N_x_starting_bound));
X(2:3,N_x_starting_bound+1:N_x,1) = nan;
Y(2:3,1:N_y_starting_bound,1) = (L * rand(2,N_y_starting_bound));
Y(2:3,N_y_starting_bound+1:N_y,1) = nan;

X(4:5,N_x_starting_bound+1:N_x,1) = 0;

% A(1, :) is a float representing the proportion of chemical X that is bound.
% A(2, :) is the time that this value is recorded.
% B(1, :) is a float representing the proportion of chemical X that is bound.
% B(2, :) is the time that this value is recorded.
A = zeros(2,T);
B = zeros(2,T);
W_1 = rand(1,T+1);
W_2 = rand(1,T+1);
W_3 = rand(1,T+1);
W_4 = rand(1,T+1);
% This odd index is due to using i-1 to gather the value of the last index in each of the arrays, and I didn't want to add a switch statement every time for i = 1. This is more computationally efficient
for i = 2:(T+1)
	% set the position to the position of the previous time.
	X(1,:,i) = X(1,:,i-1);
	X(2,:,i) = X(2,:,i-1);
	X(3,:,i) = X(3,:,i-1);
	X(4,:,i) = X(4,:,i-1);
	X(5,:,i) = X(5,:,i-1);

	% Setting the rates at the start of each loop
	r_bind = k_bind*(N_x-sum(X(1,:,i)));
	r_unbind = k_unbind*sum(X(1,:,i));
	r_feedback = k_feedback*A(1,i-1)*(N_x - sum(X(1,:,i)));
	r_sum = r_bind + r_unbind + r_feedback;

	% Generating a random number for use in Gillespie Algo

	tau = (1 / r_sum) * log(1 / W_1(i));

	% Picking the reaction using Gillespie Algo
  if (r_sum * W_2(i) < r_bind)
		unbound_list = find(X(1,:,i) == 0);
		unbound_index = unbound_list(randi(length(unbound_list)));
		X(1,unbound_index,i) = 1;
		X(2,unbound_index,i) = L * rand(1);
		X(3,unbound_index,i) = L * rand(1);
	elseif (r_bind < r_sum * W_2(i) <= r_bind + r_unbind)
		bound_list = find(X(1,:,i) == 1);
		bound_index = bound_list(randi(length(bound_list)));
		X(1,bound_index,i) = 0;
		X(2:3,bound_index,i) = nan;
	else 
		unbound_list = find(X(1,:,i) == 0);
		unbound_index = unbound_list(randi(length(unbound_list)));
		
		bound_list = find(X(1,:,i) == 1);
		bound_index = bound_list(randi(length(bound_list)));

		X(1,unbound_index,i) = 1;

		%binding slightly differently to the particle
		theta = 2 * pi * rand(1);
%		X(2,unbound_index,i) = X(2,bound_index,i) + cos(theta) * 0.001;
%		X(3,unbound_index,i) = X(3,bound_index,i) + sin(theta) * 0.001;

		%biased binding
		X(2,unbound_index,i) = X(2,bound_index,i) + 0.001;
		X(3,unbound_index,i) = X(3,bound_index,i);

	end

	% Creates masks for bound & unbound values
	unbound = (X(1,:,i) == 0);
	bound = (X(1,:,i) == 0);
	X(2,unbound,i) = X(2,unbound,i-1);

	%random angle for 2d random walk
	theta = rand(1,length(X(4,bound,i))) * 2 * pi;
	random_steps_x = cos(theta) *  a_r * tau;
	random_steps_y = sin(theta) *  a_r * tau;
	%X(2,bound,i) = mod(X(2,bound,i-1) + random_steps_x,L);
	%X(3,bound,i) = mod(X(3,bound,i-1) + random_steps_y,L);

	X(4,bound,i) = X(4,bound,i) + random_steps_x;
	X(5,bound,i) = X(5,bound,i) + random_steps_y;

	X(2,bound,i) = mod(X(2,bound,i-1) + X(4,bound,i) * tau,L);
	X(3,bound,i) = mod(X(3,bound,i-1) + X(5,bound,i) * tau,L);

	% Updating proportion and time for A
	A(:, i) = [(sum(X(1, :,i)) / N_x), (A(2, i-1) + tau)];
end
for i = 2:(T+1)
	% set the position to the position of the previous time.
	Y(1,:,i) = Y(1,:,i-1);
	Y(2,:,i) = Y(2,:,i-1);
	Y(3,:,i) = Y(3,:,i-1);

	% Setting the rates at the start of each loop
	r_bind = k_bind*(N_y-sum(Y(1,:,i)));
	r_unbind = k_unbind*sum(Y(1,:,i));
	r_feedback = k_feedback*B(1,i-1)*(N_y - sum(Y(1,:,i)));
	r_sum = r_bind + r_unbind + r_feedback;

	% Generating a random number for use in Gillespie Algo

	tau = (1 / r_sum) * log(1 / W_3(i));

	% Picking the reaction using Gillespie Algo
  if (r_sum * W_4(i) < r_bind)
		unbound_list = find(Y(1,:,i) == 0);
		unbound_index = unbound_list(randi(length(unbound_list)));
		Y(1,unbound_index,i) = 1;
		Y(2,unbound_index,i) = L * rand(1);
		Y(3,unbound_index,i) = L * rand(1);
	elseif (r_bind < r_sum * W_4(i) <= r_bind + r_unbind)
		bound_list = find(Y(1,:,i) == 1);
		bound_index = bound_list(randi(length(bound_list)));
		Y(1,bound_index,i) = 0;
		Y(2:3,bound_index,i) = nan;
	else 
		unbound_list = find(Y(1,:,i) == 0);
		unbound_index = unbound_list(randi(length(unbound_list)));
		
		bound_list = find(Y(1,:,i) == 1);
		bound_index = bound_list(randi(length(bound_list)));

		Y(1,unbound_index,i) = 1;

		theta = 2 * pi * rand(1);
		Y(2,unbound_index,i) = Y(2,bound_index,i) + cos(theta) * 0.001;
		Y(3,unbound_index,i) = Y(3,bound_index,i) + sin(theta) * 0.001;
	end

	% Creates masks for bound & unbound values
	unbound = (Y(1,:,i) == 0);
	bound = (Y(1,:,i) == 0);
	Y(2,unbound,i) = Y(2,unbound,i-1);

	%random angle for 2d random walk
	theta = rand(1) * 2 * pi;
	random_steps_x = cos(theta) * ((2 * randi(2, 1, length(Y(2,bound,i)))) - 3) * v_r * tau;
	random_steps_y = sin(theta) * ((2 * randi(2, 1, length(Y(2,bound,i)))) - 3) * v_r * tau;
	Y(2,bound,i) = mod(Y(2,bound,i-1) + random_steps_x, L);
	Y(3,bound,i) = mod(Y(3,bound,i-1) + random_steps_y, L);

	% Updating proportion and time for A
	B(:, i) = [(sum(Y(1, :,i)) / N_y), (B(2, i-1) + tau)];
end

% The steady state  f the concentration for the ODE version of this is:
a_steady_state = 1 - (k_unbind/k_feedback);


%Setting up the ODE to plot
[t, a] =  ode45(@(t, a) k_bind * (1 - a) + k_feedback * (1 - a) * a  - k_unbind * a, [0, A(2,T)], A(1,1));

%subplot(1,3,1);
%ODE plot
%plot(t,a); hold on;

%Stochastic Plot
%plot(A(2,:),A(1,:));

%Steady state comparison
%yline(a_steady_state);

%Title and labels
%title('Simple Binding-Unbinding Model w/ ODE Comparison');
%xlabel('Time');
%ylabel('Bound Proportion of Chemical A');

%Plotting a histogram of particle location
%subplot(1,3,2);
%hist(X(2,:,T),20);

%title('Histogram of Particle Position at the Ending Time');
%xlabel('Position of Particles');
%ylabel('Number of Particles');

%subplot(1,3,3);

ylabel('Y Position');
xlabel = ('X Position');

for t = 1:20:T
  title(['Particle Position: ',num2str(t)]);
	scatter(X(2,:,t),X(3,:,t),250,"m",".");
	set(gca, 'Color', 'k');
	hold on;
	scatter(Y(2,:,t),Y(3,:,t),250,"c",".");
	hold off;
  xlim([0,L]);
  ylim([0,L]);
	pause(1/5000);
end

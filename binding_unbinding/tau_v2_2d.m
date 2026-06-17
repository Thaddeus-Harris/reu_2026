% *** THIS PROGRAM DEALS W/ ONE CHEMICAL SPECIES WITH A BINDING AND UNBINDING REACTION, WITH NO FEEDBACK ***

% rates: r_bind = k_bind * (1-a)
% 			 r_unbind = k_unbind * a
x_bind_num = 0;
x_unbind_num = 0;
x_fb_num = 0;
y_bind_num = 0;
y_unbind_num = 0;
y_fb_num = 0;
tic
N_x = 2000; % The number of molecules
N_x_starting_bound = N_x/10; % The number of molecules
N_y = 2000; % The number of molecules
N_y_starting_bound = N_y/10; % The number of molecules
T = 10000; % Number of reactions
k_bind = 4e-4; % Binding rate constant
k_unbind = 0.9; % Unbinding rate constant
k_feedback = 1.0; % Feedback rate constant
L = 1; % Length of cell membrane
v_r = 0.5; % Velocity for random walk
p = 0.5; % Probability of moving to the right in random walk (0.5 for no bias)

X = zeros(3,N_x,T);
Y = zeros(3,N_y,T);
% X(1,:,t) is a bool representing whether each particle is bound or not.
X(1,1:N_x_starting_bound,1) = ones(1,N_x_starting_bound,1);

Y(1,1:N_y_starting_bound,1) = ones(1,N_y_starting_bound,1);
% X(2,:,t) is a float representing the x position of each particle.
% X(3,:,t) is a float representing the y position of each particle.
X(2:3,1:N_x_starting_bound,1) = (L * rand(2,N_x_starting_bound));
X(2:3,N_x_starting_bound+1:N_x,1) = NaN(2,N_x-N_x_starting_bound);

Y(2:3,1:N_y_starting_bound,1) = (L * rand(2,N_y_starting_bound));
Y(2:3,N_y_starting_bound+1:N_y,1) = NaN(2,N_y-N_y_starting_bound);

% A(1, :) is a float representing the proportion of chemical X that is bound.
% A(2, :) is the time that this value is recorded.
% B(1, :) is a float representing the proportion of chemical X that is bound.
% B(2, :) is the time that this value is recorded.
A = zeros(2,T);
B = zeros(2,T);
W_1 = rand(1,T+1);
W_2 = rand(1,T+1);
toc
% This odd index is due to using i-1 to gather the value of the last index in each of the arrays, and I didn't want to add a switch statement every time for i = 1. This is more computationally efficient
tic
for i = 2:(T+1)
	% set the position to the position of the previous time.
	X(1,:,i) = X(1,:,i-1);
	X(2,:,i) = X(2,:,i-1);
	X(3,:,i) = X(3,:,i-1);
	 
	Y(1,:,i) = Y(1,:,i-1);
	Y(2,:,i) = Y(2,:,i-1);
	Y(3,:,i) = Y(3,:,i-1);
	% Setting the rates at the start of each loop
	r_bind_x = k_bind*(N_x-sum(X(1,:,i)));
	r_unbind_x = k_unbind*sum(X(1,:,i));
	r_feedback_x = k_feedback*A(1,i-1)*(N_x - sum(X(1,:,i)));

	r_bind_y = k_bind*(N_x-sum(Y(1,:,i)));
	r_unbind_y = k_unbind*sum(Y(1,:,i));
	r_feedback_y = k_feedback*B(1,i-1)*(N_y - sum(Y(1,:,i)));

	rates = [r_bind_x, r_unbind_x, r_feedback_x, r_bind_y, r_unbind_y, r_feedback_y]; 
	r_sum = sum(rates);

	% Generating a random number for use in Gillespie Algo


	tau = (1 / r_sum) * log(1 / W_1(i));

	% Picking the reaction using Gillespie Algo
  if (r_sum * W_2(i) < rates(1) )
		unbound_list = find(X(1,:,i) == 0);
		unbound_index = unbound_list(randi(length(unbound_list)));
		X(1,unbound_index,i) = 1;
		X(2,unbound_index,i) = L * rand(1);
		X(3,unbound_index,i) = L * rand(1);
		x_bind_num = x_bind_num + 1;
	elseif (r_sum * W_2(i) < sum(rates(1:2)))
		bound_list = find(X(1,:,i) == 1);
		bound_index = bound_list(randi(length(bound_list)));
		X(1,bound_index,i) = 0;
		X(2,bound_index,i) = nan;
		X(3,bound_index,i) = nan;
		x_unbind_num = x_unbind_num + 1;
	elseif (r_sum * W_2(i) < sum(rates(1:3)))
		unbound_list = find(X(1,:,i) == 0);
		unbound_index = unbound_list(randi(length(unbound_list)));
		
		bound_list = find(X(1,:,i) == 1);
		bound_index = bound_list(randi(length(bound_list)));

		X(1,unbound_index,i) = 1;

		theta = 2 * pi * W_2(i);
		X(2,unbound_index,i) = X(2,bound_index,i) + cos(theta) * 0.001;
		X(3,unbound_index,i) = X(3,bound_index,i) + sin(theta) * 0.001;

		x_fb_num= x_fb_num + 1;
	elseif (r_sum * W_2(i) < sum(rates(1:4)))
		unbound_list = find(Y(1,:,i) == 0);
		unbound_index = unbound_list(randi(length(unbound_list)));
		Y(1,unbound_index,i) = 1;
		Y(2,unbound_index,i) = L * rand(1);
		Y(3,unbound_index,i) = L * rand(1);

		y_bind_num= y_bind_num + 1;
	elseif (r_sum * W_2(i) < sum(rates(1:5)))
		bound_list = find(Y(1,:,i) == 1);
		bound_index = bound_list(randi(length(bound_list)));
		Y(1,bound_index,i) = 0;
		Y(2,bound_index,i) = nan;
		Y(3,bound_index,i) = nan;
		y_unbind_num = y_unbind_num + 1;
	else 
		unbound_list = find(Y(1,:,i) == 0);
		unbound_index = unbound_list(randi(length(unbound_list)));
		
		bound_list = find(Y(1,:,i) == 1);
		bound_index = bound_list(randi(length(bound_list)));

		Y(1,unbound_index,i) = 1;

		theta = 2 * pi * W_2(i);
		Y(2,unbound_index,i) = Y(2,bound_index,i) + cos(theta) * 0.001;
		Y(3,unbound_index,i) = Y(3,bound_index,i) + sin(theta) * 0.001;
		y_fb_num= y_fb_num + 1;
	end

	% Creates masks for bound & unbound values
	bound_x = (X(1,:,i) == 1);

	%random angle for 2d random walk
	theta_1 = rand(1,length(X(2,bound_x,i))) * 2 * pi;
	random_steps_x_1 = cos(theta_1)  * v_r * tau;
	random_steps_y_1 = sin(theta_1)  * v_r * tau;
	X(2,bound_x,i) = mod(X(2,bound_x,i) + random_steps_x_1, L);
	X(3,bound_x,i) = mod(X(3,bound_x,i) + random_steps_y_1, L);

	% Creates masks for bound & unbound values
	bound_y = (Y(1,:,i) == 1);

	theta_2 = rand(1,length(Y(2,bound_y,i))) * 2 * pi;
	random_steps_x_2 = cos(theta_2) * v_r * tau;
	random_steps_y_2 = sin(theta_2) * v_r * tau;

	Y(2,bound_y,i) = mod(Y(2,bound_y,i) + random_steps_x_2, L);
	Y(3,bound_y,i) = mod(Y(3,bound_y,i) + random_steps_y_2, L);
	
	% Updating proportion and time for A
	A(:, i) = [(sum(X(1,:,i)) / N_x), (A(2, i-1) + tau)];
	B(:, i) = [(sum(Y(1,:,i)) / N_y), (B(2, i-1) + tau)];
end
toc
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

x_bind_num 
x_unbind_num
x_fb_num
y_bind_num
y_unbind_num 
y_fb_num

%subplot(1,3,3);
ylabel = ('Y Position');
xlabel = ('X Position');

for t = 1:20:T
  title(['Particle Position: ',num2str(t)]);
  set(gca,'Color','k');
	box on; grid on; axis equal;
	scatter(X(2,:,t),X(3,:,t),1000,"m",".");
	hold on;
	scatter(Y(2,:,t),Y(3,:,t),1000,"c",".");
  set(gca,'Color','k');
	box on; grid on; axis equal;
	hold off;
  xlim([0,L]);
  ylim([0,L]);
	pause(1/5000);
end

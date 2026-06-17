% *** THIS PROGRAM DEALS W/ ONE CHEMICAL SPECIES WITH A BINDING AND UNBINDING REACTION, WITH NO FEEDBACK ***

% rates: r_bind = k_bind * (1-a)
% 			 r_unbind = k_unbind * a

N_x = 1500; % The number of molecules
N_x_starting_bound = N_x/10; % The number of molecules
N_y = 1500; % The number of molecules
N_y_starting_bound = N_y/10; % The number of molecules
T = 20000; % Number of reactions
k_bind = 0.1; % Binding rate constant
k_unbind = 0.9; % Unbinding rate constant
k_feedback = 1.0; % Feedback rate constant
L = 0.1; % Length of cell membrane
D = 0.02; % Velocity for random walk
v_r = 0.02; % Velocity for PK Process
lambda = 1;
r_particle = L/50;
% X(1,:,t) is a bool representing whether each particle is bound or not.
X = zeros(5,N_x,T);
Y = zeros(5,N_y,T);
X(1,1:N_x_starting_bound,1) = 1;
Y(1,:,:) = zeros(1,N_y,T);
Y(1,1:(N_y_starting_bound),1) = 1;
Y(1,1:N_y_starting_bound,1) = 1;
% X(2,:,t) is a float representing the x position of each particle.
% X(3,:,t) is a float representing the y position of each particle.
X(2:3,1:N_x_starting_bound,1) = (L * rand(2,N_x_starting_bound));
%X(2:3,1:N_x_starting_bound,1) = L/2;
X(2:3,N_x_starting_bound+1:N_x,1) = nan;
Y(2:3,1:N_y_starting_bound,1) = (L * rand(2,N_y_starting_bound));
%Y(2:3,1:N_y_starting_bound,1) = L/2;
Y(2:3,N_y_starting_bound+1:N_y,1) = nan;

X(2:3,:,2) = X(2:3,:,1);
Y(2:3,:,2) = Y(2:3,:,1);

X(4,:,1) = cos(rand(1,N_x) * 2 * pi) * v_r;
X(5,:,1) = sin(rand(1,N_x) * 2 * pi) * v_r;

Y(4,:,1) = cos(rand(1,N_y) * 2 * pi) * v_r;
Y(5,:,1) = sin(rand(1,N_y) * 2 * pi) * v_r;
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
W_5 = rand(N_x,T+1);
W_6 = rand(N_y,T+1);

profile on;
% This odd index is due to using i-1 to gather the value of the last index in each of the arrays, and I didn't want to add a switch statement every time for i = 1. This is more computationally efficient
for i = 2:T
	% set the position to the position of the previous time.
	X(1,:,i) = X(1,:,i-1);
	X(4,:,i) = X(4,:,i-1);
	X(5,:,i) = X(5,:,i-1);

	Y(1,:,i) = Y(1,:,i-1);
	Y(4,:,i) = Y(4,:,i-1);
	Y(5,:,i) = Y(5,:,i-1);
	% Setting the rates at the start of each loop
	r_bind_x = k_bind*(N_x-sum(X(1,:,i)));
	r_unbind_x = k_unbind*sum(X(1,:,i));
	r_feedback_x = k_feedback*A(1,i-1)*(N_x - sum(X(1,:,i)));
	r_bind_y = k_bind*(N_y-sum(Y(1,:,i)));
	r_unbind_y = k_unbind*sum(Y(1,:,i));
	r_feedback_y = k_feedback*B(1,i-1)*(N_y - sum(Y(1,:,i)));
	rates = [r_bind_x, r_unbind_x,r_feedback_x, r_bind_y,r_unbind_y,r_feedback_y];
	r_sum = sum(rates);

	% Generating a random number for use in Gillespie Algo

	tau = (1 / r_sum) * log(1 / W_1(i));

	% Picking the reaction using Gillespie Algo
  if (r_sum * W_2(i) < rates(1))
		unbound_list = find(X(1,:,i) == 0);
		unbound_index = unbound_list(randi(length(unbound_list)));
		X(1,unbound_index,i) = 1;
		X(2,unbound_index,i) = L * rand(1);
		X(3,unbound_index,i) = L * rand(1);
	elseif (r_sum * W_2(i) < sum(rates(1:2)))
		bound_list = find(X(1,:,i) == 1);
		bound_index = bound_list(randi(length(bound_list)));
		X(1,bound_index,i) = 0;
		X(2:3,bound_index,i) = nan;
	elseif (r_sum * W_2(i) < sum(rates(1:3)))
		unbound_list = find(X(1,:,i) == 0);
		unbound_index = unbound_list(randi(length(unbound_list)));
		
		bound_list = find(X(1,:,i) == 1);
		bound_index = bound_list(randi(length(bound_list)));

		X(1,unbound_index,i) = 1;

		theta = 2 * pi * rand(1);
		X(2,unbound_index,i) = X(2,bound_index,i) + ((randi([0,1]) * 2) - 1) * r_particle;
		X(3,unbound_index,i) = X(3,bound_index,i) + ((randi([0,1]) * 2) - 1) * r_particle;

		X(4,unbound_index,i) = X(4,bound_index,i);
		X(5,unbound_index,i) = X(5,bound_index,i);
	elseif (r_sum * W_2(i) < sum(rates(1:4)))
		unbound_list = find(Y(1,:,i) == 0);
		unbound_index = unbound_list(randi(length(unbound_list)));
		Y(1,unbound_index,i) = 1;
		Y(2,unbound_index,i) = L * rand(1);
		Y(3,unbound_index,i) = L * rand(1);
	elseif (r_sum * W_2(i) < sum(rates(1:5)))
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
		Y(2,unbound_index,i) = Y(2,bound_index,i) + ((randi([0,1]) * 2) - 1) * r_particle;
		Y(3,unbound_index,i) = Y(3,bound_index,i) + ((randi([0,1]) * 2) - 1) * r_particle;

		Y(4,unbound_index,i) = Y(4,bound_index,i);
		Y(5,unbound_index,i) = Y(5,bound_index,i);
	end

	flip_u = W_5(:,i) < (lambda * tau);
	flip_v = W_6(:,i) < (lambda * tau);
	%random angle for 2d random walk
	theta_RW_u = rand(1,N_x) * 2 * pi;
	theta_PK_u = rand(1,length(flip_u)) * 2 * pi;
	theta_RW_v = rand(1,N_y) * 2 * pi;
	theta_PK_v = rand(1,length(flip_v)) * 2 * pi;
	random_steps_u_x = cos(theta_RW_u) *  D * tau;
	random_steps_u_y = sin(theta_RW_u) *  D * tau;
	random_steps_v_x = cos(theta_RW_v) *  D * tau;
	random_steps_v_y = sin(theta_RW_v) *  D * tau;
	random_velocity_angle_u_x = cos(theta_PK_u) *  v_r;
	random_velocity_angle_u_y = sin(theta_PK_u) *  v_r;
	random_velocity_angle_v_x = cos(theta_PK_v) *  v_r;
	random_velocity_angle_v_y = sin(theta_PK_v) *  v_r;

% Adjusting the velocity

	X(4,flip_u,i) = random_velocity_angle_u_x(flip_u);
	X(5,flip_u,i) = random_velocity_angle_u_y(flip_u);

	Y(4,flip_v,i) = random_velocity_angle_v_x(flip_v);
	Y(5,flip_v,i) = random_velocity_angle_v_y(flip_v);

	%Adjustments to posn from velocity
	X(2,:,i) = mod(X(2,:,i) + X(4,:,i) * tau,L);
	X(3,:,i) = mod(X(3,:,i) + X(5,:,i) * tau,L);

	Y(2,:,i) = mod(Y(2,:,i) + Y(4,:,i) * tau,L);
	Y(3,:,i) = mod(Y(3,:,i) + Y(5,:,i) * tau,L);


	%Adjustments from random walk
	X(2,:,i) = mod(X(2,:,i) + random_steps_u_x, L);
	X(3,:,i) = mod(X(3,:,i) + random_steps_u_y, L);

	Y(2,:,i) = mod(Y(2,:,i) + random_steps_v_x, L);
	Y(3,:,i) = mod(Y(3,:,i) + random_steps_v_y, L);

	% Updating Values
	X(2,:,i+1) = X(2,:,i); 
	X(3,:,i+1) = X(3,:,i);

	Y(2,:,i+1) = Y(2,:,i); 
	Y(3,:,i+1) = Y(3,:,i); 
	x_bound = find(X(1,:,i) == 1);   % global indices of bound X particles
	y_bound = find(Y(1,:,i) == 1);   % global indices of bound Y particles

	xp = X(2:3, x_bound, i);         % 2 x nXbound positions
	yp = Y(2:3, y_bound, i);         % 2 x nYbound positions

	% --- cross-species: X against Y ---
	dx = xp(1,:).' - yp(1,:);
	dy = xp(2,:).' - yp(2,:);
	distXY = sqrt(dx.^2 + dy.^2);
	[lxi_xy, lyi_xy] = find(distXY < r_particle);

	% --- same-species: X against X ---
	dxx = xp(1,:).' - xp(1,:);
	dyy = xp(2,:).' - xp(2,:);
	distXX = sqrt(dxx.^2 + dyy.^2);
	distXX(1:numel(x_bound)+1:end) = Inf;   % exclude each particle vs itself
	[lxa, lxb] = find(distXX < r_particle);

	% --- same-species: Y against Y ---
	dvx = yp(1,:).' - yp(1,:);
	dvy = yp(2,:).' - yp(2,:);
	distYY = sqrt(dvx.^2 + dvy.^2);
	distYY(1:numel(y_bound)+1:end) = Inf;   % exclude each particle vs itself
	[lya, lyb] = find(distYY < r_particle);

	% --- map local (bound-list) indices back to global indices, then revert ---
	X_collide = x_bound(unique([lxi_xy; lxa; lxb]));
	Y_collide = y_bound(unique([lyi_xy; lya; lyb]));

	X(2:3, X_collide, i+1) = X(2:3, X_collide, i-1);
	Y(2:3, Y_collide, i+1) = Y(2:3, Y_collide, i-1);
	% Updating proportion and time for A

	A(:, i) = [(sum(X(1, :,i)) / N_x), (A(2, i-1) + tau)];
	B(:, i) = [(sum(Y(1, :,i)) / N_y), (B(2, i-1) + tau)];
	if mod(i,100) == 0
        fprintf('step %d / %d\n', i, T+1);
    end

end
profile off;
profshow(profile('info'), 20);

% The steady state  f the concentration for the ODE version of this is:
a_steady_state = 1 - (k_unbind/k_feedback);


%Setting up the ODE to plot
%[t, a] =  ode45(@(t, a) k_bind * (1 - a) + k_feedback * (1 - a) * a  - k_unbind * a, [0, A(2,T)], A(1,1));

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
	scatter(X(2,:,t),X(3,:,t),500,"m",".");
	set(gca,'Color','k');
	hold on;
	scatter(Y(2,:,t),Y(3,:,t),500,"c",".");
	hold off;
  xlim([0,L]);
  ylim([0,L]);
	pause(1/5000);
end

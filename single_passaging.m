clc; clear;

% PARAMETERS
number_of_species = 10^5; 
birth_rate = 0.026;

% PASSAGING PARAMETERS
passaging_interval = 4;      % every 4 days

percent_of_pass = 10;        % keep 10% of cells
T_end = 4*24;                % total hours
max_passages = ceil((T_end/24)/passaging_interval);
X_old = ones(1, number_of_species);           % initial population array
P_n_t_Yule = zeros(1, number_of_species);
P_sur_sim = zeros(1, max_passages);
passage_times = zeros(1, max_passages);
mean_clone_size = zeros(1, max_passages);
clone_size_history = cell(1, max_passages);   % cell array to store distributions
flag = 0;                                     % counts number of passages
t = 0;
epsilon = 0.0003;
alpha = 0.01;
% --------consistent daily storage ------------
tvals = 0:4;
time_store = zeros(size(tvals));
growing_total_pop_store = zeros(size(tvals));
next_day = 0;
% ---------------------------------------------
% previous_pass_time=0;
while t < T_end
    t/24
    %tau-leaping growth
    RxPlus  = birth_rate.* X_old;   %division rate of clone
    positive_idx = find(X_old > 0);
    if isempty(positive_idx)
        tau = 0;  % no more growth possible
    else
        m = max(epsilon * X_old(positive_idx), 1);
        tau = min(m ./ RxPlus(positive_idx)) * alpha;
    end
    
    t = t + tau;
    births = zeros(size(RxPlus));                     % initialize
    births(positive_idx) = poissrnd(RxPlus(positive_idx) .* tau);  % number of times an event occures independently in a small interval of time, choosing a random number with mean (RxPlus*tau)
    X_old = X_old+births;
    q = percent_of_pass / 100; % proportion of sampling
    
        next_pass_time = passaging_interval*24*(1+flag);
        n0=1;
        if(t>=(passaging_interval*24)-epsilon)
            clone_dist_before_1stp=X_old; 
            n_vals=1:max(clone_dist_before_1stp);
            P_n_t_Yule =  exp(-birth_rate*n0*t).* (1 - exp(-birth_rate*t)).^(n_vals - n0);
        end
        if t >= next_pass_time
            N_total = sum(X_old);
            n_pass = floor(N_total * q);
            P_survive_theory = 1-(1-q).^X_old; % for clone vector
            P_extinct = (1-q).^X_old;
            
            X_new = binornd(X_old, q);
            flag = flag + 1;
            
            clone_size_history{flag} = X_new;   % store clone sizes
            mean_clone_size(flag)=mean(clone_size_history{flag});
            passage_times(flag) = t/24;
            num_survived = sum(X_new > 0);
            P_sur_sim(flag) = num_survived / number_of_species;
            X_old = X_new;
        end  % if t >= next_pass_time
    % -------- FIXED STORAGE (AFTER growth + passage) --------
    current_day = floor(t/24);

    while next_day <= 4 && current_day >= next_day
        time_store(next_day + 1) = next_day;
        growing_total_pop_store(next_day + 1) = sum(X_old);
        next_day = next_day + 1;
    end
    % ------------------------------------------------------
    
end %while t < T_end  


save("All_10^5_q_0.1_t_4_single_passage.mat","growing_total_pop_store","time_store","P_sur_sim","passage_times","clone_size_history","clone_dist_before_1stp","P_n_t_Yule")
figure (1);
histogram(clone_dist_before_1stp, "Normalization", "pdf")
xlabel('Clone size')
ylabel('PDF')
hold on
plot(n_vals,P_n_t_Yule)
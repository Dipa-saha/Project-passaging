clc; clear;

% PARAMETERS
number_of_species = 10^5; 
birth_rate = 0.026;

% PASSAGING PARAMETERS
passaging_interval = 4;      % every 4 days, t

percent_of_pass = 10;        % keep 10% of cells
T_end = 200*24;                % total hours
max_passages = ceil((T_end/24)/passaging_interval);
X_old = ones(1, number_of_species);           % initial population array
P_sur_sim = zeros(1, max_passages);
passage_times = zeros(1, max_passages);
mean_clone_size = zeros(1, max_passages);
clone_size_history = cell(1, max_passages);   % cell array to store distributions
flag = 0;                                     % counts number of passages
t = 0;
epsilon = 0.0003;
alpha = 0.01;
% -----------consistent daily storage ---------
tvals = 0:200;
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
        P_H2 = 1 - ((1-q)/(1 + q*(exp(birth_rate*passaging_interval*24)-1)))^n0;  %for single passage
        if t >= next_pass_time
           
            N_total = sum(X_old);
            n_pass = floor(N_total * q);
         
            % Random sampling pool hypergeometric
             pool = repelem(1:length(X_old), X_old);   %repeats each clone ID according to how many cells it has
             picked = randsample(pool, n_pass, false); % randomly picking n_pass without replacement
             X_new = histcounts(picked,1:length(X_old)+1); % counting how many cells picked from each clone
            
            flag = flag + 1;
            clone_size_history{flag} = X_new;   % store clone sizes
            mean_clone_size=mean(clone_size_history{flag});
            
            passage_times(flag) = t/24;
            num_survived = sum(X_new > 0);
            P_sur_sim(flag) = num_survived / number_of_species;
                  
            X_old = X_new;
        end  % if t >= next_pass_time
    % -------- FIXED STORAGE (AFTER growth + passage) ------
    current_day = floor(t/24);

    while next_day <= 200 && current_day >= next_day
        time_store(next_day + 1) = next_day;
        growing_total_pop_store(next_day + 1) = sum(X_old);
        next_day = next_day + 1;
    end
    % ------------------------------------------------------
    
end %while t < T_end 
%Theory for mean clone size
mean_x = zeros(1,10);
for n=1:50
x0 = 1;
mean_x(n) = x0*q^n * exp(birth_rate*passaging_interval*24*n);
end
%------
save("All_10^5_q_0.1_t_4_50_passages_HG.mat","growing_total_pop_store","time_store","P_sur_sim","passage_times","clone_size_history","mean_clone_size","mean_x")
figure (1);
plot(time_store, growing_total_pop_store, 's--', 'LineWidth',1.5)
figure(2);
semilogy(passage_times,mean_clone_size, 'o', 'LineWidth',1.5)
xlabel('t (days)')
ylabel('\langle x \rangle')
hold on
semilogy(passage_times,mean_x, '-', 'LineWidth',1.5)
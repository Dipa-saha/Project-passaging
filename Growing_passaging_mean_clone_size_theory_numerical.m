clc; clear;

% PARAMETERS
number_of_species = 10^5; 
birth_rate = 0.026;

% PASSAGING PARAMETERS
passaging_interval = 4;      % every 4 days
percent_of_pass = 10;        % keep 10% of cells

T_end = 40*24;               
max_passages = ceil((T_end/24 )/passaging_interval);

X_old = ones(1, number_of_species);

P_sur_sim = zeros(1, max_passages);
P_sur_theory = zeros(1, max_passages);
passage_times = zeros(1, max_passages);
clone_size_history = cell(1, max_passages);

flag = 0;
t = 0;

epsilon = 0.0003;
alpha = 0.01;
% -------- NEW: consistent daily storage --------
tvals = 0:40;
time_store = zeros(size(tvals));
growing_total_pop_store = zeros(size(tvals));
next_day = 0;
% ---------------------------------------------

while t < T_end
t/24
    % tau-leaping growth
    RxPlus  = birth_rate .* X_old;
    positive_idx = find(X_old > 0);

    if isempty(positive_idx)
        tau = 0;
    else
        m = max(epsilon * X_old(positive_idx), 1);
        tau = min(m ./ RxPlus(positive_idx)) * alpha;
    end
    
    t = t + tau;

    births = zeros(size(RxPlus));
    births(positive_idx) = poissrnd(RxPlus(positive_idx) .* tau);
    X_old = X_old + births;

    % PASSAGING
    q = percent_of_pass / 100;
    next_pass_time = 4*24 + passaging_interval*24*(flag);

    if t >= next_pass_time
        N_total = sum(X_old);
        n_pass = floor(N_total * q);

        P_survive_theory = 1 - (1 - q).^X_old;
        % Random sampling binomial
        X_new = binornd(X_old, q);

        flag = flag + 1;
        clone_size_history{flag} = X_new;
        passage_times(flag) = t/24;

        num_survived = sum(X_new > 0);
        P_sur_sim(flag) = num_survived / number_of_species;
        P_sur_theory(flag) = mean(P_survive_theory);
        X_old = X_new;   % AFTER PASSAGING
    end

    % -------- FIXED STORAGE (AFTER growth + passage) --------
    current_day = floor(t/24);

    while next_day <= 40 && current_day >= next_day
        time_store(next_day + 1) = next_day;
        growing_total_pop_store(next_day + 1) = sum(X_old);
        next_day = next_day + 1;
    end
    % ------------------------------------------------------

end
theory_vals = zeros(size(time_store));
for i = 1:length(time_store)
    t = time_store(i);
    x0 = growing_total_pop_store(i);
    %if use_downscaling
        theory_vals(i) = x0 * exp(birth_rate * t) * (0.975 ^ floor(t));    
end
save("All_10^5_N_t_q_0.1_40_binom.mat","theory_vals","time_store","growing_total_pop_store")
figure (1);
plot(time_store, growing_total_pop_store, 'o-', 'LineWidth',1.5)
hold on
plot(time_store, theory_vals, 'k-', 'LineWidth', 2)
xlabel('t(days)')
ylabel('N(t)')
legend('Simulation', 'Theory')

%Theory & simulation results for mean clone size----
mean_x = zeros(1,10);
mean_clone_size = zeros(1,10);
for n=1:10
mean_clone_size(n) = mean(clone_size_history{n});
x0 = 1;
mean_x(n) = x0*q^n * exp(birth_rate*passaging_interval*24*n);
end
figure(2);
semilogy(passage_times,mean_clone_size, 'o','MarkerSize',10)
xlabel('t (days)')
ylabel('\langle x \rangle')
hold on
semilogy(passage_times,mean_x, 'k-', 'LineWidth',1.5)
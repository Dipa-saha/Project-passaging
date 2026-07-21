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
P_n_t_Yule_before = zeros(1, number_of_species);
P_n_t_Yule_after = zeros(1, number_of_species);
clone_size_history = cell(1, max_passages);   % cell array to store distributions
flag = 0;                                     % counts number of passages
t = 0;
epsilon = 0.0003;
alpha = 0.01;

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
        if flag==0 && t >= passaging_interval*24
        t_before = t;
        clone_dist_before_1stp=X_old;
        surv_before = clone_dist_before_1stp(clone_dist_before_1stp>0);
        [counts_before, bins_before] = histcounts(surv_before,...
       'BinMethod','integers','Normalization','probability');
        r = exp(-birth_rate * t_before);
        n_before = bins_before(1:end-1);
        P_before = counts_before;
        x_before = n_before * abs(log(1-r));
        n_vals=1:max(clone_dist_before_1stp);
        P_n_t_Yule_before =  r.* (1 - r).^(n_vals - 1);
        x_theory_before = n_vals * abs(log(1-r));
        end
        if t >= next_pass_time
            t_after = t;
            N_total = sum(X_old);
            n_pass = floor(N_total * q);
            P_survive_theory = 1-(1-q).^X_old; % for clone vector
            P_extinct = (1-q).^X_old;
           
            X_new = binornd(X_old, q);
            flag = flag + 1;
            
            clone_size_history{flag} = X_new;   % store clone sizes
            X_old = X_new;
            surv_after  = clone_size_history{flag}(clone_size_history{flag}>0);
            [counts_after, bins_after] = histcounts(surv_after,...
                'BinMethod','integers','Normalization','probability');
            n_after = bins_after(1:end-1) + diff(bins_after)/2;
            P_after = counts_after;
            r_eff = r/(r + q - q*r);
            x_after = n_after * abs(log(1-r_eff));
            n_vals_=1:max(clone_size_history{flag});
            P_n_t_Yule_after =  r_eff.* (1 - r_eff).^(n_vals_ - 1);
            x_theory_after = n_vals_ * abs(log(1-r_eff));
        end  % if t >= next_pass_time
end %while t < T_end  

y_before_scaled = P_before * (1-r)/r;
y_after_scaled  = P_after  * (1-r_eff)/r_eff;

y_theory_before_scaled = P_n_t_Yule_before * (1-r)/r;
y_theory_after_scaled  = P_n_t_Yule_after  * (1-r_eff)/r_eff;
save("All_10^5_q_0.1_t_4_single_passage_similar_shape of dist.mat","x_before","y_before_scaled","x_after","y_after_scaled","x_theory_before","y_theory_before_scaled","x_theory_after","y_theory_after_scaled")
figure(4)
semilogy(x_before(1:4:end),...
       y_before_scaled(1:4:end),'bo')
hold on
semilogy(x_theory_before,...
       y_theory_before_scaled,'b--','LineWidth',2)
semilogy(x_after,...
         y_after_scaled,'rs')
semilogy(x_theory_after,...
       y_theory_after_scaled,'r--','LineWidth',2)
xlabel('Scaled Population Size $(x)$','Interpreter','latex')
ylabel('Scaled Probability $(P_x)$','Interpreter','latex')
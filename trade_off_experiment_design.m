clc; clear;

% PARAMETERS
number_of_species = 10^5; 
birth_rate = 0.026;

% PASSAGING PARAMETERS
passaging_interval = 4;      % every 4 days

percent_of_pass = 10;        % keep 10% of cells
T_end = 80*24;                % total hours
max_passages = ceil((T_end/24)/passaging_interval);
x0=[1,2,4,6,8,10,20,60,100];
P_sur_sim = zeros(1,length(x0));
P_sur_theory = zeros(1,length(x0));
P_s_theory = zeros(1,length(x0));
P_ext = zeros(1, length(x0));
ratio_theory = zeros(1,length(x0));
ratio = zeros(1,length(x0));
initial_cell_stock = zeros(1,length(x0));
K_20 = zeros(1,length(x0));
for i=1:length(x0)
X_old = ones(1, number_of_species)*x0(i);           % initial population array
initial_cell_stock(i) = sum(X_old);
t = 0;
epsilon = 0.0003;
alpha = 0.01;
flag=0;
q = percent_of_pass / 100; % proportion of sampling
while t < T_end
    [i,t/24]
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
    
    
        next_pass_time = passaging_interval*24*(1+flag);
        if t >= next_pass_time
            flag=flag+1;
            % Random sampling 
            X_new = binornd(X_old, q);
            X_old = X_new;
        end  % if t >= next_pass_time
    
end %while t < T_end 
num_survived = sum(X_new > 0);
P_sur_sim(i) = num_survived / number_of_species;
ratio(i)=num_survived/initial_cell_stock(i);
end
% Theory
for i=1:length(x0)
    r = exp(-birth_rate*passaging_interval*24);
    P_ext(i) = min(1,r*(1-q)/(q*(1-r)));
    P_sur_theory(i) = 1-(P_ext(i)).^x0(i);
    X_old = ones(1, number_of_species)*x0(i);           % initial population array
    initial_cell_stock(i) = sum(X_old);
    ratio_theory(i) = P_sur_theory(i)*(number_of_species/initial_cell_stock(i));
end
save("trade_off_10^5_all_new.mat","x0","P_sur_sim","ratio","ratio_theory","P_sur_theory")
%Plot
figure(4);
semilogx(x0,ratio,'o', 'MarkerSize', 8)
hold on
semilogx(x0,ratio_theory,'k-','LineWidth', 2)
xlabel('x_0')
ylabel('clone survives/initial cell stock')

figure(5);
semilogx(x0,P_sur_sim,'o','MarkerSize', 8)
hold on
semilogx(x0,P_sur_theory,'k-','LineWidth', 2)
xlabel('x_0')
ylabel('fraction of clone survives')
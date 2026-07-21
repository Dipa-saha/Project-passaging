clc; clear

q = 0.1;
number_of_species = 100000;

% --- Simulation ---
k0 = 5;
X_old_k_5 = ones(1, number_of_species)*k0;
X_new_k_5_binom = binornd(X_old_k_5, q);

n0 = 10000;
X_old_k_10k = ones(1, number_of_species)*n0;
X_new_k_10k_binom = binornd(X_old_k_10k, q);

% --- Analytical ---
x = 0:k0;
P_k_5 = binopdf(x, k0, q);

mu = n0*q; %mean
sigma = sqrt(n0*q*(1-q)); %variance
y = (mu - 4*sigma):(mu + 4*sigma);
P_k_10k = normpdf(y, mu, sigma);

% --- Plot ---
subplot(1,2,1)
hold on
xline(k0*q,'k--','LineWidth',2)
histogram(X_new_k_5_binom,'Normalization','pdf')
hold on
plot(x, P_k_5,'k','LineWidth',2)
xlabel('Clone sizes')
ylabel('PDF')
title('k = 5')

subplot(1,2,2)
hold on
xline(n0*q,'k--','LineWidth',2)
histogram(X_new_k_10k_binom,'Normalization','pdf')
hold on
plot(y, P_k_10k,'k','LineWidth',2)
xlabel('Clone sizes')
ylabel('PDF')
title('k = 10000')
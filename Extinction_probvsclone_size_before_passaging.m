clc; clear;

birth_rate = 0.026;   
q = 0.1;
number_of_species = 1000;

epsilon = 0.0003;
alpha = 0.01;
T_end= 4*24;

k0_vals = [2,4,6,8,10,15,20,30,40,50,60,70,80,90,100];
p_ext_vals = zeros(size(k0_vals));
for i = 1:length(k0_vals)
    
    % Reset time for each k0
    t = 0;
    k0 = k0_vals(i);
    X_old = ones(1, number_of_species) * k0;
   
    % binomial sampling 
    X_new = binornd(X_old, q);
    
    % extinction probability
    num_ext = sum(X_new==0);
    p_ext_vals(i) = (num_ext / number_of_species);
   
end
p_ext_theory_vals=(1-q).^k0_vals;

% Plot
figure
plot(k0_vals, p_ext_vals, 'o','LineWidth',2)
xlabel('k')
ylabel('P_{ext}(k)')
hold on
plot(k0_vals, p_ext_theory_vals, '-','LineWidth',2)
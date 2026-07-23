clc; clear;

number_of_species = 1000;
birth_rate = 0.026;
Numruns=1;
q_vals = 0.025:0.05:1;
    for run=1:Numruns
    t_vals = 3;
    P0 = zeros(Numruns, length(q_vals));
    epsilon = 0.0003;
    for j = 1:length(q_vals)
    r = exp(-birth_rate*t_vals*24);
    T_end = (60*t_vals)*24;
        t = 0;
        flag = 0;
        q = q_vals(j);
        X_old = ones(1, number_of_species);
        
        P_sur_sim = [];

        while t < T_end
            [run,q,t/24]
            RxPlus = birth_rate .* X_old;

            positive_idx = find(X_old > 0);

            if isempty(positive_idx)
                break;
            end

            m = max(epsilon * X_old(positive_idx), 1);
            tau = min(m ./ RxPlus(positive_idx)) * 0.01;

            t = t + tau;

            births = zeros(size(X_old));
            births(positive_idx) = poissrnd(RxPlus(positive_idx) .* tau);

            X_old = X_old + births;

            next_pass_time = t_vals*24*(1+flag);
            if t >= next_pass_time
                X_new = binornd(X_old, q);  %Binomial
                X_old = X_new;

                flag = flag + 1;

                P_sur_sim(flag) = sum(X_new > 0) / number_of_species;

                % -------------------------
                % CHECK CONVERGENCE
                % -------------------------
                if flag > 2
                    differ = P_sur_sim(flag-1)-P_sur_sim(flag);
                    b = abs(differ)/P_sur_sim(flag);
                    disp(['deltaP/P = ', num2str(b)])
                    if abs(differ)/P_sur_sim(flag) <= 0.001
                        break;
                    end
                end
            end
        end

        % store final plateau value
        if ~isempty(P_sur_sim)
            P0(run,j) = P_sur_sim(end);
        end
    
    
    end
    end
    %theory
    A = zeros(Numruns, length(q_vals));
    P0_theory = zeros(Numruns, length(q_vals));
for j = 1:length(q_vals)
    r = exp(-birth_rate*t_vals*24);
    q = q_vals(j);
    A(j) = min(1,r*(1-q)/(q*(1-r))); %extinction
    P0_theory(j) = 1-A(j);
end
save("q_t_Ps_saturation_t_3_L0_1k.mat",'q_vals','t_vals','P0','P0_theory')
%Plot
figure(1)
plot(q_vals, P0,'o')
hold on
plot(q_vals, P0_theory,'-')
xlabel('q')
ylabel('1-P0*')  %saturated survival prob

% dummy legend
blue = [0 0.45 0.74];
h1 = plot(nan, nan, '--s', ...
    'Color', blue, ...
    'MarkerFaceColor', blue, ...
    'MarkerEdgeColor', blue, ...
    'LineWidth', 2, ...
    'MarkerSize', 10);

hold on

red = [0.85 0.33 0.1];
h2 = plot(nan, nan, '--^', ...
    'Color', red, ...
    'MarkerFaceColor', red, ...
    'MarkerEdgeColor', red, ...
    'LineWidth', 2, ...
    'MarkerSize', 10);

green = [0.47 0.67 0.19];
h3 = plot(nan, nan, '--o', ...
    'Color', green, ...
    'MarkerFaceColor', green, ...
    'MarkerEdgeColor', green, ...
    'LineWidth', 2, ...
    'MarkerSize', 10);

legend([h1 h2 h3], {'t = 1','t = 2','t = 3'})


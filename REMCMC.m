function output = REMCMC(data, phi_measure, tau)

% -----------------------------------------------------------------------
% This code implements the Replica Exchange Markov Chain Monte Carlo method
% for finding the MIB described in Kitozono et al (2018), "Efficient Algorithms
% for Searching the Minimum Information Partition in Integrated
% Information Theory." This algorithm differs from theirs in that it looks
% for a bipartiton that minimizes *normalized* integrated information,
% rather than non-normalized integrated information. In principle this
% algorithm can be used for large systems, but it's very slow to converge
% for any size system. For this reason, we've written this code such that
% if you terminate it while running, it will store its current best guesses
% for the MIB in the variable 'output' in your workspace
% -----------------------------------------------------------------------


% If you interrupt this function with ctrl+c, it will save the current
% partitions and integrated information values in the variable "output":
disp('If you terminate this script with ctrl+c, its interim results will be saved in a variable called "output"')

cleanupObj = onCleanup(@return_current_output);
global interrupt_output;

function return_current_output()
    assignin('base', 'output', interrupt_output);
end

N=size(data,1); % number of nodes
randn('seed',12345);


% INITIALIZE TEMPERATURES
for i = 1:6
    S_t{i} = randsample(N,randsample(N-1,1)); % random bipartition
    
    % change the community assignment of a random node
    e = randsample(N,1);
    if any(S_t{i}==e)
        S_c{i}=setdiff(S_t{i},e);
    else
        S_c{i}=[S_t{i}; e];
    end
    
    % compute normalized integrated information for the old bipartition
    S_c_part=ones(1,N);
    S_c_part(S_c{i})=2;
    S_c_out = phi_comp(data, S_c_part, phi_measure, tau,0);
    S_c_phi = S_c_out.phi_norm;
    
    % compute normalized integrated information for the new bipartition
    S_t_part=ones(1,N);
    S_t_part(S_t{i})=2;
    S_t_out = phi_comp(data, S_t_part, phi_measure, tau,0);
    S_t_phi = S_t_out.phi_norm;
    
    % magnitude of changes in integrated information as a result of
    % bipartition changes
    changes(i) = abs(S_t_phi-S_c_phi);
    
end

temps = 1:100000:1000000; % try out a large range of candidate temperatures

% fit acceptance ratio as a function of temperature
for i = 1:length(temps)
    rat(i) = mean(exp(-temps(i)*changes));
end

% set the highest temperature so that its average acceptance ratio is 0.01
% (parameter set in Kitazono et al)
top_temp_f=fit(temps',rat'-.01,'pchipinterp');
new_temps(1) = fzero(top_temp_f,[0 max(temps)]);

% set the lowest temperature so that its average acceptance ratio is 0.5
% (parameter set in Kitazono et al)
bottom_temp_f=fit(temps',rat'-.5,'pchipinterp');
new_temps(6) = fzero(bottom_temp_f,[0 max(temps)]);

% set intermediate temperatures through a geometric progression
clear temps
temps=new_temps;
for m = 2:5
    temps(m) = temps(1)*(temps(6)/temps(1))^((m-1)/5);
end
clear new_temps


N=size(data,1);

% start with random bipartitions for each of the 6 temperatures
for i = 1:length(temps)
    S_t{i} = randsample(N,randsample(N-1,1));
end


figure; % figure will monitor minimization progress
converge=0; % script will terminate when the convergence criterion flips to 1
t=1;
while converge ==0 
    
    for i = 1:length(temps)
        
        % change the community assignment of a random node
        e = randsample(N,1);
        
        if any(S_t{i}==e)
            S_c{i}=setdiff(S_t{i},e);
        else
            S_c{i}=[S_t{i}; e];
        end
        
        % new bipartition
        S_c_part=ones(1,N);
        S_c_part(S_c{i})=2;
        
        % old bipartition
        S_t_part=ones(1,N);
        S_t_part(S_t{i})=2;
        
        % normalized integrated information across the new bipartition
        S_c_out = phi_comp(data, S_c_part, phi_measure, tau,0);
        S_c_phi = S_c_out.phi_norm;
        S_c_phi_full = S_c_out.phi;
        
        % normalized integrated information across the old bipartition
        S_t_out = phi_comp(data, S_t_part, phi_measure, tau,0);
        S_t_phi = S_t_out.phi_norm;
        S_t_phi_full = S_t_out.phi;
        
        % changes in integrated information
        phi_diff(i,t) = S_t_phi-S_c_phi;
        
        % compute the bipartition swap probability r
        r(i,t) = exp(temps(i)*(phi_diff(i,t))); % (Eq. 41 in Kitazono et al)
        
        % accept or reject new bipartition?
        alpha = min([1, r(i,t)]);
        u = rand;
        if u < alpha
            S_t{i} = S_c{i};
            phi(i) = S_c_phi;
            full_phi(i) = S_c_phi_full;
        else
            phi(i) = S_t_phi;
            full_phi(i) = S_t_phi_full;
        end
        all_phi(i,t) = phi(i);
        
    end
    % plot current values of normalized integrated information
    plot(all_phi');
    xlabel('Iteration')
    ylabel('\Phi (normalized, non-extrapolated)');
    title({'Minimization of \Phi Across 6 Sequences, Each With Its Own Temperature',...
        '(this will pause every 5 iterations, for the first 200, to update temps)'});
    drawnow;
    
    % every 5 steps, exchange temperatures according to the Metropolis
    % criterion
    if mod(t,5)==0
        for i = 1:length(temps)-1
            r_prime = exp((temps(i+1)-temps(i))*(phi(i+1)-phi(i))); % (Eq. 43 in Kitazono et al)
            p = min([1, r_prime]);
            u = rand;
            
            % swap or not swap temperatures across sequences?
            if u<p
                temp1 = S_t{i};
                temp2 = S_t{i+1};
                S_t{i} = temp2;
                S_t{i+1} = temp1;
            end
        end
    end
    
    % update temperatures every 5 steps, until the 200th step
    if mod(t,5)==0 && t<200
        old_temps=temps;
        
        % regress mean of normalized integrated information as a function
        % of temperature
        mean_func = fit(temps', nanmean(all_phi(:,(t-4):t),2), 'pchipinterp');
        
        % regress variance of normalized integrated information as a function
        % of temperature (add a small positive constant to avoid dividing by
        % zero in case there were no partition changes at a particular
        % temperature)
        var_func = fit(temps', nanstd(all_phi(:,(t-4):t),[],2).^2+.00001, 'pchipinterp');
        
        % store changes in normalized integrated information
        changes = abs(phi_diff(phi_diff<=0));
        
        % calculate new low and high temperatures, so that the average
        % acceptance ratio for the high temperature is 0.01 and the average
        % acceptance ratio for the low temperature is 0.5
        clear rat
        temps = 1:100000:1000000;
        for i = 1:length(temps)
            rat(i) = mean(exp(-temps(i)*changes));
        end
        top_temp_f=fit(temps',rat'-.01,'pchipinterp');
        new_temps(1) = fzero(top_temp_f,[0 max(temps)]);
        
        bottom_temp_f=fit(temps',rat'-.5,'pchipinterp');
        new_temps(6) = fzero(bottom_temp_f,[0 max(temps)]);
        clear temps
        
        % set the intermediate temperatures by minimizing the following
        % cost function (Eqs. 45 & 46 in Kitazono et al):
        e_func = @(x)nansum([(.5*erfc((mean_func(x(2))-mean_func(x(1)))/...
            (sqrt(2*(var_func(x(2))+var_func(x(1))))))...
            + (1-.5*erfc((mean_func(x(2))-mean_func(x(1)))/...
            (sqrt(2*(var_func(x(2))+var_func(x(1)))))))...
            *exp((x(1)-x(2))*(mean_func(x(1))-mean_func(x(2)))))^(-4)...
            ...
            (.5*erfc((mean_func(x(3))-mean_func(x(2)))/...
            (sqrt(2*(var_func(x(3))+var_func(x(2))))))...
            + (1-.5*erfc((mean_func(x(3))-mean_func(x(2)))/...
            (sqrt(2*(var_func(x(3))+var_func(x(2)))))))...
            *exp((x(2)-x(3))*(mean_func(x(2))-mean_func(x(3)))))^(-4)...
            ...
            (.5*erfc((mean_func(x(4))-mean_func(x(3)))/...
            (sqrt(2*(var_func(x(4))+var_func(x(3))))))...
            + (1-.5*erfc((mean_func(x(4))-mean_func(x(3)))/...
            (sqrt(2*(var_func(x(4))+var_func(x(3)))))))...
            *exp((x(3)-x(4))*(mean_func(x(3))-mean_func(x(4)))))^(-4)...
            ...
            (.5*erfc((mean_func(x(5))-mean_func(x(4)))/...
            (sqrt(2*(var_func(x(5))+var_func(x(4))))))...
            + (1-.5*erfc((mean_func(x(5))-mean_func(x(4)))/...
            (sqrt(2*(var_func(x(5))+var_func(x(4)))))))...
            *exp((x(4)-x(5))*(mean_func(x(4))-mean_func(x(5)))))^(-4)...
            ...
            (.5*erfc((mean_func(x(6))-mean_func(x(5)))/...
            (sqrt(2*(var_func(x(6))+var_func(x(5))))))...
            + (1-.5*erfc((mean_func(x(6))-mean_func(x(5)))/...
            (sqrt(2*(var_func(x(6))+var_func(x(5)))))))...
            *exp((x(5)-x(6))*(mean_func(x(5))-mean_func(x(6)))))^(-4)]);
        
        % set the lower and upper bounds for the minimization. These
        % bounds are either set according to the highest and lowest
        % temperatures set above, or set so as to avoid temperatures
        % that, when fed into the variance function defined above,
        % yield negative variance)
        a = new_temps(1);
        b = new_temps(6);
        
        if var_func(b)>0
            lower_bound = repmat(b, 1,6);
        else
            first_vals = var_func(0:old_temps(1));
            if any(first_vals<0)
                first_ind = find(first_vals<0,1)-1;
                lower_bound = repmat(fzero(var_func,[first_ind old_temps(1)]), 1,6);
            else
                lower_bound = zeros(1,6);
            end
        end
        
        if var_func(a)>0
            upper_bound = repmat(a,1,6);
        else
            if any(var_func((lower_bound(1)+1):(a/100):(a*100))<=0)
                upper_bound = repmat(fzero(var_func,[lower_bound(1)+1, a*100])-1,1,6);
            else
                upper_bound = [];
            end
        end
        
        % set new temperatures by minimizing the cost function:
        evalc('temps = fmincon(e_func, old_temps,[],[],[],[], lower_bound, upper_bound);');
        
        temps(1)=new_temps(1);
        temps(6)=new_temps(6);
    end
    
    % after the 300th iteration, start calculating the convergence
    % criterion defined in Brooks & Gelman (1998). Terminate the script if
    % all sequences reach a convergence criterion Rc of 1.01
    if t>300
        new_phi=all_phi(:,201:end);
        n=size(new_phi,2);
        
        for i = 1:size(new_phi,1)
            s1 = new_phi(i,1:round(n/2));
            s2 = new_phi(i,round(n/2)+1:end);
            
            B=(((nanmean(s1)-nanmean([s1 s2]))^2)+((nanmean(s2)-nanmean([s1 s2]))^2))*n;
            W = (nansum((s1-nanmean(s1)).^2)+nansum((s2-nanmean(s2)).^2))/(2*(n-1));
            
            var = (n-1)/n*W + B/n;
            V(i) = var + B/(2*n);
            
            R(i) = V(i)/W;
        end
        
        d = 2*V/(std(V).^2);
        Rc = (d+3)/(d+1).*R;
        if all(Rc<1.01)
            converge = 1;
        end
    end
    
    
    % store output
    output.phi = full_phi;
    output.phi_norm = phi;
    output.partitions = S_t;
    % store interrupt output in case you terminate this script through
    % ctrl+c
    interrupt_output = output;
    t=t+1;
end
end

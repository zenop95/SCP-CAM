function tcaEstimate()
%  Compute the Time of Closest Approach from the orthogonality between the 
%  relative position and the relative velocity.
%     
%      INPUT:
%      - mean_state_1: primary object mean state in ECI reference frame at
%                      TCA0_et
%      - mean_state_2: secondary object mean state in ECI reference frame at
%                      TCA0_et
%       - TCA0_et     : TCA first guess (in ET)
%  
%       OUTPUT:
%       - TCA_et_out: computed TCA (in ET)
% 
%       Author:  Andrea De Vittori, Politecnico di Milano, 21 March 2022
%               e-mail: andrea.devittori@polimi.it
    % Lists definition
    err_tot = [];
    time_list = [];
    name = [];
    
    % Path for CDMs
    path = 'CDM_Python_Perturbed';
    CDMs = dir(fullfile(path, '*.txt'));

    % Load spice kernel
    cspice_furnsh('latest_leapseconds.tls');

    % Parameters definition
    t0_et = 6e8;
    mu = 398600.4415;
    t0_utc = t0_et - cspice_deltet(t0_et, 'ET');
    
    flag_quartic = 1;

    % Method for picardlindelof
    coeff_extractor = picardlindelof();
   
    % Loop over the CDMs 
    for i = 1:length(CDMs)

        % Extract array for the Primary and Secondary states
        arr = load(fullfile(path, CDMs(i).name));
        dim = size(arr, 1);

        % Find the initial delta_t
        array = strsplit(CDMs(i).name, '_');
        delta_t = str2double(strrep(array{2}, '.txt', ''));
        name = [name, delta_t];

        % Define the TCA guess
        t0_utc_guess = t0_utc - delta_t;
        t0_et_guess = t0_utc_guess + cspice_deltet(t0_utc_guess, 'UTC');
        start = tic();

        % ridefine err at each loop on i
        err = [];
        for j = 1:dim

            % extract the primary and secondary states
            primary = arr(j, 1:6);
            secondary = arr(j, 7:end);
            if flag_quartic

                % Find TCA_et_new/TCA_UTC_new and Primary/Secondary state with the Quartic formula approximation
                [a0, a1, a2, a3, a4] = coeff_extractor.extractcoefficients(primary, secondary);
                [z1, z2, z3, z4] = quartic_solver(a0/a4, a1/a4, a2/a4, a3/a4);
                solutions = [z1, z2, z3, z4];
                real_solutions = solutions(imag(solutions) == 0);
                [~, index] = min(abs(real_solutions));
                computed_tcas = real_solutions(index);
                elems_1 = cspice_oscelt(primary, t0_et_guess, 398600.4415);
                elems_2 = cspice_oscelt(secondary, t0_et_guess, 398600.4415);
                t0_utc_guess_new = t0_utc_guess + computed_tcas;
                t0_et_guess_new = t0_utc_guess_new + cspice_deltet(t0_utc_guess_new, 'UTC');
                primary = elems2state(elems_1, t0_et_guess_new, computed_tcas);
                secondary = elems2state(elems_2, t0_et_guess_new, computed_tcas);

            end

            % keplerian TCA finding
            TCA_out = t0_utc_guess_new; %find_TCA_ET(primary, secondary,  t0_et_guess_new);

            % Append the TCA error
            err = [err, abs(TCA_out - t0_utc)];
        end

        % Append the overall error and computational time
        err_tot = [err_tot; err];
        time_list = [time_list, toc(start)];
    end

    % Plot the computational time vs TCA_guess
    figure;
    scatter(name, time_list, 1000, 'black', 'filled');
    xlabel('TCA guess [s]', 'FontSize', 30);
    ylabel('Computational time for 2170 cases [s]', 'FontSize', 30);
    set(gca, 'FontSize', 30);
    grid on;
    set(gca, 'GridAlpha', 0.5);

    % Plot the mean error vs TCA_guess
    figure;
    for i = 1:length(CDMs)
        array = strsplit(CDMs(i).name, '_');
        delta_t = str2double(strrep(array{2}, '.txt', ''));
        disp(delta_t);
        scatter(delta_t, mean(err_tot(i, :)), 500, 'black', 'filled');
        hold on;
    end
    hold off;
    xlabel('TCA guess [s]', 'FontSize', 30);
    ylabel('TCA mean estimation error [s]', 'FontSize', 30);
    set(gca, 'FontSize', 30);
    grid on;
    set(gca, 'GridAlpha', 0.5);
end

% ... (Include the other functions such as picardlindelof, fsolve_fun, elems2state, etc.)


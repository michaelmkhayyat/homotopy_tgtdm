import casadi.*

% In this planning technique, we just want each point to lie outside of the
% convex region described by the 6 inequalities. No homotopies are defined.
% We are just seeking a motion planning in a nonconvex space.

% Allocating the variables
bigM = 10000;

nohomotopy_setup = [];

nohomotopy_setup.vars.sigma.sigma_A = cell(1,size(playersPairs,1));
nohomotopy_setup.vars.sigma.sigma_B = cell(1,size(playersPairs,1));
nohomotopy_setup.vars.sigma.sigma_C = cell(1,size(playersPairs,1));
nohomotopy_setup.vars.sigma.sigma_D = cell(1,size(playersPairs,1));
nohomotopy_setup.vars.sigma.sigma_E = cell(1,size(playersPairs,1));
nohomotopy_setup.vars.sigma.sigma_F = cell(1,size(playersPairs,1));
nohomotopy_setup.discreteflags.sigma.sigma_A = cell(1,size(playersPairs,1));
nohomotopy_setup.discreteflags.sigma.sigma_B = cell(1,size(playersPairs,1));
nohomotopy_setup.discreteflags.sigma.sigma_C = cell(1,size(playersPairs,1));
nohomotopy_setup.discreteflags.sigma.sigma_D = cell(1,size(playersPairs,1));
nohomotopy_setup.discreteflags.sigma.sigma_E = cell(1,size(playersPairs,1));
nohomotopy_setup.discreteflags.sigma.sigma_F = cell(1,size(playersPairs,1));

nohomotopy_setup.constraints.general.constraints =  cell(1,size(playersPairs,1));
nohomotopy_setup.constraints.general.ub =  cell(1,size(playersPairs,1));
nohomotopy_setup.constraints.general.lb =  cell(1,size(playersPairs,1));

nohomotopy_setup.constraints.vars.constraints =  cell(1,size(playersPairs,1));
nohomotopy_setup.constraints.vars.ub =  cell(1,size(playersPairs,1));
nohomotopy_setup.constraints.vars.lb =  cell(1,size(playersPairs,1));

for i = 1:1:size(playersPairs,1)

    tempconstraints = [];
    tempub = [];
    templb = [];

    varsconstraints = [];
    varsub = [];
    varslb = [];

    if intersectionBetweenEnvelopesOfPlayers(i) > 0
        nohomotopy_setup.vars.h{i} = binvar(1,1);

        nohomotopy_setup.vars.sigma.sigma_A{i} = binvar(players{playersPairs(i,2)}.opt.params.N+1,1);

        nohomotopy_setup.vars.sigma.sigma_B{i} = binvar(players{playersPairs(i,2)}.opt.params.N+1,1);
        
        nohomotopy_setup.vars.sigma.sigma_C{i} = binvar(players{playersPairs(i,2)}.opt.params.N+1,1);

        nohomotopy_setup.vars.sigma.sigma_D{i} = binvar(players{playersPairs(i,2)}.opt.params.N+1,1);

        nohomotopy_setup.vars.sigma.sigma_E{i} = binvar(players{playersPairs(i,2)}.opt.params.N+1,1);
        
        nohomotopy_setup.vars.sigma.sigma_F{i} = binvar(players{playersPairs(i,2)}.opt.params.N+1,1);

        x_1 = players{playersPairs(i,1)}.opt.vars.x;
        u_1 = players{playersPairs(i,1)}.opt.vars.u;

        x_2 = players{playersPairs(i,2)}.opt.vars.x;
        u_2 = players{playersPairs(i,2)}.opt.vars.u;

        sigma_a = nohomotopy_setup.vars.sigma.sigma_A{i};
        sigma_b = nohomotopy_setup.vars.sigma.sigma_B{i};
        sigma_c = nohomotopy_setup.vars.sigma.sigma_C{i};
        sigma_d = nohomotopy_setup.vars.sigma.sigma_D{i};
        sigma_e = nohomotopy_setup.vars.sigma.sigma_E{i};
        sigma_f = nohomotopy_setup.vars.sigma.sigma_F{i};

        for k = 1:1:players{playersPairs(i,1)}.opt.params.N+1
            
            %%%% CONSTRAINTS EQ.13 %%%%%
            % 13a
            tempconstraints = [tempconstraints; x_1(k,1) + bigM*sigma_d(k)];
            tempub = [tempub; s1_l_entry{i} + bigM];
            templb = [templb; -inf];

            % 13b
            tempconstraints = [tempconstraints; x_1(k,1) + bigM*sigma_d(k)];
            tempub = [tempub; inf];
            templb = [templb; s1_l_entry{i}];

            % 13c
            tempconstraints = [tempconstraints; x_2(k,1) + bigM*sigma_a(k)];
            tempub = [tempub; s2_l_entry{i} + bigM];
            templb = [templb; -inf];

            % 13d
            tempconstraints = [tempconstraints; x_2(k,1) + bigM*sigma_a(k)];
            tempub = [tempub; inf];
            templb = [templb; s2_l_entry{i}];


            % 13e
            tempconstraints = [tempconstraints; x_1(k,1) - bigM*sigma_c(k)];
            tempub = [tempub; inf];
            templb = [templb; s1_h_exit{i} - bigM];

            % 13f
            tempconstraints = [tempconstraints; x_1(k,1) - bigM*sigma_c(k)];
            tempub = [tempub; s1_h_exit{i}];
            templb = [templb; -inf];

            % 13g
            tempconstraints = [tempconstraints; x_2(k,1) - bigM*sigma_f(k)];
            tempub = [tempub; inf];
            templb = [templb; s2_h_exit{i} - bigM];

            % 13h
            tempconstraints = [tempconstraints; x_2(k,1) - bigM*sigma_f(k)];
            tempub = [tempub; s2_h_exit{i}];
            templb = [templb; -inf];

            % 13i
            tempconstraints = [tempconstraints; x_1(k,1) - x_2(k,1) + bigM*sigma_e(k)];
            tempub = [tempub; bigM - (s2_h_entry{i} - s1_l_entry{i})];
            templb = [templb; -inf];

            % 13j
            tempconstraints = [tempconstraints; x_1(k,1) - x_2(k,1) + bigM*sigma_e(k)];
            tempub = [tempub; inf];
            templb = [templb; - (s2_h_entry{i} - s1_l_entry{i})];
            
            % 13k
            tempconstraints = [tempconstraints; x_2(k,1) - x_1(k,1) + bigM*sigma_b(k)];
            tempub = [tempub; bigM + (s2_l_entry{i} - s1_h_entry{i})];
            templb = [templb; -inf];

            % 13l
            tempconstraints = [tempconstraints; x_2(k,1) - x_1(k,1) + bigM*sigma_b(k)];
            tempub = [tempub; inf];
            templb = [templb; (s2_l_entry{i} - s1_h_entry{i})];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%

            tempconstraints = [tempconstraints; sigma_a(k)  + sigma_b(k) + sigma_c(k) + sigma_d(k)  + sigma_e(k) +  sigma_f(k)];
            tempub = [tempub; inf];
            templb = [templb; 1];

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if k < players{playersPairs(i,1)}.opt.params.N+1
                tempconstraints = [tempconstraints; sigma_a(k+1) - sigma_a(k)];
                tempub = [tempub; 0];
                templb = [templb; -1];
    
                tempconstraints = [tempconstraints; sigma_c(k+1) - sigma_c(k)];
                tempub = [tempub; 1];
                templb = [templb; 0];
                
                tempconstraints = [tempconstraints; sigma_f(k+1) - sigma_f(k)];
                tempub = [tempub; 1];
                templb = [templb; 0];
                
                tempconstraints = [tempconstraints; sigma_d(k+1) - sigma_d(k)];
                tempub = [tempub; 0];
                templb = [templb; -1];
            end

        end

        % Constraints on sigmas
        varsconstraints = [varsconstraints; sigma_a(:)];
        varsub = [varsub; ones(size(sigma_a(:)))];
        varslb = [varslb; zeros(size(sigma_a(:)))];

        varsconstraints = [varsconstraints; sigma_b(:)];
        varsub = [varsub; ones(size(sigma_b(:)))];
        varslb = [varslb; zeros(size(sigma_b(:)))];

        varsconstraints = [varsconstraints; sigma_c(:)];
        varsub = [varsub; ones(size(sigma_c(:)))];
        varslb = [varslb; zeros(size(sigma_c(:)))];

        varsconstraints = [varsconstraints; sigma_d(:)];
        varsub = [varsub; ones(size(sigma_d(:)))];
        varslb = [varslb; zeros(size(sigma_d(:)))];

        varsconstraints = [varsconstraints; sigma_e(:)];
        varsub = [varsub; ones(size(sigma_e(:)))];
        varslb = [varslb; zeros(size(sigma_e(:)))];

        varsconstraints = [varsconstraints; sigma_f(:)];
        varsub = [varsub; ones(size(sigma_f(:)))];
        varslb = [varslb; zeros(size(sigma_f(:)))];

        nohomotopy_setup.constraints.general.constraints{i} = tempconstraints;
        nohomotopy_setup.constraints.general.ub{i} = tempub;
        nohomotopy_setup.constraints.general.lb{i} = templb;

        nohomotopy_setup.constraints.vars.constraints{i} = varsconstraints;
        nohomotopy_setup.constraints.vars.ub{i} = varsub;
        nohomotopy_setup.constraints.vars.lb{i} = varslb;
    end
    

end
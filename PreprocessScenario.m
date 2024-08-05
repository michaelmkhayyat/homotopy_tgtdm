playersPairs = nchoosek(1:1:n,2);
intersectionBetweenEnvelopesOfPlayers = zeros(size(playersPairs,1),1);
pointsGlobalCoordinates = cell(1,size(playersPairs,1));

p_star1 = cell(1,size(playersPairs,1));
p_star2 = cell(1,size(playersPairs,1));
s1_l_entry = cell(1,size(playersPairs,1));
s1_h_entry = cell(1,size(playersPairs,1));
s2_l_entry = cell(1,size(playersPairs,1));
s2_h_entry = cell(1,size(playersPairs,1));
s1_l_exit = cell(1,size(playersPairs,1));
s1_h_exit = cell(1,size(playersPairs,1));
s2_l_exit = cell(1,size(playersPairs,1));
s2_h_exit = cell(1,size(playersPairs,1));
sgn = cell(1,size(playersPairs,1));

for i = 1:1:size(playersPairs,1)

    pointsToCheckForIntersections = [];


    pointsUpperBoundaryPlayer1 = players{playersPairs(i,1)}.pathInfo.upperBound;
    pointsLowerBoundaryPlayer1 = players{playersPairs(i,1)}.pathInfo.lowerBound;

    pointsUpperBoundaryPlayer2 = players{playersPairs(i,2)}.pathInfo.upperBound;
    pointsLowerBoundaryPlayer2 = players{playersPairs(i,2)}.pathInfo.lowerBound;

    % Upper and Upper
    [pointsToCheckForIntersections{1}(:,1),pointsToCheckForIntersections{1}(:,2)] = intersections(pointsUpperBoundaryPlayer1(:,1), pointsUpperBoundaryPlayer1(:,2), pointsUpperBoundaryPlayer2(:,1), pointsUpperBoundaryPlayer2(:,2));
    % Upper and Lower
    [pointsToCheckForIntersections{2}(:,1),pointsToCheckForIntersections{2}(:,2)] = intersections(pointsUpperBoundaryPlayer1(:,1), pointsUpperBoundaryPlayer1(:,2), pointsLowerBoundaryPlayer2(:,1), pointsLowerBoundaryPlayer2(:,2));
    % Lower and Upper
    [pointsToCheckForIntersections{3}(:,1),pointsToCheckForIntersections{3}(:,2)] = intersections(pointsLowerBoundaryPlayer1(:,1), pointsLowerBoundaryPlayer1(:,2), pointsUpperBoundaryPlayer2(:,1), pointsUpperBoundaryPlayer2(:,2));
    % Lower and Lower
    [pointsToCheckForIntersections{4}(:,1),pointsToCheckForIntersections{4}(:,2)] = intersections(pointsLowerBoundaryPlayer1(:,1), pointsLowerBoundaryPlayer1(:,2), pointsLowerBoundaryPlayer2(:,1), pointsLowerBoundaryPlayer2(:,2));

    intersectionBetweenEnvelopesOfPlayers(i) = any([~isempty(pointsToCheckForIntersections{1}),~isempty(pointsToCheckForIntersections{2}),~isempty(pointsToCheckForIntersections{3}),~isempty(pointsToCheckForIntersections{4})]);

    if intersectionBetweenEnvelopesOfPlayers(i) > 0

        pointsToCheckForIntersections = pointsToCheckForIntersections(~cellfun('isempty',pointsToCheckForIntersections));
        
        tempPointsCoordinates = [];

        for j = 1:1:length(pointsToCheckForIntersections)
           
            if size(pointsToCheckForIntersections{j},1) == 1
                tempPointsCoordinates = [tempPointsCoordinates; pointsToCheckForIntersections{j}(1,:)];
            else
                tempPointsCoordinates = [tempPointsCoordinates;pointsToCheckForIntersections{j}(1,:); pointsToCheckForIntersections{j}(end,:)];
            end

        end

        pointsGlobalCoordinates{i} = tempPointsCoordinates;

        temp_1 = [];
        temp_2 = [];
        
        for j = 1:1:size(pointsGlobalCoordinates{i},1)
                t = global2frenet(players{playersPairs(i,1)}.referencePath, [pointsGlobalCoordinates{i}(j,1) pointsGlobalCoordinates{i}(j,2) 0 0 0 0]);
                temp_1(j) = t(1);
                t = global2frenet(players{playersPairs(i,2)}.referencePath, [pointsGlobalCoordinates{i}(j,1) pointsGlobalCoordinates{i}(j,2) 0 0 0 0]);
                temp_2(j) = t(1);                
        end

        [~,idxMin_1] = min(temp_1);
        [~,idxMax_1] = max(temp_1);

        [~,idxMin_2] = min(temp_2);
        [~,idxMax_2] = max(temp_2);

        if mod(size(pointsGlobalCoordinates{i},1),2) == 0
            p_star1{i} = pointsGlobalCoordinates{i}(idxMin_1,:);
            p_star2{i} = pointsGlobalCoordinates{i}(idxMax_1,:);

            if abs(temp_2(idxMin_1) - temp_2(idxMax_1)) < 0.1
                p_star2{i} = pointsGlobalCoordinates{i}(idxMax_2,:);
            end

        else
            p_star1{i} = pointsGlobalCoordinates{i}(idxMin_1,:);
            p_star2{i} = [];
        end
        playersPairs(i,1)
        playersPairs(i,2)

        if isempty(p_star2{i})

            temp1_pstar1 = global2frenet(players{playersPairs(i,1)}.referencePath, [p_star1{i}(1) p_star1{i}(2) 0 0 0 0]);
            temp2_pstar1 = global2frenet(players{playersPairs(i,2)}.referencePath, [p_star1{i}(1) p_star1{i}(2) 0 0 0 0]);

            s1_l_entry{i} = temp1_pstar1(1)-0.5*players{playersPairs(i,1)}.params.safetyLength;
            s1_h_entry{i} = temp1_pstar1(1)+0.5*players{playersPairs(i,1)}.params.safetyLength;
            s2_l_entry{i} = temp2_pstar1(1)-0.5*players{playersPairs(i,1)}.params.safetyLength;
            s2_h_entry{i} = temp2_pstar1(1)+0.5*players{playersPairs(i,2)}.params.safetyLength;

            s1_l_exit{i} = 1000;%players{playersPairs(i,1)}.referencePath.PathLength; %- 0.5*players{playersPairs(i,1)}.params.safetyLength;
            s1_h_exit{i} = 1000;%players{playersPairs(i,1)}.referencePath.PathLength; %+ 0.5*players{playersPairs(i,1)}.params.safetyLength;
            s2_l_exit{i} = 1000;%players{playersPairs(i,2)}.referencePath.PathLength; %- 0.5*players{playersPairs(i,2)}.params.safetyLength;
            s2_h_exit{i} = 1000;%players{playersPairs(i,2)}.referencePath.PathLength; %+ 0.5*players{playersPairs(i,2)}.params.safetyLength;

            sgn{i} = 1;

        else

            temp1_pstar1 = global2frenet(players{playersPairs(i,1)}.referencePath, [p_star1{i}(1) p_star1{i}(2) 0 0 0 0]);
            temp2_pstar1 = global2frenet(players{playersPairs(i,2)}.referencePath, [p_star1{i}(1) p_star1{i}(2) 0 0 0 0]);

            temp1_pstar2 = global2frenet(players{playersPairs(i,1)}.referencePath, [p_star2{i}(1) p_star2{i}(2) 0 0 0 0]);
            temp2_pstar2 = global2frenet(players{playersPairs(i,2)}.referencePath, [p_star2{i}(1) p_star2{i}(2) 0 0 0 0]);

            if sign(temp1_pstar1(1) - temp1_pstar2(1)) == sign(temp2_pstar1(1) - temp2_pstar2(1))
                
                s1_l_entry{i} = temp1_pstar1(1)-0.5*players{playersPairs(i,1)}.params.safetyLength;
                s1_h_entry{i} = temp1_pstar1(1)+0.5*players{playersPairs(i,1)}.params.safetyLength;
                s2_l_entry{i} = temp2_pstar1(1)-0.5*players{playersPairs(i,2)}.params.safetyLength;
                s2_h_entry{i} = temp2_pstar1(1)+0.5*players{playersPairs(i,2)}.params.safetyLength;

                s1_l_exit{i} = temp1_pstar2(1)-0.5*players{playersPairs(i,1)}.params.safetyLength;
                s1_h_exit{i} = temp1_pstar2(1)+0.5*players{playersPairs(i,1)}.params.safetyLength;
                s2_l_exit{i} = temp2_pstar2(1)-0.5*players{playersPairs(i,2)}.params.safetyLength;
                s2_h_exit{i} = temp2_pstar2(1)+0.5*players{playersPairs(i,2)}.params.safetyLength;

                if abs(temp1_pstar1(1) - temp1_pstar2(1)) < players{playersPairs(i,1)}.params.safetyLength
                    s1_l_exit{i} = s1_l_entry{i};                    
                    s1_h_entry{i} = s1_h_exit{i};
                end

                if abs(temp2_pstar1(1) - temp2_pstar2(1)) < players{playersPairs(i,2)}.params.safetyLength
                    s2_h_entry{i} = s2_h_exit{i};
                    s2_l_exit{i} = s2_l_entry{i};
                end

                sgn{i} = 1;
            else

                sgn{i} = -1;

                sorted_p1 = sort([temp1_pstar1; temp1_pstar2]);
                sorted_p2 = sort([temp2_pstar1; temp2_pstar2]);

                s1_l_entry{i} = sorted_p1(1)-0.5*players{playersPairs(i,1)}.params.safetyLength;
                s1_h_entry{i} = sorted_p1(1)+0.5*players{playersPairs(i,1)}.params.safetyLength;
                s2_l_entry{i} = sorted_p2(1)-0.5*players{playersPairs(i,2)}.params.safetyLength;
                s2_h_entry{i} = sorted_p2(1)+0.5*players{playersPairs(i,2)}.params.safetyLength;

                s1_l_exit{i} = sorted_p1(2)-0.5*players{playersPairs(i,1)}.params.safetyLength;
                s1_h_exit{i} = sorted_p1(2)+0.5*players{playersPairs(i,1)}.params.safetyLength;
                s2_l_exit{i} = sorted_p2(2)-0.5*players{playersPairs(i,2)}.params.safetyLength;
                s2_h_exit{i} = sorted_p2(2)+0.5*players{playersPairs(i,2)}.params.safetyLength;

                s1_l_exit{i} = s1_l_entry{i};
                s1_h_entry{i} = s1_h_exit{i};

                s2_l_exit{i} = s2_l_entry{i};
                s2_h_entry{i} = s2_h_exit{i};
                
            end


        end

        t = global2frenet(players{playersPairs(i,2)}.referencePath, [pointsGlobalCoordinates{i}(j,1) pointsGlobalCoordinates{i}(j,2) 0 0 0 0]);

    end

                temp1_pstar1
            temp2_pstar1

            temp1_pstar2 
            temp2_pstar2 

end
%%

%/* Now, we initialize the homotopy classes variables and the corresponding
%   inequalities /*

% Constraints Method 1 setups up the constraints as described in Eq.13
constraints_method1()

% Constraints Method 2 setups up the constraints as described in Eq.17
%constraints_method2()

import casadi.*

J = 0;

state_vars = [];
state_vars_ub = [];
state_vars_lb = [];
state_vars_discreteflags = [];

control_vars = [];
control_vars_ub = [];
control_vars_lb = [];
control_vars_discreteflags = [];

parameter_vars =[];

dynamics = [];
dynamics_ub = [];
dynamics_lb = [];

homotopy_vars = [];
homotopy_vars_ub = [];
homotopy_vars_lb = [];
homotopy_vars_discreteflags = [];

homotopy_constraints = [];
homotopy_constraints_ub = [];
homotopy_constraints_lb = [];

x0_init = [];
x0_states = [];
x0_controls = [];

rf_vars = [];
rf_lb = [];
rf_ub = [];

for i = 1:1:n

    % Summing all personal cost functions
    J = J + players{i}.opt.costFunction;

    % Concatenating the state variables and their bounds for all players 

    state_vars = [state_vars; players{i}.opt.constraints.state.vars(:)];
    state_vars_ub = [state_vars_ub; players{i}.opt.constraints.state.ub(:)];
    state_vars_lb = [state_vars_lb; players{i}.opt.constraints.state.lb(:)];
    state_vars_discreteflags = [state_vars_discreteflags; players{i}.opt.discreteflag.x(:)];

    % Concatenating the control variables and their bounds  for all players 
    control_vars = [control_vars; players{i}.opt.constraints.control.vars(:)];
    control_vars_ub = [control_vars_ub; players{i}.opt.constraints.control.ub(:)];
    control_vars_lb = [control_vars_lb; players{i}.opt.constraints.control.lb(:)];
    control_vars_discreteflags = [control_vars_discreteflags; players{i}.opt.discreteflag.u(:)];

    % Concatenating the parameter variables
    parameter_vars = [parameter_vars; players{i}.opt.vars.p(:)];
    
    % Concatenating the dynamics constraints and their bounds for all players.
    dynamics = [dynamics; players{i}.opt.constraints.dynamics(:)];
    dynamics_ub = [dynamics_ub; zeros(length(players{i}.opt.constraints.dynamics(:)),1)];
    dynamics_lb = [dynamics_lb; zeros(length(players{i}.opt.constraints.dynamics(:)),1)];

    x0_init = [x0_init;  players{i}.initialStateFrenet(:)];

    tempGuess_states =  repmat(players{i}.initialStateFrenet,players{i}.opt.params.N+1,1);
    tempGuess_controls = zeros(players{i}.opt.params.N,1);

    x0_states = [x0_states; tempGuess_states(:)];
    x0_controls = [x0_controls; tempGuess_controls(:)];

  %  rf_vars = [rf_vars;players{i}.opt.recursiveFeasabilityConstraints.vars];
 %   rf_lb = [rf_lb;players{i}.opt.recursiveFeasabilityConstraints.lb];
 %   rf_ub = [rf_ub;players{i}.opt.recursiveFeasabilityConstraints.ub];
end

x0_homotopy = [];

h_variables =  [];
sigma_variables = [];
slack_variables = [];


varsFieldNames = fieldnames(homotopy_setup.vars);
sigmaFieldNames = fieldnames(homotopy_setup.vars.sigma);
slacksExist = any(strcmp(varsFieldNames,'slack'));
if slacksExist
    slackFieldNames = fieldnames(homotopy_setup.vars.slack);
end
for i = 1:1:size(playersPairs,1)

    if intersectionBetweenEnvelopesOfPlayers(i) > 0

        h_variables = [h_variables; homotopy_setup.vars.h{i}];

        for k = 1:1:length(sigmaFieldNames)

            tmp = homotopy_setup.vars.sigma.(sigmaFieldNames{k})(i);
            
            sigma_variables = [sigma_variables; tmp{:}];

        end
        
        if slacksExist
            for k = 1:1:length(slackFieldNames)
                 tmp = homotopy_setup.vars.slack.(slackFieldNames{k})(i);
                slack_variables = [slack_variables; tmp{:}];
            end
        end
        
    
        homotopy_vars = [homotopy_vars; homotopy_setup.constraints.vars.constraints{i}(:)];
        homotopy_vars_ub = [homotopy_vars_ub; homotopy_setup.constraints.vars.ub{i}(:) ];
        homotopy_vars_lb = [homotopy_vars_lb; homotopy_setup.constraints.vars.lb{i}(:)];

        homotopy_constraints = [homotopy_constraints; homotopy_setup.constraints.general.constraints{i}(:) ];
        homotopy_constraints_ub = [homotopy_constraints_ub; homotopy_setup.constraints.general.ub{i}(:)];
        homotopy_constraints_lb = [homotopy_constraints_lb; homotopy_setup.constraints.general.lb{i}(:)];

        tempGuess_homotopy = zeros(length(homotopy_setup.constraints.vars.constraints{i}(:)),1);

        x0_homotopy = [x0_homotopy; tempGuess_homotopy];

    end
end

allVariables = [state_vars; control_vars; h_variables; sigma_variables];
if slacksExist
    allVariables = [allVariables; slack_variables];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Tactical Game-theoretic Decision-making with Homotopy Class Constraints %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc
%% Setting optimization parameters
N = 35;
dT = 0.1;
totalSimulationTime = 40;

% plot(scenario)

%% Initializing the players
% Vector of nice colors (red, blue, green, orange, yellow, purple, cyan)
col = {[0.6350 0.0780 0.1840], [0 0.4470 0.7410], [0.4660 0.6740 0.1880], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560], [0.3010 0.7450 0.9330]};

load("exportedTracks_00.mat")

% Note these must be modified according to the scenario. These specify the
% entry points to the roundabout.

%for Map 00
entryPointsForPlayers = [56,20,54,27];
exitPointsForPlayers = [86, 69, 102, 42];

players{1} = Player_rounD(1, [exportedtracks{1}.xC(1:2:end), exportedtracks{1}.yC(1:2:end)], [40 2.5], [40 4.5], 1.5, 3.6, 1.6, col{1}, N, dT, meterPerPixel);
players{2} = Player_rounD(2, [exportedtracks{2}.xC(1:2:end), exportedtracks{2}.yC(1:2:end)], [8 3], [40 3.5], 1.5, 3.6, 1.6, col{2}, N, dT, meterPerPixel);
players{3} = Player_rounD(3, [exportedtracks{3}.xC(1:2:end), exportedtracks{3}.yC(1:2:end)], [45 1], [40 1], 1.5, 3.6, 1.6,col{3}, N, dT, meterPerPixel);
players{4} = Player_rounD(4, [exportedtracks{4}.xC(1:2:end), exportedtracks{4}.yC(1:2:end)], [15 3], [40 3], 1.5, 3.6, 1.6, col{4}, N, dT, meterPerPixel);

n = length(players);

% Preprocessing the intersection of reference path evelopes of the players

PreprocessScenario()

%% Displaying the scenario and the players

% Plotting the scenario and the important points for the intersection of 
% the envelopes of the reference paths

fullfig
axis tight
set(gca,'ydir','reverse')
backgroundImage = imread(backgroundImagePath);
backgroundImage = backgroundImage(:, :, :);
h = image([.0 .0], [.0 .0], backgroundImage);
uistack(h,'bottom')
xlim([400 1200])
ylim([50 850])
axis square
hold on
for i = 1:1:n
    plot(players{i}.pathInfo.centerPoints(:,1)/meterPerPixel, players{i}.pathInfo.centerPoints(:,2)/meterPerPixel, 'Color',players{i}.params.col, 'LineWidth',2)
    plot(players{i}.pathInfo.lowerBound(:,1)/meterPerPixel, players{i}.pathInfo.lowerBound(:,2)/meterPerPixel, 'Color',players{i}.params.col)
    plot(players{i}.pathInfo.upperBound(:,1)/meterPerPixel, players{i}.pathInfo.upperBound(:,2)/meterPerPixel, 'Color',players{i}.params.col)
    players{i}.drawPlayer([players{i}.currentStateGlobal(1)/meterPerPixel,players{i}.currentStateGlobal(2)/meterPerPixel, players{i}.currentStateGlobal(3)])

end
%saveas(gcf,'scenario.pdf')
hold on
for i = 1:1:length(p_star1)
    if ~isempty(p_star1{i})
        scatter(p_star1{i}(:,1)/meterPerPixel, p_star1{i}(:,2)/meterPerPixel,'k')
    end
    if ~isempty(p_star2{i})
        scatter(p_star2{i}(:,1)/meterPerPixel, p_star2{i}(:,2)/meterPerPixel,'k')
    end
end

% Predefine Homtopy Class (h_flag = 1) or Solve for Homotopy Class (h_flag = 0)
h_flag = 1;

%%
if h_flag == 1

    % We iterate the possible homotopy classes
    % We define vector h_possible as the possible homotopy classes between
    % players. Note that -1 indicates no interaction between players is present
    % (this is distinguish homotopy class 0)

    nonZeroElements = nnz(intersectionBetweenEnvelopesOfPlayers);

    possibleCombinations = ff2n(nonZeroElements);
    
    totalControlEffort = zeros(1,size(possibleCombinations,1));
    totalProgress = zeros(1,size(possibleCombinations,1)); 
    totalTime = cell(1,size(possibleCombinations,1));
    
    x = [state_vars; control_vars; homotopy_vars];
    lbx = [state_vars_lb; control_vars_lb; homotopy_vars_lb];
    ubx = [state_vars_ub; control_vars_ub; homotopy_vars_ub];
    p = parameter_vars;

    g = [dynamics; homotopy_constraints;scenario_specific_vars];
    lbg = [dynamics_lb; homotopy_constraints_lb;scenario_specific_lb];
    ubg = [dynamics_ub; homotopy_constraints_ub;scenario_specific_ub];

    constraints = [lbg <= g <= ubg];
    constraints = [constraints; lbx <= x <= ubx];
    initialGuessAllVariables = [x0_states; x0_controls];
    tmp = zeros(1,length(allVariables)-length(initialGuessAllVariables));
    initialGuessAllVariables = [initialGuessAllVariables(:); tmp(:)];

    outputVariables = {state_vars, control_vars,h_variables,sigma_variables,parameter_vars};
    if slacksExist
        outputVariables = {state_vars, control_vars,h_variables,sigma_variables,slack_variables};
    end

    inputVariables = {p,h_variables};
    options = sdpsettings('solver','gurobi');
    options = sdpsettings(options,'warmstart',1);
    options = sdpsettings(options,'usex0',1);
    controller = optimizer(constraints,J,options,inputVariables,outputVariables,[N n]);
%%
    for i = 1:1:size(possibleCombinations,1)
        
        finishFlag = zeros(size(players));

        homotopy_vars_values = possibleCombinations(i,:);

        fprintf(['Solving For homotopy vector ' repmat(' %1.0f ',1,numel(homotopy_vars_values)) ' in MPC Fashion \n'],homotopy_vars_values);
        
        x0_init = [];
        for j = 1:1:n
            players{j}.stateHistory{i}(1,:) = players{j}.initialStateFrenet;
            x0_init = [x0_init;  players{j}.initialStateFrenet(:)];
        end

        x0hvalues_init = {x0_init, homotopy_vars_values'};
        
        for k = 1:1:totalSimulationTime/dT

            %x0_init
            initialTime = tic;
            [output,mess,~,~,controller,c] = controller(x0hvalues_init);
            finalTime = toc(initialTime);
            totalTime{i}(k) =finalTime;

            states_optimal = output{1};
            controls_optimal = output{2};
            homotopy_optimal = output{3};
            sigma_optimal = output{4};

            if slacksExist
                slacks_optimal = output{5};
            end    

            players_states = [];
            players_controls = [];
            x0_init = [];
            x0_states = [];
            x0_controls = [];
            x0_homotopy = [];

            for j = 1:1:n
        
                tempStates = states_optimal(1:2*(N+1));
                states_optimal(1:2*(N+1)) = [];
            
                tempControls = controls_optimal(1:N);
                controls_optimal(1:N) = [];
            
                tempStates = reshape(tempStates,[N+1,2]);
            
                players_states{j} = tempStates;
                players_controls{j} = tempControls;

                players{j}.currentStateFrenet = tempStates(2,:);
                x0_init = [x0_init; players{j}.currentStateFrenet(:)];

                players{j}.controlHistory{i}(k,:) = tempControls(1);
                players{j}.stateHistory{i}(k+1,:) = tempStates(2,:);

                if players{j}.currentStateFrenet(1) > exitPointsForPlayers(j) + 2 % A margin of 2 meters is added
                    finishFlag(j) = 1;
                end

            end
            
            x0hvalues_init = {x0_init, homotopy_vars_values'};

           if all(finishFlag)
               break;
          end

        end

    end

else

    totalProgress = 0;
    totalControlEffort = 0;
    totalTime = [];
    
    x = [state_vars; control_vars; homotopy_vars];
    p = parameter_vars;

    lbx = [state_vars_lb; control_vars_lb; homotopy_vars_lb];
    ubx = [state_vars_ub; control_vars_ub; homotopy_vars_ub];
 
    g = [dynamics; homotopy_constraints;scenario_specific_vars];
    lbg = [dynamics_lb; homotopy_constraints_lb;scenario_specific_lb];
    ubg = [dynamics_ub; homotopy_constraints_ub;scenario_specific_ub];

    constraints = [lbg <= g <= ubg];
    constraints = [constraints; lbx <= x <= ubx];
    initialGuessAllVariables = [x0_states; x0_controls];

    outputVariables = {state_vars, control_vars, h_variables,sigma_variables,parameter_vars};
    if slacksExist
        outputVariables = {state_vars, control_vars, h_variables,sigma_variables,slack_variables};
    end
    inputVariables = {p};

    options = sdpsettings('solver','gurobi');
    if slacksExist
        options.gurobi.MIQCPMethod = 1;
    end
    options = sdpsettings(options,'verbose',0);
    options = sdpsettings(options,'warning',1);
    options = sdpsettings(options,'usex0',1);

    controller = optimizer(constraints,J,options,inputVariables,outputVariables,[N n]);
    
    x0_init = [];

    for j = 1:1:n
        players{j}.stateHistory{1}(1,:) = players{j}.initialStateFrenet;
        x0_init = [x0_init;  players{j}.initialStateFrenet(:)];
    end

    for k = 1:1:totalSimulationTime/dT
        
        finishFlag = zeros(size(players));
        
        initialTime = tic;
        [output,mess,~,~,controller,c] = controller(x0_init);
        finalTime = toc(initialTime);
        totalTime{1}(k) = c.solvertime;
        
        states_optimal = output{1};
        controls_optimal = output{2};
        homotopy_optimal = output{3};
        sigma_optimal = output{4};

        if slacksExist
            slacks_optimal = output{5};
        end

       fprintf(['Optimal homotopy ' repmat(' %1.0f ',1,numel(homotopy_optimal)) '\n'],homotopy_optimal);
        players_states = [];
        players_controls = [];
        x0_init = [];

        x0_states = [];
        x0_controls = [];
        x0_homotopy = [];

        for j = 1:1:n
    
            tempStates = states_optimal(1:2*(N+1));
            states_optimal(1:2*(N+1)) = [];
        
            tempControls = controls_optimal(1:N);
            controls_optimal(1:N) = [];
        
            tempStates = reshape(tempStates,[N+1,2]);
        
            players_states{j} = tempStates;
            players_controls{j} = tempControls;

            players{j}.currentStateFrenet = tempStates(2,:);
            x0_init = [x0_init; players{j}.currentStateFrenet(:)];

            players{j}.controlHistory{1}(k,:) = tempControls(1);
            players{j}.stateHistory{1}(k+1,:) = tempStates(2,:);

            tmp1 = [tempStates(2:N+1,:); [0 0]];
            tmp2 = [tempControls(2:N);0];
            x0_states = [x0_states; tmp1(:)];
            x0_controls = [x0_controls; tmp2(:)];
            
                            if players{j}.currentStateFrenet(1) > exitPointsForPlayers(j) + 2 % A margin of 2 meters is added
                    finishFlag(j) = 1;
                end
        end

        x0_homotopy = [homotopy_optimal];

        if all(finishFlag)
            break;
        end
    end
    
    for kk = 1:1:n
        totalProgress = totalProgress + (players{kk}.currentStateFrenet(1) - players{kk}.initialStateFrenet(1));
        totalControlEffort = totalControlEffort + norm(players{kk}.controlHistory{1},2);
    end

end
%%
% Visualize solution animation
 
homotopyClassToVisualize =1;
lengthOfSimulation = length(players{1}.stateHistory{homotopyClassToVisualize});
wayPoints = [];
speed = [];
for i = 1:1:n

    statesHistory  = players{i}.stateHistory{homotopyClassToVisualize};

    for k = 1:1:size(statesHistory,1)

        s_temp = frenet2global(players{i}.referencePath, [statesHistory(k,1) statesHistory(k,2) 0 0 0 0]);

        wayPoints{i}(k,:) = [s_temp(1) s_temp(2) s_temp(3)];
        speed{i}(k,:) = [s_temp(5)];

    end

end 

figure
set(gca,'ydir','reverse')
axis tight
for j = 1:1:size(statesHistory,1)
        pause(0.01)
    clf
    backgroundImage = imread(backgroundImagePath);
    backgroundImage = backgroundImage(:, :, :);
    h = image([.0 .0], [.0 .0], backgroundImage);
    uistack(h,'bottom')
    hold on
    for i = 1:1:n
        players{i}.drawPlayer([wayPoints{i}(j,1)/meterPerPixel wayPoints{i}(j,2)/meterPerPixel wayPoints{i}(j,3)])
    end
end

%% Visualize Output
time = [0:dT:(lengthOfSimulation-1)*dT];
figure
title('Acceleration Command of Players')
hold on
for i = 1:1:n
    plot(time(1:end-1), players{i}.controlHistory{homotopyClassToVisualize}(:)', 'Color',players{i}.params.col)
end
xlabel('Time [s]');
ylabel('Acceleration m/s^2');

figure
title('Velocity of Players')
hold on
for i = 1:1:n
    plot(time', players{i}.stateHistory{homotopyClassToVisualize}(:,2), 'Color',players{i}.params.col)
end
xlabel('Time [s]');
ylabel('Velocity m/s');
%%
% Plotting the collisions area between the players
for i = 1:1:size(playersPairs,1)
    if intersectionBetweenEnvelopesOfPlayers(i) ~= 0
        s1 = players{playersPairs(i,1)}.stateHistory{homotopyClassToVisualize}(:,1);
        s2 = players{playersPairs(i,2)}.stateHistory{homotopyClassToVisualize}(:,1);
        figure
        axis equal        
        axis tight
        xlim([0 100])
        ylim([0 100])
        xpoints = [s1_l_entry{i}, s1_h_entry{i}, s1_h_exit{i}, s1_h_exit{i}, s1_l_exit{i}, s1_l_entry{i}, s1_l_entry{i}];
        ypoints = [s2_l_entry{i}, s2_l_entry{i}, s2_l_exit{i}, s2_h_exit{i}, s2_h_exit{i}, s2_h_entry{i}, s2_l_entry{i}];
        plot(xpoints,ypoints)
        hold on
        scatter(s1,s2)
        xlim([0 100])
        ylim([0 100])
        xlabel(strcat('s_',num2str(players{playersPairs(i,1)}.ID,'%d')))
        ylabel(strcat('s_',num2str(players{playersPairs(i,2)}.ID,'%d')))
        title(strcat('Collisions Area Between Player ',num2str(players{playersPairs(i,1)}.ID,'%d'),' and Player ',num2str(players{playersPairs(i,2)}.ID,'%d')))
    end
end
%% Results Analysis
% The heuristic we use is the time spent in the roundabout as well as the
% total task completion time.

tableCompletionTime = [];
tableSolvingTime = [];
tableNetProgress = [];

if h_flag == 1
    totalHomotopy = size(possibleCombinations,1);
else
    totalHomotopy = 1;
end

for j = 1:1:totalHomotopy
    
    homotopyClassToAnalyze =j;
    computationalTimesOfInterest = totalTime{homotopyClassToAnalyze};
    TaskCompletionTime = 0;
    for i = 1:1:n
        sValuesOfInterest = players{i}.stateHistory{homotopyClassToAnalyze}(2:end,1);
        controlValuesOfInterest =  players{i}.controlHistory{homotopyClassToAnalyze};
    
        [ind1, ~] = find(sValuesOfInterest>=entryPointsForPlayers(i));
        indx = ind1(sValuesOfInterest(ind1,:)<=exitPointsForPlayers(i));
        ind{i} = [min(indx); max(indx)];
        range{i} = max(indx) - min(indx);
        controlEffortDuringRoundAboutTraversal{i} = norm(controlValuesOfInterest(min(indx):max(indx)),2);
        TaskCompletionTime = max(TaskCompletionTime, max(indx));
    end  
    
    SolvingTime = sum(computationalTimesOfInterest(1:TaskCompletionTime));
    controlEffortUntilTaskCompletionTime = [];
    netProgress = 0;
    for i = 1:1:n
        controlValuesOfInterest =  players{i}.controlHistory{homotopyClassToAnalyze};
        controlEffortUntilTaskCompletionTime(i) = norm(controlValuesOfInterest(1:TaskCompletionTime),2);
        netProgress = netProgress +  players{i}.stateHistory{homotopyClassToAnalyze}(TaskCompletionTime,1) - players{i}.stateHistory{homotopyClassToAnalyze}(1,1);
    end
    
    netControlEffort(j) = sum(controlEffortUntilTaskCompletionTime);

    tableCompletionTime(j) = dT*TaskCompletionTime;
    tableSolvingTime(j) = SolvingTime;
    tableNetProgress(j) = netProgress;
    clc
end

if h_flag == 1
    resultsAnalysisTable = table(possibleCombinations,tableCompletionTime',tableSolvingTime',netControlEffort',tableNetProgress',tableNetProgress'./tableCompletionTime','VariableNames', {'Possible Combinations', 'Task Completion Time', 'Solving Time', 'Net Control Effort', 'Net Progress','Progress/Task Completion Time'});
else
     resultsAnalysisTable = table(tableCompletionTime',tableSolvingTime',netControlEffort',tableNetProgress',tableNetProgress'./tableCompletionTime','VariableNames', { 'Task Completion Time', 'Solving Time', 'Net Control Effort', 'Net Progress','Progress/Task Completion Time'});
end

%%
fullfig
set(gca,'ydir','reverse')
backgroundImage = imread(backgroundImagePath);
backgroundImage = backgroundImage(:, :, :);
h = image([.0 .0], [.0 .0], backgroundImage);
uistack(h,'bottom')
hold on

statesplayers{1} = frenet2global(players{1}.referencePath,[76 0 0 0 0 0]);
statesplayers{2} = frenet2global(players{2}.referencePath,[78 0 0 0 0 0]);
statesplayers{3} = frenet2global(players{3}.referencePath,[105 0 0 0 0 0]);
statesplayers{4} = frenet2global(players{4}.referencePath,[74 0 0 0 0 0]);

for i = 1:1:4
        plot(players{i}.pathInfo.centerPoints(:,1)/meterPerPixel, players{i}.pathInfo.centerPoints(:,2)/meterPerPixel, 'Color',players{i}.params.col)
    plot(players{i}.pathInfo.lowerBound(:,1)/meterPerPixel, players{i}.pathInfo.lowerBound(:,2)/meterPerPixel, 'Color',players{i}.params.col)
    plot(players{i}.pathInfo.upperBound(:,1)/meterPerPixel, players{i}.pathInfo.upperBound(:,2)/meterPerPixel, 'Color',players{i}.params.col)

players{i}.drawPlayer([statesplayers{i}(1)/meterPerPixel,statesplayers{i}(2)/meterPerPixel,statesplayers{i}(3)]);
end
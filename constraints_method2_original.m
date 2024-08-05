import casadi.*

% Allocating the variables
bigM = 10000;

homotopy_setup = [];

homotopy_setup.vars.h = cell(1,size(playersPairs,1));

homotopy_setup.vars.sigma.sigma_AF = cell(1,size(playersPairs,1));
homotopy_setup.vars.sigma.sigma_BE = cell(1,size(playersPairs,1));
homotopy_setup.vars.sigma.sigma_CD = cell(1,size(playersPairs,1));

homotopy_setup.vars.sigma.sigma_AFp = cell(1,size(playersPairs,1));
homotopy_setup.vars.sigma.sigma_BEp = cell(1,size(playersPairs,1));
homotopy_setup.vars.sigma.sigma_CDp = cell(1,size(playersPairs,1));


homotopy_setup.vars.slack.slack_AF = cell(1,size(playersPairs,1));
% homotopy_setup.vars.slack.slack_BE = cell(1,size(playersPairs,1));
homotopy_setup.vars.slack.slack_CD = cell(1,size(playersPairs,1));

homotopy_setup.constraints.general.constraints =  cell(1,size(playersPairs,1));
homotopy_setup.constraints.general.ub =  cell(1,size(playersPairs,1));
homotopy_setup.constraints.general.lb =  cell(1,size(playersPairs,1));

homotopy_setup.constraints.vars.constraints =  cell(1,size(playersPairs,1));
homotopy_setup.constraints.vars.ub =  cell(1,size(playersPairs,1));
homotopy_setup.constraints.vars.lb =  cell(1,size(playersPairs,1));

for i = 1:1:size(playersPairs,1)

    tempconstraints = [];
    tempub = [];
    templb = [];

    varsconstraints = [];
    varsub = [];
    varslb = [];


    if intersectionBetweenEnvelopesOfPlayers(i) > 0
        homotopy_setup.vars.h{i} = binvar(1,1);

        homotopy_setup.vars.sigma.sigma_AF{i} = binvar(players{playersPairs(i,2)}.opt.params.N+1,1);

        homotopy_setup.vars.sigma.sigma_BE{i} = binvar(players{playersPairs(i,2)}.opt.params.N+1,1);
        
        homotopy_setup.vars.sigma.sigma_CD{i} = binvar(players{playersPairs(i,2)}.opt.params.N+1,1);

        homotopy_setup.vars.slack.slack_AF{i} = sdpvar(players{playersPairs(i,2)}.opt.params.N+1,1);

        homotopy_setup.vars.slack.slack_CD{i} = sdpvar(players{playersPairs(i,2)}.opt.params.N+1,1);
        
        homotopy_setup.vars.sigma.sigma_AFp{i} = sdpvar(players{playersPairs(i,2)}.opt.params.N+1,1);

        homotopy_setup.vars.sigma.sigma_BEp{i} = sdpvar(players{playersPairs(i,2)}.opt.params.N+1,1);
        
        homotopy_setup.vars.sigma.sigma_CDp{i} = sdpvar(players{playersPairs(i,2)}.opt.params.N+1,1);


        h = homotopy_setup.vars.h{i};
        x_1 = players{playersPairs(i,1)}.opt.vars.x;
        u_1 = players{playersPairs(i,1)}.opt.vars.u;

        x_2 = players{playersPairs(i,2)}.opt.vars.x;
        u_2 = players{playersPairs(i,2)}.opt.vars.u;

        sigma_af = homotopy_setup.vars.sigma.sigma_AF{i};
        sigma_be = homotopy_setup.vars.sigma.sigma_BE{i};
        sigma_cd = homotopy_setup.vars.sigma.sigma_CD{i};
        
        
        sigma_afp = homotopy_setup.vars.sigma.sigma_AFp{i};
        sigma_bep = homotopy_setup.vars.sigma.sigma_BEp{i};
        sigma_cdp = homotopy_setup.vars.sigma.sigma_CDp{i};

        slack_af = homotopy_setup.vars.slack.slack_AF{i};
        slack_cd = homotopy_setup.vars.slack.slack_CD{i};

        for k = 1:1:players{playersPairs(i,1)}.opt.params.N+1
            
            %%%% CONSTRAINTS EQ.17 %%%%%

            % 17a
            tempconstraints = [tempconstraints; x_2(k,1) + bigM*sigma_af(k) - bigM*h*sigma_af(k)];
            tempub = [tempub; s2_l_entry{i} + bigM];
            templb = [templb; -inf];

            % 17b
            tempconstraints = [tempconstraints; x_2(k,1) - x_1(k,1) + bigM*sigma_be(k) - bigM*h*sigma_be(k)];
            tempub = [tempub; bigM + (s2_l_entry{i} - s1_h_entry{i})];
            templb = [templb; -inf];
            
            % 17c
            tempconstraints = [tempconstraints; x_1(k,1) - bigM*sigma_cd(k) + bigM*h*sigma_cd(k)];
            tempub = [tempub; inf];
            templb = [templb; s1_h_exit{i} - bigM];
                        
            % 17d
            tempconstraints = [tempconstraints; x_1(k,1) + bigM*sigma_cd(k)*h];
            tempub = [tempub; s1_l_entry{i} + bigM];
            templb = [templb; -inf];
 
            % 17e
            tempconstraints = [tempconstraints; x_2(k,1) - x_1(k,1) - bigM*sigma_be(k)*h];
            tempub = [tempub; inf];
            templb = [templb; (s2_h_entry{i} - s1_l_entry{i}) - bigM];

            % 17f
            tempconstraints = [tempconstraints; x_2(k,1) - bigM*sigma_af(k)*h];
            tempub = [tempub; inf];
            templb = [templb; s2_h_exit{i} - bigM];

            tempconstraints = [tempconstraints;  sigma_af(k) + sigma_be(k) + sigma_cd(k)];
            tempub=[tempub; 1];
            templb=[templb;1];
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            tempconstraints = [tempconstraints; slack_af(k) - bigM*h];
            tempub = [tempub; 1];
            templb = [templb;-inf];

            tempconstraints = [tempconstraints; -slack_af(k) - bigM*h];
            tempub = [tempub; 0];
            templb = [templb;-inf];                

            tempconstraints = [tempconstraints; slack_af(k) + bigM*h];
            tempub = [tempub; bigM];
            templb = [templb;-inf];
            
            tempconstraints = [tempconstraints; -slack_af(k) + bigM*h];
            tempub = [tempub; 1 + bigM];
            templb = [templb;-inf];
            
            
            tempconstraints = [tempconstraints; slack_cd(k) - bigM*h];
            tempub = [tempub; 1];
            templb = [templb;-inf];

            tempconstraints = [tempconstraints; -slack_cd(k) - bigM*h];
            tempub = [tempub; 0];
            templb = [templb;-inf];                

            tempconstraints = [tempconstraints; slack_cd(k) + bigM*h];
            tempub = [tempub; bigM];
            templb = [templb;-inf];
            
            tempconstraints = [tempconstraints; -slack_cd(k) + bigM*h];
            tempub = [tempub; 1 + bigM];
            templb = [templb;-inf];

           
            %%% CONSTRAINTS EQ.18 %%%%%
            if k < players{playersPairs(i,1)}.opt.params.N+1
                tempconstraints = [tempconstraints; sigma_af(k) - sigma_af(k+1) - slack_af(k+1)];
                tempub = [tempub; 0];
                templb = [templb; 0];
    
                tempconstraints = [tempconstraints; sigma_cd(k+1) - sigma_cd(k) -  slack_cd(k+1)];
                tempub = [tempub; 0];
                templb = [templb; 0];

            end

        end

        discflags = [];
        varsconstraints = [varsconstraints; h];
        varsub = [1];
        varslb = [0];

        % Constraints on sigmas
        % 
%%%%%%%%%%%%%%%%%%%

        varsconstraints = [varsconstraints; sigma_af(:)];
        varsub = [varsub; ones(size(sigma_af(:)))];
        varslb = [varslb; zeros(size(sigma_af(:)))];

        varsconstraints = [varsconstraints; sigma_be(:)];
        varsub = [varsub; ones(size(sigma_be(:)))];
        varslb = [varslb; zeros(size(sigma_be(:)))];

        varsconstraints = [varsconstraints; sigma_cd(:)];
        varsub = [varsub; ones(size(sigma_cd(:)))];
        varslb = [varslb; zeros(size(sigma_cd(:)))];

        varsconstraints = [varsconstraints; slack_af(:)];
        varsub = [varsub; ones(size(slack_af(:)))];
        varslb = [varslb; -ones(size(slack_af(:)))];

%             varsconstraints = [varsconstraints; slack_be(k)];
%             varsub = [varsub; 1];
%             varslb = [varslb; -1];

        varsconstraints = [varsconstraints; slack_cd(:)];
        varsub = [varsub; ones(size(slack_cd(:)))];
        varslb = [varslb; -ones(size(slack_cd(:)))];

        homotopy_setup.constraints.general.constraints{i} = tempconstraints;
        homotopy_setup.constraints.general.ub{i} = tempub;
        homotopy_setup.constraints.general.lb{i} = templb;

        homotopy_setup.constraints.vars.constraints{i} = varsconstraints;
        homotopy_setup.constraints.vars.ub{i} = varsub;
        homotopy_setup.constraints.vars.lb{i} = varslb;

    end
    

end
%% Now, we impose scenario-specific constraints
n = length(players);
maxLength = 0;

for i = 1:1:n
    maxLength = max(maxLength, players{i}.referencePath.PathLength);
end

figure
ylim([0,n+1])
xlim([0,maxLength])
yticks([1 2 3  4])
ytickformat('Player %.0f')
ytickangle(45)
hold on
for i = 1:1:n
    yline(i,'LineWidth',2,'Color','k')
end

hVec = {0,1,[],0,0,[]};

for j = 1:1:length(intersectionBetweenEnvelopesOfPlayers)
    if intersectionBetweenEnvelopesOfPlayers(j) >0
        player1Index = playersPairs(j,1);
        player2Index = playersPairs(j,2);

        scatter(s1_l_entry{j},player1Index,'k','filled');
        text(s1_l_entry{j},player1Index, strcat('EP_',num2str(player2Index)),'VerticalAlignment','top');

        scatter(s2_l_entry{j},player2Index,'k','filled');
        text(s2_l_entry{j},player2Index, strcat('EP_',num2str(player1Index)),'VerticalAlignment','top');

         scatter(s1_h_exit{j},player1Index,'k','filled');
         text(s1_h_exit{j},player1Index, strcat('XP_',num2str(player2Index)),'VerticalAlignment','top');

         scatter(s2_h_exit{j},player2Index,'k','filled');
         text(s2_h_exit{j},player2Index, strcat('XP_',num2str(player1Index)),'VerticalAlignment','top');
    end
end

geometricPrecedenceList = cell(n,9);

% 1 displays Player Number
% 2 displays Entry point to different intersection areas along the player's
% path;
% 3 displays the indices of other players (sorted)
% 4 displays the indices in pairsPlayers corresponding to the players
% 5 displays whether we are interested in the normal or inverted
% constraints
% 6 displays Exit point to different intersection areas along the player's
% path;
% 7 displays the indices of other players (sorted)
% 8 displays the indices in pairsPlayers corresponding to the players
% 9 displays whether we are interested in the normal or inverted
% constraints


for i = 1:1:n
    geometricPrecedenceList{i} = i;
end

for j = 1:1:length(intersectionBetweenEnvelopesOfPlayers)
    if intersectionBetweenEnvelopesOfPlayers(j) >0
        
        player1Index = playersPairs(j,1);
        player2Index = playersPairs(j,2);
        
        geometricPrecedenceList{player1Index,3} = [geometricPrecedenceList{player1Index,3}, player2Index];
        geometricPrecedenceList{player2Index,3} = [geometricPrecedenceList{player2Index,3}, player1Index];

        geometricPrecedenceList{player1Index,4} = [geometricPrecedenceList{player1Index,4}, j];
        geometricPrecedenceList{player2Index,4} = [geometricPrecedenceList{player2Index,4}, j];
        
        geometricPrecedenceList{player1Index,7} = [geometricPrecedenceList{player1Index,7}, player2Index];
        geometricPrecedenceList{player2Index,7} = [geometricPrecedenceList{player2Index,7}, player1Index];

        geometricPrecedenceList{player1Index,8} = [geometricPrecedenceList{player1Index,8}, j];
        geometricPrecedenceList{player2Index,8} = [geometricPrecedenceList{player2Index,8}, j];


        s1e_p1 = s1_l_entry{j};
        s2e_p2 = s2_l_entry{j};
        
        s1x_p1 = s1_h_exit{j};
        s2x_p2 = s2_h_exit{j};

        geometricPrecedenceList{player1Index,2} = [geometricPrecedenceList{player1Index,2}, s1e_p1];
        geometricPrecedenceList{player2Index,2} = [geometricPrecedenceList{player2Index,2}, s2e_p2];

        geometricPrecedenceList{player1Index,6} = [geometricPrecedenceList{player1Index,6}, s1x_p1];
        geometricPrecedenceList{player2Index,6} = [geometricPrecedenceList{player2Index,6}, s2x_p2];
    end

end

for j = 1:1:length(players)
    [B,I] = sort(geometricPrecedenceList{j,2});
    geometricPrecedenceList{j,2} = B;
    geometricPrecedenceList{j,3} = geometricPrecedenceList{j,3}(I);
    geometricPrecedenceList{j,4} = geometricPrecedenceList{j,4}(I);    
    
    [B,I] = sort(geometricPrecedenceList{j,6});
    geometricPrecedenceList{j,6} = B;
    geometricPrecedenceList{j,7} = geometricPrecedenceList{j,7}(I);
    geometricPrecedenceList{j,8} = geometricPrecedenceList{j,8}(I);    
end

for j = 1:1:length(players)

    % Part for Entry Points
    playerIndex = geometricPrecedenceList{j,1};
    otherPlayersIndices = geometricPrecedenceList{j,3};

    temp = [];

    for i = 1:1:length(otherPlayersIndices)

        if otherPlayersIndices(i) > playerIndex
            temp(i) =1;
        else
            temp(i) = -1;
        end
    end

    geometricPrecedenceList{j,5} = temp;

    % Part for Exit Points
    otherPlayersIndices = geometricPrecedenceList{j,7};

    temp = [];

    for i = 1:1:length(otherPlayersIndices)

        if otherPlayersIndices(i) > playerIndex
            temp(i) =1;
        else
            temp(i) = -1;
        end
    end

    geometricPrecedenceList{j,9} = temp;

end

%%
scenario_specific_vars = [];
scenario_specific_lb = [];
scenario_specific_ub = [];

for j = 1:1:length(players)
    temp_vars = [];
    temp_lb = [];
    temp_ub = [];
    player1Index = geometricPrecedenceList{j,1};
    sprintf('Now we do for Player %d',player1Index)
    otherPlayersIndices_Entry = geometricPrecedenceList{j,3};
    pairsPlayersIndices_Entry = geometricPrecedenceList{j,4};
    orderOfSigmas_Entry = geometricPrecedenceList{j,5};

    otherPlayersIndices_Exit = geometricPrecedenceList{j,7};
    pairsPlayersIndices_Exit = geometricPrecedenceList{j,8};
    orderOfSigmas_Exit = geometricPrecedenceList{j,9};

    if length(otherPlayersIndices_Entry) > 1
    sprintf('Entry Points')
        % First Type of Constraint -> Moving from left to right of
        % s_player1Index

        temp_otherPlayersIndices = otherPlayersIndices_Entry;
        temp_pairsPlayersIndices = pairsPlayersIndices_Entry;
        temp_orderOfSigmas = orderOfSigmas_Entry;

        while length(temp_otherPlayersIndices)>1

            player2Index = temp_otherPlayersIndices(1);
            sprintf('Now we do between Player %d and Player %d',player1Index, player2Index)

            if temp_orderOfSigmas(1) == 1
                sprintf('This is the usual order of sigmas, we select D')
                sigmaOfInterest =  homotopy_setup.vars.sigma.sigma_CD{temp_pairsPlayersIndices(1)};
                % sprintf('Just to double check, the sigmaOfInterest corresponds to indices of Players %d and %d', playersPairs(temp_pairsPlayersIndices(1),1), playersPairs(temp_pairsPlayersIndices(1),2))
            else
                sprintf('This is the reverse order of sigmas, we select A')
                sigmaOfInterest =  homotopy_setup.vars.sigma.sigma_AF{temp_pairsPlayersIndices(1)};
                % sprintf('Just to double check, the sigmaOfInterest corresponds to indices of Players %d and %d', playersPairs(temp_pairsPlayersIndices(1),1), playersPairs(temp_pairsPlayersIndices(1),2))
            end
            
            for i = 2:1:length(temp_otherPlayersIndices)

                sprintf('Now we apply the precendence constraints indices of Players %d and %d', playersPairs(temp_pairsPlayersIndices(i),1), playersPairs(temp_pairsPlayersIndices(i),2))
                h_varx = homotopy_setup.vars.h{temp_pairsPlayersIndices(i)};

                    sigmaOfInterest_2_1 =  homotopy_setup.vars.sigma.sigma_CD{temp_pairsPlayersIndices(i)};
                    sigmaOfInterest_2_2 =  homotopy_setup.vars.sigma.sigma_AF{temp_pairsPlayersIndices(i)};

%                 for k = 1:1:length(sigmaOfInterest_2_1)
%                     temp_vars = [temp_vars; sigmaOfInterest_2_1(k) + sigmaOfInterest_2_3(k)];
%                     temp_lb = [temp_lb; 1 - bigM*(1-sigmaOfInterest(k)) - bigM*(1-h_varx)];
%                     temp_ub = [temp_ub; 1 + bigM*(1-sigmaOfInterest(k)) + bigM*(1-h_varx)];
% 
%                     temp_vars = [temp_vars; sigmaOfInterest_2_2(k)];
%                     temp_lb = [temp_lb; 1 - bigM*(1-sigmaOfInterest(k)) - bigM*(h_varx)];
%                     temp_ub = [temp_ub; 1 + bigM*(1-sigmaOfInterest(k)) + bigM*(h_varx)];                   
%                 end

                for k = 1:1:length(sigmaOfInterest_2_1)
                    temp_vars = [temp_vars; sigmaOfInterest_2_1(k) +  sigmaOfInterest_2_2(k)];
                    temp_lb = [temp_lb; 0.9 - bigM*(1-sigmaOfInterest(k))];
                    temp_ub = [temp_ub; 1.1 + bigM*(1-sigmaOfInterest(k))];              
                end

            end
            
            temp_otherPlayersIndices(1) = [];
            temp_pairsPlayersIndices(1) = [];
            temp_orderOfSigmas(1) = [];

        end


        % Second Type of Constraint -> Moving from right to left of
        % s_player1Index

        temp_otherPlayersIndices = flip(otherPlayersIndices_Exit);
        temp_pairsPlayersIndices = flip(pairsPlayersIndices_Exit);
        temp_orderOfSigmas = flip(orderOfSigmas_Exit);
        sprintf('Exit Points')
        while length(temp_otherPlayersIndices)>1

            player2Index = temp_otherPlayersIndices(1);
            sprintf('Now we do between Player %d and Player %d',player1Index, player2Index)
            if temp_orderOfSigmas(1) == 1
                sprintf('This is the usual order of sigmas, we select C')
                sigmaOfInterest =  homotopy_setup.vars.sigma.sigma_C{temp_pairsPlayersIndices(1)};
            else
                sprintf('This is the usual order of sigmas, we select F')
                sigmaOfInterest =  homotopy_setup.vars.sigma.sigma_F{temp_pairsPlayersIndices(1)};
            end
            
            for i = 2:1:length(temp_otherPlayersIndices)
                sprintf('Now we apply the precendence constraints indices of Players %d and %d', playersPairs(temp_pairsPlayersIndices(i),1), playersPairs(temp_pairsPlayersIndices(i),2))
                h_varx = homotopy_setup.vars.h{temp_pairsPlayersIndices(i)};
                sigmaOfInterest_2_1 =  homotopy_setup.vars.sigma.sigma_AF{temp_pairsPlayersIndices(i)};
                sigmaOfInterest_2_2 =  homotopy_setup.vars.sigma.sigma_CD{temp_pairsPlayersIndices(i)};
% 
%                 for k = 1:1:length(sigmaOfInterest_2_1)
%                     temp_vars = [temp_vars; sigmaOfInterest_2_1(k)];
%                     temp_lb = [temp_lb; 1 - bigM*(1-sigmaOfInterest(k)) - bigM*(1-h_varx)];
%                     temp_ub = [temp_ub; 1 + bigM*(1-sigmaOfInterest(k)) + bigM*(1-h_varx)];
% 
%                     temp_vars = [temp_vars; sigmaOfInterest_2_2(k) + sigmaOfInterest_2_3(k)];
%                     temp_lb = [temp_lb; 1 - bigM*(1-sigmaOfInterest(k)) - bigM*(h_varx)];
%                     temp_ub = [temp_ub; 1 + bigM*(1-sigmaOfInterest(k)) + bigM*(h_varx)];                   
%                 end

                for k = 1:1:length(sigmaOfInterest_2_1)
                    temp_vars = [temp_vars; sigmaOfInterest_2_1(k)+sigmaOfInterest_2_2(k)];
                    temp_lb = [temp_lb; 0.9 - bigM*(1-sigmaOfInterest(k))];
                    temp_ub = [temp_ub; 1.1 + bigM*(1-sigmaOfInterest(k))];
                
                end

            end
            
            temp_otherPlayersIndices(1) = [];
            temp_pairsPlayersIndices(1) = [];
            temp_orderOfSigmas(1) = [];

        end
        

    end
    
    scenario_specific_vars = [scenario_specific_vars; temp_vars];
    scenario_specific_lb = [scenario_specific_lb; temp_lb];
    scenario_specific_ub = [scenario_specific_ub; temp_ub];

end
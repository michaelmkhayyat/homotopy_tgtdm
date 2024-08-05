intersectionBetweenEnvelopesOfPlayers;
playersPairs;

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
        text(s1_l_entry{j},player1Index, strcat('E_L P_',num2str(player2Index)),'VerticalAlignment','top');

        scatter(s2_l_entry{j},player2Index,'k','filled');
        text(s2_l_entry{j},player2Index, strcat('E_L P_',num2str(player1Index)),'VerticalAlignment','top');

%         scatter(s1_h_entry{j},player1Index,'k','filled');
%         text(s1_h_entry{j},player1Index, strcat('E_H P_',num2str(player2Index)),'VerticalAlignment','top');
% 
%         scatter(s2_h_entry{j},player2Index,'k','filled');
%         text(s2_h_entry{j},player2Index, strcat('E_H P_',num2str(player1Index)),'VerticalAlignment','top');
% 
         scatter(s1_h_exit{j},player1Index,'k','filled');
         text(s1_h_exit{j},player1Index, strcat('X_H P_',num2str(player2Index)),'VerticalAlignment','top');

         scatter(s2_h_exit{j},player2Index,'k','filled');
         text(s2_h_exit{j},player2Index, strcat('X_H P_',num2str(player1Index)),'VerticalAlignment','top');
%          
%          scatter(s1_l_exit{j},player1Index,'k','filled');
%          text(s1_l_exit{j},player1Index, strcat('X_L P_',num2str(player2Index)),'VerticalAlignment','top');
% 
%          scatter(s2_l_exit{j},player2Index,'k','filled');
%          text(s2_l_exit{j},player2Index, strcat('X_L P_',num2str(player1Index)),'VerticalAlignment','top');
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
% How to determine if a a homotopy class is feasible or not


nonZeroElements = nnz(intersectionBetweenEnvelopesOfPlayers);
possibleCombinations = ff2n(nonZeroElements);

hvec = {};


% We need to add some order to each vector. We represent entry point by E,
% and exit point by X
orderedList = cell(n,4);
for i = 1:1:length(players)
    player1Index = geometricPrecedenceList{i,1};
    EntryPointsWithPlayers = geometricPrecedenceList{i,2};
    EntryPointPlayersIndices = geometricPrecedenceList{i,3};
    EntryPointsIndicesPairsPlayers = geometricPrecedenceList{i,4};
    
    ExitPointsWithPlayers = geometricPrecedenceList{i,6};
    ExitPointPlayersIndices = geometricPrecedenceList{i,7};
    ExitPointsIndicesPairsPlayers = geometricPrecedenceList{i,8};

    EX_str =[repmat(["e"],1,length(EntryPointsWithPlayers)),repmat(["x"],1,length(EntryPointsWithPlayers))];
    EX_vec = [EntryPointsWithPlayers,ExitPointsWithPlayers];
    EX_ind = [EntryPointPlayersIndices, ExitPointPlayersIndices];
    EX_pp = [EntryPointsIndicesPairsPlayers, ExitPointsIndicesPairsPlayers];
    
    [B,I] = sort(EX_vec);

    EX_vec = B;
    EX_str = EX_str(I);
    EX_ind = EX_ind(I);
    EX_pp = EX_pp(I);

    orderedList{i,1} = player1Index;
    orderedList{i,2} = EX_ind;
    orderedList{i,3} = EX_str;
    orderedList{i,4} = EX_pp;
    orderedList{i,5} = B;
end
%%

overLapMatrix = eye(n);
copyOrderedList = orderedList;

cc = copyOrderedList;
for kk=9%1:1:length(possibleCombinations)
    orderedList = copyOrderedList;
    homotopyClassToTest = possibleCombinations(kk,:);
    j = 1;
    
    for i = 1:1:length(intersectionBetweenEnvelopesOfPlayers)
    
        if intersectionBetweenEnvelopesOfPlayers(i) >0
            hvec{i} = homotopyClassToTest(j);
            j = j+1;
        else
            hvec{i} = [];
        end
    
    end
    
    progressMade = ones(n,1);
    emptyFlag = zeros(n,1);
    
    playerIndex = 1;
    
    while any(progressMade) && ~all(emptyFlag)
        
        progressMade = zeros(n,1);
        
        for i = 1:1:n
            
            if emptyFlag(i)
                continue;
            end
            
            player1Index = orderedList{i,1};
            otherPlayersIndices = orderedList{i,2};
            entryexitStrings = orderedList{i,3};
            pairsPlayersIndices = orderedList{i,4};
            lengthIndices = length(otherPlayersIndices);
    
            for j = 1:1:lengthIndices
    
                cond1 = orderedList{i,2}(1) > player1Index && orderedList{i,3}(1) == "e" && hvec{orderedList{i,4}(1)} == 0;
                cond2 = orderedList{i,2}(1) < player1Index && orderedList{i,3}(1) == "e" && hvec{orderedList{i,4}(1)} == 1;
                cond3 = orderedList{i,2}(1) > player1Index && orderedList{i,3}(1) == "x" && hvec{orderedList{i,4}(1)} == 0;
                cond4 = orderedList{i,2}(1) < player1Index && orderedList{i,3}(1) == "x" && hvec{orderedList{i,4}(1)} == 1;
    
                if cond1||cond2||cond3||cond4
    
                    progressMade(i) = 1;
                    
                    idx = find(orderedList{orderedList{i,2}(1),2}==player1Index,1,"first");
    
                    orderedList{orderedList{i,2}(1),2}(idx) = [];
                    orderedList{orderedList{i,2}(1),3}(idx) = [];
                    orderedList{orderedList{i,2}(1),4}(idx) = [];


                    if cond1 || cond2
                       overLapMatrix(i,orderedList{i,2}(1)) = 1;
                    end

                    if cond3 || cond4
                        overLapMatrix(i,orderedList{i,2}(1)) = 0; 
                    end
    
                    orderedList{i,2}(1) = [];
                    orderedList{i,3}(1) = [];
                    orderedList{i,4}(1) = [];
        
                else
                    break;
                end
    
            end
            
            for k = 1:1:n
                
                if isempty(orderedList{k,2})
                    emptyFlag(k) = 1;
                
                end
    
            end
                
        end
        
        if all(emptyFlag) || all(~progressMade)
            break;
        end
    
    
    end
    
    if all(emptyFlag)
        disp('Feasible!')
    else
        disp('Infeasible!')
    end
end

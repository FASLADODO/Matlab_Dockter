% Simulation based on Cracker Barrels 15 peg logic game
% Rod Dockter
% December 2014

%% initialize
numPegs = 15;
numRows = 5;

% {row}.column
pegData{1}.plotPos = [1];
pegData{2}.plotPos = [1,2];
pegData{3}.plotPos = [1,2,3];
pegData{4}.plotPos = [1,2,3,4];
pegData{5}.plotPos = [1,2,3,4,5];
% only for plotting in triangle
pegData{1}.plotOffset = 1;
pegData{2}.plotOffset = 1.5;
pegData{3}.plotOffset = 2;
pegData{4}.plotOffset = 2.5;
pegData{5}.plotOffset = 3;
% 1 = peg, 0 = hole
pegData{1}.pegIO = [1];
pegData{2}.pegIO = [1,1];
pegData{3}.pegIO = [1,1,1];
pegData{4}.pegIO = [1,1,1,1];
pegData{5}.pegIO = [1,1,1,1,1];
% move lookup table [jumpx,jumpy,openx,openy]
pegData{1}.moveLUT{1}.t = [2,1,3,1;
                        2,2,3,3];
pegData{2}.moveLUT{1}.t = [3,1,4,1;
                        3,2,4,3];
pegData{2}.moveLUT{2}.t = [3,2,4,2;
                        3,3,4,4];
pegData{3}.moveLUT{1}.t = [2,1,1,1;
                        3,2,3,3;
                        4,1,5,1;
                        4,2,5,3];
pegData{3}.moveLUT{2}.t = [4,2,5,2;
                        4,3,5,4];
pegData{3}.moveLUT{3}.t = [2,2,1,1;
                        3,2,3,1;
                        4,3,5,3;
                        4,4,5,5];
pegData{4}.moveLUT{1}.t = [3,1,2,1;
                        4,2,4,3];
pegData{4}.moveLUT{2}.t = [3,2,2,2;
                        4,3,4,4];
pegData{4}.moveLUT{3}.t = [3,2,2,1;
                        4,2,4,1];
pegData{4}.moveLUT{4}.t = [3,3,2,2;
                        4,3,4,2];
pegData{5}.moveLUT{1}.t = [4,1,3,1;
                        5,2,5,3];
pegData{5}.moveLUT{2}.t = [4,2,3,2;
                        5,3,5,4];
pegData{5}.moveLUT{3}.t = [4,2,3,1;
                        4,3,3,3;
                        5,2,5,1;
                        5,4,5,5];
pegData{5}.moveLUT{4}.t = [4,3,3,2;
                        5,3,5,2];
pegData{5}.moveLUT{5}.t = [4,4,3,3;
                        5,4,5,3];
                    
%% run sim

% user entered start peg
randRow = 2;
randCol = 2;

% reset stuff
clear bestRun bestStart
minLeft = 15;
numberTrials = 0;

tic
while 1
    clear allMoves
    moveIndex = 1;
    
    % reset pegs all in
    pegData{1}.pegIO = [1];
    pegData{2}.pegIO = [1,1];
    pegData{3}.pegIO = [1,1,1];
    pegData{4}.pegIO = [1,1,1,1];
    pegData{5}.pegIO = [1,1,1,1,1];

    %pegs left
    numLeft = numPegs;

    %initially remove one peg
%     randRow = randi([1 5],1,1);
%     randCol = randi([1 randRow],1,1);
    pegData{randRow}.pegIO(randCol) = 0;
    numLeft = numLeft - 1;
    startPeg = [randRow,randCol];

    while 1

        %reset currently available moves
        clear availableMoves
        indA = 0;

        %get all currently available moves
        for rr=1:5
           for cc = 1:rr
               [moveR, moveC] = size(pegData{rr}.moveLUT{cc}.t);
               for oo = 1:moveR
                  % check for all move options with a peg then an opening via the
                  % lookup table
                  if ( pegData{rr}.pegIO(cc) == 1) %check that current peg is in
                      if( pegData{ pegData{rr}.moveLUT{cc}.t(oo,1) }.pegIO( pegData{rr}.moveLUT{cc}.t(oo,2) ) == 1 ) % check that hop peg is in
                          if (pegData{ pegData{rr}.moveLUT{cc}.t(oo,3) }.pegIO( pegData{rr}.moveLUT{cc}.t(oo,4) ) == 0 ) %check that move space is open
                              indA = indA + 1;
                              availableMoves{indA}.moves = [rr,cc,pegData{rr}.moveLUT{cc}.t(oo,:) ]; 
                          end
                      end
                  end
               end
           end
        end

        %if no moves left break
        if (indA == 0)
           break; 
        end

        %choose random move
        randMove = randi([1 length(availableMoves)],1,1);
        %changle relevant peg info
        pegData{availableMoves{randMove}.moves(1)}.pegIO(availableMoves{randMove}.moves(2)) = 0;
        pegData{availableMoves{randMove}.moves(3)}.pegIO(availableMoves{randMove}.moves(4)) = 0;
        pegData{availableMoves{randMove}.moves(5)}.pegIO(availableMoves{randMove}.moves(6)) = 1;

        %decrement amount left
        numLeft = numLeft - 1;

        allMoves(moveIndex,:) = availableMoves{randMove}.moves;
        moveIndex = moveIndex + 1;


    end
    
    % keep trying until a run with 1 left
    if (numLeft < minLeft)
        minLeft = numLeft;
        bestRun = allMoves;
        bestStart = startPeg;
    end
    if(minLeft == 1)
       break; 
    end
    
    %count up the number of trials
    numberTrials = numberTrials + 1;
    
end
toc

numberTrials

%% Plot best run

% reset pegs all in
pegData{1}.pegIO = [1];
pegData{2}.pegIO = [1,1];
pegData{3}.pegIO = [1,1,1];
pegData{4}.pegIO = [1,1,1,1];
pegData{5}.pegIO = [1,1,1,1,1];

% length best run
[lBest,garb] = size(bestRun);

%best start peg
pegData{startPeg(1)}.pegIO(startPeg(2)) = 0;
    
% Plot initial peg
figure(1)
for ii = 1:numRows
    for jj = 1:ii
        if ii == startPeg(1) && jj == startPeg(2)
            plot(pegData{ii}.plotPos(jj)-pegData{ii}.plotOffset,1-ii,'o','MarkerEdgeColor',[1,0,0],'MarkerSize',10);
        elseif pegData{ii}.pegIO(jj) == 1
            plot(pegData{ii}.plotPos(jj)-pegData{ii}.plotOffset,1-ii,'x','MarkerEdgeColor',[0,0,0],'MarkerSize',10);
        else
            plot(pegData{ii}.plotPos(jj)-pegData{ii}.plotOffset,1-ii,'o','MarkerEdgeColor',[0,0,0],'MarkerSize',10);
        end
        hold on
    end
end
hold off

title('Peg Position (x=peg, o=hole)')
axis([-3,3,-5,1]);  

pause(2)


for kk = 1:lBest
    % reset in out
    pegData{bestRun(kk,1)}.pegIO(bestRun(kk,2)) = 0;
    pegData{bestRun(kk,3)}.pegIO(bestRun(kk,4)) = 0;
    pegData{bestRun(kk,5)}.pegIO(bestRun(kk,6)) = 1;
    
    figure(1)
    for ii = 1:numRows
        for jj = 1:ii
            %plot each move
            if ii == bestRun(kk,1) && jj == bestRun(kk,2) % start peg
                plot(pegData{ii}.plotPos(jj)-pegData{ii}.plotOffset,1-ii,'o','MarkerEdgeColor',[1,0,0],'MarkerSize',10);
            elseif ii == bestRun(kk,3) && jj == bestRun(kk,4) % jump peg
                plot(pegData{ii}.plotPos(jj)-pegData{ii}.plotOffset,1-ii,'o','MarkerEdgeColor',[0,1,0],'MarkerSize',10);
            elseif ii == bestRun(kk,5) && jj == bestRun(kk,6) % ending peg
                plot(pegData{ii}.plotPos(jj)-pegData{ii}.plotOffset,1-ii,'x','MarkerEdgeColor',[1,0,0],'MarkerSize',10);
            elseif pegData{ii}.pegIO(jj) == 1
                plot(pegData{ii}.plotPos(jj)-pegData{ii}.plotOffset,1-ii,'x','MarkerEdgeColor',[0,0,0],'MarkerSize',10);
            else
                plot(pegData{ii}.plotPos(jj)-pegData{ii}.plotOffset,1-ii,'o','MarkerEdgeColor',[0,0,0],'MarkerSize',10);
            end
            hold on
        end
    end
    hold off

    title('Peg Position (x=peg, o=hole)')
    axis([-3,3,-5,1]);
    
    pause(2)

end

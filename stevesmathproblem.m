%steves math problem

yearin = 75000;
interest = 0.03;
numyears = 10;

total = 0;
prevyear = 0;

stash = [];
% on principal + earnings
for yy = 1:numyears
    valtemp = prevyear + yearin;
    earnings = valtemp*interest
    total = earnings + valtemp
    prevyear = total;
    stash = [stash; total, earnings];
end

total

figure
plot(stash(:,1),'b');
hold on
plot(stash(:,2),'r');
hold off
%% only on principal

total = 0;
prevyear = 0;
sumearnings = 0;

stash = [];
for yy = 1:numyears
    valtemp = prevyear + yearin;
    earnings = valtemp*interest;
    sumearnings = sumearnings + earnings;
    total = sumearnings + valtemp;
    prevyear = valtemp;
    stash = [stash; total, earnings];
end
prevyear

total

figure
plot(stash(:,1),'b');
hold on
plot(stash(:,2),'r');
hold off

%% Other way

yy = numyears;
fact = 1:numyears;
principal = yearin*numyears
earnings = fact.*(yearin*interest)
total = principal + sum(earnings)
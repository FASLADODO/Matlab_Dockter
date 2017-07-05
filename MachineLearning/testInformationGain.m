%edibility data from banerjee
%Lecture 10 _ Decision Tree

keyskin = {'smooth';'scaly';'rough'};
Skin = [1,1,2,3,3,2,1,1,3,1,2,2,3,3]';
keyskin(Skin)

keycolor = {'pink';'purple';'orange'};
Color = [1,1,1,2,3,3,2,3,2,2,2,1,2,3]';

Thorny = [0,0,0,0,1,1,0,1,1,1,0,1,0,1]';
Flowering = [1,0,1,1,1,0,1,1,1,0,0,1,0,0]';
Edible = [1,1,0,0,0,0,1,0,0,0,0,0,1,1]';

DataAll = [Skin,Color,Thorny,Flowering,Edible]

%% check entropy 
%Compare
H_Edible = EntropySample(Edible)

H_Skin = EntropySample(Skin)

%% check conditional entropy 

H_Edible_Skin = ConditionalEntropy(Edible,Skin)

%% check information gain

%how much we know about Edibility given knowledge of Skin condition
IG_Edible_Skin = InformationGain(Edible,Skin)   %(0.2467)

%how much we know about Edibility given knowledge of Color
IG_Edible_Color = InformationGain(Edible,Color)   %(0.0292)

%how much we know about Edibility given knowledge of Thorniness
IG_Edible_Thorny = InformationGain(Edible,Thorny) %(0.1519)

%how much we know about Edibility given knowledge of Flowering
IG_Edible_Flowering = InformationGain(Edible,Flowering) %(0.0481)

%This means that Skin has the highest information gain, classify on this

%% Classify using skin

figure
gscatter(ones(length(Skin),1),Skin,Edible)
ylabel('Skin category')
title('Categorical data')

[Thresh,Direction,IG] = DecisionStumpCategorical(Skin,Edible)

EdibleEst = Direction((Skin == Thresh) + 1)';

corr = EdibleEst == Edible;
acc=  mean(corr)


%% Tennis data from Jure

colz.Day = 1;
colz.Outlook = 2;
colz.Temperature = 3;
colz.Humidity = 4;
colz.Wind = 5;
colz.Play = 6;

%create data set
Day = [1:14]';
keyoutlook = {'sunny';'overcast';'rain'};
Outlook = [0,0,1,2,2,2,1,0,0,2,0,1,1,2]';
keytemp = {'hot';'mild';'cool'};
Temperature = [0,0,0,1,2,2,2,1,2,1,1,1,0,1]';
keyhumidity = {'high';'normal'};
Humidity = [0,0,0,0,1,1,1,0,1,1,1,0,1,0]';
keywind = {'weak';'strong'};
Wind = [0,1,0,0,0,1,1,0,0,0,1,1,0,1]';
Play = [0,0,1,1,1,0,1,0,1,1,1,1,1,0]';

%combine all data
Data_Tennis = [Day,Outlook,Temperature,Humidity,Wind,Play]

%get subset where outlook was rain
S_Rain = Data_Tennis(Data_Tennis(:,colz.Outlook) == 2,:);

%compute new information gains
IG_S_outlook = InformationGain(S_Rain(:,colz.Play),S_Rain(:,colz.Temperature))
IG_S_humidity = InformationGain(S_Rain(:,colz.Play),S_Rain(:,colz.Humidity))
IG_S_wind = InformationGain(S_Rain(:,colz.Play),S_Rain(:,colz.Wind))









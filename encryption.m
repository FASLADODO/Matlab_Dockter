charstr = 'what the hell are you doing here?';



%%
str3 = ' ';

str3 = strcat(str3,num2str(45));
str3 = strcat(str3,num2str(69));
str3

%% encrypt

alphabet = 'abcdefghijklmnopqrstuvwxyz .!?';
alphsize = length(alphabet);
pattern = 1:alphsize;
pattern = pattern + 69;

messageencrypt=[];

for ii = 1:length(charstr)
    messageencrypt(ii) = pattern(strfind(alphabet,charstr(ii)));
end

messageencrypt

%% decrypt

messagedecrypt='';

for ii = 1:length(charstr)
    messagedecrypt(ii) = alphabet(find(pattern == messageencrypt(ii)));
end

messagedecrypt
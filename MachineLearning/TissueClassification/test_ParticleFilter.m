%% load up da data

load tissueData.mat

NonLinear_Fun_1 = @(Param,X) Param(1)*X(:,1) + Param(2)*X(:,2) + Param(3)*exp(-Param(4)*X(:,3) );

%% figure out true nonlinear params

for gg = 1:length(PatientData{1}.tissue)
    DataTrue{gg} = [];
    UTrue{gg} = [];
end

labelMaker = [1;2];

%get training and validation data
for pp = 1:length(PatientData)
    for gg = 1:length(PatientData{pp}.tissue)
        for ss = 1:length(PatientData{pp}.tissue{gg}.grasp)
            struct = PatientData{pp}.tissue{gg}.grasp{ss};

            DataTrue{gg} = [DataTrue{gg}; struct.State];
            UTrue{gg} = [UTrue{gg}; struct.Input];
        end
    end
    
end

%% get true nlin parameters

beta0{1} = [1, 0.02, 0.1, 0.01];
beta0{2} = [1, 0.9, 0.1, 0.4];

for gg = 1:length(PatientData{1}.tissue)
    ParamsTrue{gg} = nlinfit(DataTrue{gg},UTrue{gg},NonLinear_Fun_1,beta0{gg});
end

%% Particle Filter

NP = 40;
sdev_init = 0.5;
sdev = sdev_init;
nparams = length(ParamsTrue{1});
nclasses = length(PatientData{1}.tissue);
Particles = [];
Error = [];

%seed initial particles
for ccg = 1:nclasses
    ParamBase = ParamsTrue{ccg};
    for nn = 1:(NP/nclasses)
        ptemp = ParamBase + randn(1,nparams)*sdev;
        Particles = [Particles; ptemp];
    end
end

ParamStore = [];

% figure
% gg = 1;
% scatter(DataTrue{gg}(:,key.theta),DataTrue{gg}(:,key.thetadot),'ro')
% hold on
% gg = 2;
% scatter(DataTrue{gg}(:,key.theta),DataTrue{gg}(:,key.thetadot),'b+')
% hold off


ClassStore = [];

for pp = 1:length(PatientData)
    for gg = 1:length(PatientData{pp}.tissue)
        %new grasp reset particles
        for ccg = 1:nclasses
            ParamBase = ParamsTrue{ccg};
            for nn = 1:(NP/nclasses)
                ptemp = ParamBase + randn(1,nparams)*sdev;
                Particles = [Particles; ptemp];
            end
        end
        for ss = 1:length(PatientData{pp}.tissue{gg}.grasp)
            struct = PatientData{pp}.tissue{gg}.grasp{ss};
            timesteps = length(struct.Input);
            
            %go through all timesteps
            for tt = 1:timesteps;
                datanow = struct.State(tt,:);
                intputnow = struct.Input(tt,:);
                
                %try all parameters, figure out errors
                for nn = 1:NP
                    f_est = NonLinear_Fun_1(Particles(nn,:),datanow);
                    Error(nn,:) = abs(f_est-intputnow);
                end
                
                %resample
                [~,bestid] = min(Error,[],1);
                ParamBase = Particles(bestid,:);
                sdev = mean(Error)/5;
                thresherr = median(Error);
                
                %all new particles
                for nn = 1:(NP/nclasses)
                    ptemp = ParamBase + randn(1,nparams)*sdev;
                    if(Error(nn,:) > thresherr)
                        Particles(nn,:) = ptemp;
                    end
                end
                
                %classify based on best params
                for ccg = 1:nclasses
                    errParam(ccg) = mean(abs(ParamBase - ParamsTrue{ccg}));
                end
                [~,cid] = min(errParam);
                ClassStore = [ClassStore; cid,gg];
            end
            ParamStore = [ParamStore; ParamBase];
        end
    end
end


corr = ClassStore(:,1) == ClassStore(:,2);
acc = mean(corr)



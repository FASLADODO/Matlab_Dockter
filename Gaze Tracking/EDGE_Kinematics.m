% computes the Fwd Kinematics of EDGE for snsr.dataRAWfixed formant and
% units.  D is (n x 13) where n is number of data samples and D has format
%       [Time, Q1L Q2L, DL, Q3L, QgL, FgL, Q1R Q2R, DR, Q3R, QgR, FgR]
%   kVer: the kinematics version number:  -1 is default(EDGE2)
%       0=LabView Version; 1=EDGE1; 2=EDGE2(default) and 3=EDGE3; for Gen1
%       'G2' = 2nd Genration Edge
%   figPlot: if 1 or greater integer,plots 3D toolpath into figure(figPlot)
%       Othewise, no plot is created. 
function [TipL TipR] = EDGE_Kinematics(D, kVer, figPlot)
% kVer=2;

fprintf(1, '\nComputing Kinematics ... ');


if (size(D,2)~=13)
    error('Incorrect data format, cannot compute kinematics');
end

%% Set up constants 
Time=1; Q1L=2; Q2L=3; dL=4; Q3L=5 ; QgL=6 ; FgL=7 ; 
        Q1R=8; Q2R=9; dR=10;Q3R=11; QgR=12; FgR=13; 
        
c=pi/180;
% These are the link parameters...these MUST be accurate.  
alpha12L = deg2rad(75);
alpha23L = deg2rad(59.89); % should be 60, Solidworks mishap made it 59.89
alpha12R = -deg2rad(75);
alpha23R = -deg2rad(59.89); % should be 60, " " ""
% These are the anticipated initial pose (mostly wrong for non LabView
% versions of EDGE);  The errors here are fixed by the icpFix below.
t1L_hat = -0.535; %-0.49;% rad; %Found by rough measrument with physical system (approx)
t2L_hat = -1.43;% rad; %Found by rough measrument with physical system (VERY approx)
t3L_hat = -.050;% rad; %Total guess based on how people's hand actually fits
d4L_hat = -23.57;% cm;  %Found by rough measrument with physical system (approx)
t1R_hat = 0.535;% rad; %Found by rough measrument with physical system (approx)
t2R_hat = 1.43;% rad; %Found by rough measrument with physical system (VERY approx)
t3R_hat = .050;% rad; %Total guess based on how people's hand actually fits
d4R_hat = -23.57;% cm;  %Found by rough measrument with physical system (approx)
axisOffset_ClipApplier = 0; %allocate
clipApplierKinematics = false;

% Select the icpFix and other EDGE version&device-specific params
switch kVer
    % version 2 (G2, 2nd Generation EDGE's had different homebase pos.
    case {'G2','g2', 'G2clipApplier', 'g2clipApplier'} % Note, you'll need to recheck these at the bottom of the script
        % hbPosG2=[NaN -29.68 26.37 -7.5 0 NaN NaN -29.68 26.37 -7.5 0 NaN NaN];
%         t1L_hat = t1L_hat + -29.68*c;
%         t2L_hat = t2L_hat +  26.37*c;
%         t3L_hat = t3L_hat +   0.00*c;
%         d4L_hat = d4L_hat +  7.50 ;
% 
%         t1R_hat = t1R_hat +  -29.68*c;
%         t2R_hat = t2R_hat +  -26.37*c;
%         t3R_hat = t3R_hat +   0.00*c;
%         d4R_hat = d4R_hat +  7.50;
    
               
        % new from LLoyd
        k = pi/180.0; % for coversion between deg to radians
        alpha12L =  k*75.0;
        alpha23L =  k*60.0;  %k*59.89; % should be 60, Solidworks mishap made it 59.89; Lloyd says Gen2 EDGE fixed it
        alpha12R = -k*75.0 *(-1);
        alpha23R = -k*60.0 *(-1);  %-k*59.89% should be 60, " " ""
        % New values: Determined directly, carefully from sw models (joint angles at home base position, Generation 2)
        % All length units are in cm, all angles degrees converted to radians(i.e., *k) for faster processing later on.
        ballTipOffset = -0.776; % As measured in soldiworks model: the distance from ball tip center to plane that tool tip butts against; All tool block ref traces should use the ball tip
        toolTipDelta_Grasper = 0.0; % There reference for differences to tool length: Maryland grasper tool insert in its tool; Found by careful physical measurement
        toolTipDelta_CurvedShear = -0.14351; % So small, we'll just neglect it.
        toolTipDelta_NeedlDriver = 2.01295; % cm
        toolTipDelta_ClipApplier = 1.58496 ;     % currently unknown
        axisOffset_ClipApplier = 0.30051502;
        
        
        % If 'clip' is found, use clip-applier kinematics values and offset
        if ~isempty (findstr('clip', lower(kVer)))
           clipApplierKinematics = true; 
           ballTipOffset = ballTipOffset + toolTipDelta_ClipApplier;
           warning('Using Clip Applier Kinematics ...');
        else
            clipApplierKinematics = false;
            
            % if 'needle' is found, apply the needle driver length offset
            if ~isempty (findstr('needle', lower(kVer)))
                ballTipOffset = ballTipOffset + toolTipDelta_NeedlDriver;                
            end            
        end
            
        t1L_hat = -68.79*k;
        t2L_hat = -54.12*k;
        d4L_hat = -17.7828194  + ballTipOffset; % ref block always uses ball tips
        t1R_hat = -68.79*k*(1);
        t2R_hat = -54.12*k*(1);
        d4R_hat = -17.7828194 + ballTipOffset; % ref block always uses ball tips
       

        
%         TcomboL =[
%         -0.1313    0.5449   -0.8281   -6.7284
%         -0.8735    0.3314    0.3565  -13.6276
%          0.4687    0.7702    0.4325   13.0491
%          0         0         0        1.0000];
%         
%         icpFix.RfitL = [
%             TcomboL(1:3, 1:3);%eye(3)
%             ];
%         icpFix.TfitL = [
%             TcomboL(1:3, 4);
%             ];
        
        % Use identity matrix only for G2; the correction matrix is found
        % later elsewhere.
        icpFix.RfitR = [
            eye(3)
            ];
        icpFix.TfitR = [
            0
            0
            0
            ];
        
        icpFix.RfitL = icpFix.RfitR;
        icpFix.TfitL = icpFix.TfitR;
        fprintf(1,'Using Gen2 ver...');
        
    % EDGE1 for UMN/NOLA studies
    case 1
        icpFix.RfitL = [
            0.999644 -0.00570911 -0.0260549
            0.00555772 0.999967 -0.00587899
            0.0260876 0.0057321 0.999643
            ];
        icpFix.TfitL = [
            1.30839
            1.5174
            -0.0845075
            ];
        icpFix.RfitR = [
            0.998876 0.0381344 -0.028162
            -0.0370169 0.998546 0.0391892
            0.0296155 -0.0381027 0.998835
            ];
        icpFix.TfitR = [
            0.0926862
            1.29623
            -0.254456
            ];
        fprintf(1,'Using EDGE1 UMN/NOLA ver...');

    % EDGE3 for UMN/NOLA studies
    case 3
        icpFix.RfitL = [
            0.999447 0.0255452 0.021308
            -0.025908 0.999521 0.0169271
            -0.0208654 -0.0174698 0.99963
            ];
        icpFix.TfitL = [
            0.866233
            0.493905
            0.214628
            ];
        icpFix.RfitR = [
            0.99948 0.0149005 -0.0285802
            -0.0148515 0.999888 0.00192586
            0.0286057 -0.0015004 0.99959
            ];
        icpFix.TfitR = [
            0.395366
            0.760924
            0.134333
            ];
        % homebases were a little offset in EDGE 3
        d4L_hat = d4L_hat +.5;
        d4R_hat = d4R_hat +.5;
        fprintf(1,'Using EDGE3 UMN/NOLA ver...');
    
    % EDGE LabView version (used in UWS);  This should induce NO added
    % Rot or Translation since the baseline kinematics should be used.  
    % NOTE: in the LabView version (UWS), there was significant play in the
    % homebase switch, which induced about 1cm of error in dL/dR initial
    % values.  This causes inconsistencies between start/stop posititions
    % of tool tip in tasks or between tasks.  DO NOT use ICP to "fix" this.
    % The baseline kinematics are the best possible version. 
    case 0
        icpFix.RfitL = [eye(3)];
        icpFix.TfitL =[
            0
            0
            0];
        icpFix.RfitR =[ eye(3)];
        icpFix.TfitR =[
            0 
            0
            0];
        fprintf(1,'Using LabView ver...');
   
        
    % EDGE2 for UMN/NOLA studies and default
    %case 2      
    otherwise
        icpFix.RfitL = [
            0.9990    0.0429    0.0090
            -0.0430    0.9989    0.0158
            -0.0083   -0.0162    0.9998];
        icpFix.TfitL =[
            0.9966
            1.0044
            -0.3750];
        icpFix.RfitR =[
            0.9999    0.0078   -0.0121
            -0.0071    0.9984    0.0557
            0.0125   -0.0556    0.9984];
        icpFix.TfitR =[
            0.0985
            1.2173
            -0.1108];
        if(kVer == 2)
            fprintf(1,'Using Edge2 UMN/NOLA ver ...');
        else
            fprintf(1,'Using Default ver (Edge2)...');
        end
end

% ICP fix (found via EdgeRefBlocks.m and related files):
% QUOTE: 
%     Conclusion: Best data is from:
%     EDGE2_FinalAnalysis_DiagRefBlock_ConstRotTipKnuckledDown_TimK_11.28.2011-
%     22.58.46-trimmed.txt
%     The resulting fixes are:
%     icpFix.RfitL =
%         0.9990    0.0429    0.0090
%        -0.0430    0.9989    0.0158
%        -0.0083   -0.0162    0.9998
%     icpFix.TfitL =
%         0.9966
%         1.0044
%        -0.3750
%     icpFix.RfitR =
%         0.9999    0.0078   -0.0121
%        -0.0071    0.9984    0.0557
%         0.0125   -0.0556    0.9984
%     icpFix.TfitR =
%         0.0985
%         1.2173
%        -0.1108
%     In terms of roll/pitch/yaw:  (use '[y p r] =dcm2angle( RfitR{f} )' )
%        LeftFix:   0.907638 -0.515495 2.45611 
%        RightFix:  3.19208   0.692078 0.44414 
%     These are the fixes to be use in EDGE_Kinematics.m: icpFix.

% For EDGE 3 
% //////////////////////////////////////////////////////////
% % For File: EDGE3_FinalAnalysis_DiagRefBlock_TimK_BeforeDecommissioningEdge3FIXEDHomebases_10.13.2011-14.06.22-trimmed.txt
% icpFix.RfitL = [
%     0.999993 0.000413843 -0.00369288
%     -0.000134262 0.997154 0.0753894
%     0.00371357 -0.0753884 0.997147
%     ];
% icpFix.TfitL = [
%     0.8995
%     -0.307596
%     0.610852
%     ];
% icpFix.RfitR = [
%     0.999464 0.0141929 -0.0295028
%     -0.0141717 0.999899 0.000927003
%     0.0295129 -0.000508401 0.999564
%     ];
% icpFix.TfitR = [
%     0.508455
%     0.422948
%     0.490037
%     ];
% % In terms of roll/pitch/yaw:  (use [y p r] =dcm2angle( RfitR{f} ) )
% %     LeftFix:   4.32363 0.211587 0.0237116 % [r p  y]
% %     RightFix:  0.0531365 1.69063 0.813574
% % These are the fixes to be use in EDGE_Kinematics.m: icpFix.
% % //////////////////////////////////////////////////////////
% 
% 
% % //////////////////////////////////////////////////////////
% % For File: EDGE1_UMN_Subj11_1_RefBlock__12.10.2010-04.53.22-trimmed.txt
% icpFix.RfitL = [
%     0.999644 -0.00570911 -0.0260549
%     0.00555772 0.999967 -0.00587899
%     0.0260876 0.0057321 0.999643
%     ];
% icpFix.TfitL = [
%     1.30839
%     1.5174
%     -0.0845075
%     ];
% icpFix.RfitR = [
%     0.998876 0.0381344 -0.028162
%     -0.0370169 0.998546 0.0391892
%     0.0296155 -0.0381027 0.998835
%     ];
% icpFix.TfitR = [
%     0.0926862
%     1.29623
%     -0.254456
%     ];
% In terms of roll/pitch/yaw:  (use [y p r] =dcm2angle( RfitR{f} ) )
%     LeftFix:   -0.336958 1.493 -0.327221 % [r p  y] 
%     RightFix:  2.24685 1.61378 2.18634 
% These are the fixes to be use in EDGE_Kinematics.m: icpFix.
% //////////////////////////////////////////////////////////





%[hbL hbR] = EDGE_Kinematics(zeros(1,13))

% Original guess for homebase position) ???
% hbL =[   -3.8740   22.0380   -5.8924];
% hbR =[    1.3386   21.9417   -6.2529];


% These were generated by a Solidworks API vba macro.  They are the Homogeneouse 
% Transformation Matrix between EDGE origin and Sphere center for either
% hand.
T_EDGE_ORG_to_SphereCenterL = 1.0e+002 * [
    -0.00000000000000   0.00565982635388  -0.00824417161660  -0.67184679157301
    -0.00829523437902   0.00460413820792   0.00316085399212  -0.85607950145680
    0.00558471902580   0.00683873358205   0.00469495861500   1.76693066043015
    0                  0                  0   0.01000000000000];
T_EDGE_ORG_to_SphereCenterR = 1.0e+002 * [
    0.00064166752366   0.00564102729416   0.00823207591411   0.67127560851876
    -0.00826610767782  -0.00432117300377   0.00360540257533  -0.85647124629220
    0.00559103985389  -0.00703606956603   0.00438565826464   1.76556649172967
    0                  0                  0                 0.01000000000000];

% Preallocate:
TipL = zeros(size(D,1) ,3);
TipR = zeros(size(D,1) ,3);

% Preset commonly computed trig vals...
c12L = cos(alpha12L);
c23L = cos(alpha23L);
s12L = sin(alpha12L);
s23L = sin(alpha23L);
c12R = cos(alpha12R);
c23R = cos(alpha23R);
s12R = sin(alpha12R);
s23R = sin(alpha23R);
a =  axisOffset_ClipApplier;
E = a/c23L; % error term due to clip applier offset.  NOTE, c23R and c23L must be equal except sign.

% Compute fwd Kinematics for each point
for i=1:size(D,1);
    %introduce OFFSETS!!!!
    % % The following are from the original LabView format...
    % %         t1L=  ALL_DATA(i,5) + t1L_hat;
    % %         t2L= -ALL_DATA(i,6) + t2L_hat;
    % %         t3L=  ALL_DATA(i,7) + t3L_hat;
    % %         d4L= -ALL_DATA(i,8) + d4L_hat;
    % %         t1R=  ALL_DATA(i,14) + t1R_hat;
    % %         t2R= -ALL_DATA(i,15) + t2R_hat;
    % %         t3R=  ALL_DATA(i,16) + t3R_hat;
    % %         d4R=  ALL_DATA(i,17) + d4R_hat;
    %         t1L=  ALL_DATA(i,5) + t1L_hat;
    %         t2L= -ALL_DATA(i,6) + t2L_hat;
    %         t3L=  ALL_DATA(i,7) + t3L_hat;
    %         d4L= -ALL_DATA(i,8) + d4L_hat;
    %         t1R=  ALL_DATA(i,14) + t1R_hat;
    %         t2R= -ALL_DATA(i,15) + t2R_hat;
    %         t3R=  ALL_DATA(i,16) + t3R_hat;
    %         d4R=  ALL_DATA(i,17) + d4R_hat;
    t1L=  D(i,Q1L)*c + t1L_hat; %Q1L
    t2L=  D(i,Q2L)*c + t2L_hat; %Q2L
    t3L=  D(i,Q3L)*c + t3L_hat; %Q3L (Rot)
    d4L= -D(i,dL )*1 + d4L_hat; %dL
    t1R= (-1)*-D(i,Q1R)*c + t1R_hat;%Q1R    ...same sign as Left
    t2R= (-1)*-D(i,Q2R)*c + t2R_hat;%Q2R    ...same sign
    t3R=  D(i,Q3R)*c + t3R_hat;%Q3R (Rot)
    d4R= -D(i,dR )*1 + d4R_hat;%dR          ...INVERT sign for right hand
    d4R = -d4R;

    % for trig shorthand
    c1L  = cos(t1L);
    c2L  = cos(t2L);
    s1L  = sin(t1L);
    s2L  = sin(t2L);
    c1R  = cos(t1R);
    c2R  = cos(t2R);
    s1R  = sin(t1R);
    s2R  = sin(t2R);
    
    
    % do the computation based on clip applier or standard forward
    % kinematics
    if ~clipApplierKinematics

        % original working implementation in EDGE/Jay software (no shorthand) 
        % P_L = [ ...
        %     [ (-cos(t1L)*sin(t2L)-sin(t1L)*cos(t2L)*cos(alpha12L))*sin(alpha23L)*d4L-sin(t1L)*sin(alpha12L)*cos(alpha23L)*d4L]
        %     [ (-sin(t1L)*sin(t2L)+cos(t1L)*cos(t2L)*cos(alpha12L))*sin(alpha23L)*d4L+cos(t1L)*sin(alpha12L)*cos(alpha23L)*d4L]
        %     [                                       -cos(t2L)*sin(alpha12L)*sin(alpha23L)*d4L+cos(alpha12L)*cos(alpha23L)*d4L]
        %     [                                                                                                               1]];
        % P_R = [ ...
        %     [ (-cos(t1R)*sin(t2R)-sin(t1R)*cos(t2R)*cos(alpha12R))*sin(alpha23R)*d4R-sin(t1R)*sin(alpha12R)*cos(alpha23R)*d4R]
        %     [ (-sin(t1R)*sin(t2R)+cos(t1R)*cos(t2R)*cos(alpha12R))*sin(alpha23R)*d4R+cos(t1R)*sin(alpha12R)*cos(alpha23R)*d4R]
        %     [                                       -cos(t2R)*sin(alpha12R)*sin(alpha23R)*d4R+cos(alpha12R)*cos(alpha23R)*d4R]
        %     [                                                                                                               1]];
        
        % trig shorthand ... (Note, left should be identical to right
        % except for L / R joint angle substitutions ... t1L --> t1R , etc.
        
        P_L = [ ...
            (- c1L*s2L - c12L*c2L*s1L)*d4L*s23L + (-c23L*s12L*s1L)*d4L 
            (c12L*c1L*c2L - s1L*s2L)  *d4L*s23L + ( c23L*c1L*s12L)*d4L
            (-c2L*s12L)*d4L*s23L                + ( c12L*c23L    )*d4L
            1                                                           ];
        P_R =[ ...
            (- c1R*s2R - c12R*c2R*s1R)*d4R*s23R + (-c23R*s12R*s1R)*d4R 
            (c12R*c1R*c2R - s1R*s2R)  *d4R*s23R + ( c23R*c1R*s12R)*d4R
            (-c2R*s12R)*d4R*s23R                + ( c12R*c23R    )*d4R
            1                                                           ];
    else
        
        % Using Clip Applier Kinematics (only applies to right hand) ...
        P_L = [ ...
            (- c1L*s2L - c12L*c2L*s1L)*d4L*s23L + (-c23L*s12L*s1L)*d4L 
            (c12L*c1L*c2L - s1L*s2L)  *d4L*s23L + ( c23L*c1L*s12L)*d4L
            (-c2L*s12L)*d4L*s23L                + ( c12L*c23L    )*d4L
            1                                                           ];
        P_R =[ ...
            (- c1R*s2R - c12R*c2R*s1R)*d4R*s23R + (-c23R*s12R*s1R)*d4R - (s12R*s1R)*E %+ (c1R*c2R - c12R*s1R*s2R)*a
            (c12R*c1R*c2R - s1R*s2R  )*d4R*s23R + ( c23R*c1R*s12R)*d4R + (c1R*s12R)*E %+ (c2R*s1R + c12R*c1R*s2R)*a
            (-c2R*s12R)*d4R*s23R                + ( c12R*c23R    )*d4R +  c12R*E      %- (s12R*s2R)*a
            1                                                           ];
        fprintf('%g    %g   %g \n', [ (c1R*c2R - c12R*s1R*s2R)*a    (c2R*s1R + c12R*c1R*s2R)*a   (s12R*s2R)*a ]);
    end

    XL = T_EDGE_ORG_to_SphereCenterL * P_L;
    XR = T_EDGE_ORG_to_SphereCenterR * P_R;
    %XR = T_EDGE_ORG_to_SphereCenterR*[-1 0 0 0; 0 -1 0 0; 0 0 1 0 ; 0 0 0 1] * P_R; % try a mirror
    %myT=[  0.6332    0.1265   -0.7636    9.4498
    %0.7738   -0.0807    0.6283  -19.4620
    %0.0178   -0.9887   -0.1490    1.1157
    %     0         0         0    1.0000];
    %XR = myT*P_R;

    %         XL =  P_L;
    %         XR =  P_R;

    %T_TipL(i, :)=XL(1:3)' + [  61. 75 -159];
    %T_TipR(i, :)=XR(1:3)' + [ -60 75.7 -159.3];
%     TipL(i,:) = XL(1:3)' + [  60 75 -159];
%     TipR(i,:) = XR(1:3)' + [ -60 75 -159];

%     TipL(i,:) = XL(1:3)' + [  67 85 -177]; 
%     TipR(i,:) = XR(1:3)' + [ -67 85 -177];   
         
%     TipL(i,:) = XL(1:3)' + [  60 75 -159] + -hbL;
%     TipR(i,:) = XR(1:3)' + [ -60 75 -159] + -hbR;

%     TipL(i,:) = XL(1:3)' + [  66 79 -172]; 
%     TipR(i,:) = XR(1:3)' + [ -66 79 -172];

    % Subtract supposed position of sphere centers relative to ref-frame.
    TipL(i,:) = XL(1:3)' ;%+ [  60 75 -159];
    TipR(i,:) = XR(1:3)' ;%+ [ -60 75 -159];
    
    % Apply IcpFix to account for all "guess" errors above...  
     TipL(i,:) = [(icpFix.RfitL)*TipL(i,:)' + icpFix.TfitL]';
     TipR(i,:) = [(icpFix.RfitR)*TipR(i,:)' + icpFix.TfitR]';
     
     %TTT = [ icpFix.RfitL icpFix.TfitL; [0 0 0 1] ] * [eye(3) [  60 75 -159]'; [0 0 0 1 ] ]  *  T_EDGE_ORG_to_SphereCenterL 
     %[TipL(i,:)' ; 1] - TTT*P_L
     XrawL(i,:) = P_L;
     XrawR(i,:) = P_R;

     XsphereL(i,:) = XL ;%+ [  60 75 -159];
     XsphereR(i,:) = XR ;%+ [ -60 75 -159];
end % of for i



% % From ESPcontroller.cpp (->loadDefaultSettings()) ...
% % // Should add the delta here between 1st and 2nd generation homebase...
% % // the two generations had homebase in different locations.
% % sensorCal[i].homePoseOrigin[LEFT][J1 ] = -29.68;//0.0;
% % sensorCal[i].homePoseOrigin[LEFT][J2 ] =  26.37;//0.0;
% % sensorCal[i].homePoseOrigin[LEFT][LIN] = -7.5;//0.0;
% % sensorCal[i].homePoseOrigin[RGHT][J1 ] = -29.68;//0.0;
% % sensorCal[i].homePoseOrigin[RGHT][J2 ] =  26.37;//0.0;
% % sensorCal[i].homePoseOrigin[RGHT][LIN] = -7.5;//0.0;
% % hbPosG2=[NaN -29.68 26.37 -7.5 0 NaN NaN -29.68 26.37 -7.5 0 NaN NaN];
% % [hbDeltaG2L hbDeltaG2R] = EDGE_Kinematics(hbPosG2);
% 
% hbDeltaG2L =[   -2.8376    9.2681    2.8317];
% hbDeltaG2R =[    1.8716    9.2880    2.6167];
% 
% for i=1:3
%     TipL(:,i) = TipL(:,i) + hbDeltaG2L(i);
%     TipR(:,i) = TipR(:,i) + hbDeltaG2R(i);
% end
disp(' DONE.');

% Plot 3-D toolpath if requested
if( figPlot>0)
    colordef white;   
    figure(figPlot);
    
    %plot both hands
    ms=8;
    
    figure(figPlot); clf
    plot3(XrawR(:,1),XrawR(:,2),XrawR(:,3), '.r', 'markersize', ms); hold on;
    plot3(XsphereR(:,1),XsphereR(:,2),XsphereR(:,3), '.m', 'markersize', ms); hold on;
    plot3(TipR(:,1),TipR(:,2),TipR(:,3), 'r', 'markersize', ms); hold on;
    legend( 'Raw-R', 'SphereCtrAdj-R', 'Final-R')
    xlabel('x');     ylabel('y');    zlabel('z')
    grid on;    axis vis3d
    title('EDGE\_Kinematics Internal Variables (Right)')
    
    figure(figPlot+1);clf
    plot3(XrawL(:,1),XrawL(:,2),XrawL(:,3), '.b', 'markersize', ms); hold on;
    plot3(XsphereL(:,1),XsphereL(:,2),XsphereL(:,3), '.c', 'markersize', ms); hold on;
    plot3(TipL(:,1),TipL(:,2),TipL(:,3), 'b', 'markersize', ms); hold on;
    legend( 'Raw-L', 'SphereCtrAdj-L', 'Final-L')
    xlabel('x');     ylabel('y');    zlabel('z')
    grid on;    axis vis3d
    title('EDGE\_Kinematics Internal Variables (Left)')
    
    
    
    
    
end

% Select what format to output
switch kVer
    % version 2 (G2, 2nd Generation EDGE's )
    case {'G2','g2', 'G2clipApplier', 'g2clipApplier'} 
        %only spit out raw kinematics; no extra base frame Tmtx
        TipL = XrawL(:,1:3);
        TipR = XrawR(:,1:3);
    otherwise
        % return corrected tip position 
        return
end
function RBF_Z = getRBFdata(SegData,group,leaveout,SepHands)
%getRBFdata: for each grasp indepently compute the RBF from derivatives
%SegData{}: Edge data structure for grasps
%Group(): which elements of SegData to use 1:3 (nov,int,exp)
%leaveout: Which surgeon to leave out (0 will result in all being used)
%rbfstates: 'accjerk' or 'velacc'

if(SepHands)
    RBF_Z{1} = [];
    RBF_Z{2} = [];
else
   RBF_Z = []; 
end

    for ii = 1:length(SegData{group}.Trial)  %loop through all expert surgeons (keep some)
        %only use surgeon who is not being left out
        if(ii ~= leaveout)
            for hh = 1:length(SegData{group}.Trial{ii}.Hand)
                
                for ss = 1:length(SegData{group}.Trial{ii}.Hand{hh}.Segment) %all segments from that surgeon
                    %grab the state information

                    if(SepHands)
                        RBF_Z{hh} = [RBF_Z{hh}; SegData{group}.Trial{ii}.Hand{hh}.Segment{ss}.RBF];
                    else
                        RBF_Z = [RBF_Z; SegData{group}.Trial{ii}.Hand{hh}.Segment{ss}.RBF];
                    end
                end
            end
        end
    end
    
end
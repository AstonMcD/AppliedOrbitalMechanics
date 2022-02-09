%
GraceA;
GraceB;
    j = 0;
    k=0;
    l=0;
    m=0;
    C = cell(length(GraceA),4);
for i = 1:length(GraceA)
    DateString = datestr(datenum([GraceA(i,1),GraceA(i,2)+1,GraceA(i,3)]));
%Both Activated
    if GraceA(i,4) == 1 && GraceB(i,4) == 1  
        j = j + 1;    
        BothActivated(j) = DateString;
        C{j,1} = DateString;
%Only Grace A Activated
    elseif GraceA(i,4) == 1 && GraceB(i,4) == 0
        k = k + 1;
        OnlyGraceA_Activated(k) = DateString;
        C{k,2} = DateString;
%Only Grace B Activated
    elseif GraceA(i,4) == 0 && GraceB(i,4) == 1
        l = l + 1; 
         OnlyGraceB_Activated(l) = DateString;
        C{l,3} = DateString;
%Neither Activated
    else
        m = m + 1;
        C{m,4} = DateString;
        NeitherActivated{m} = DateString;
    end
end
writecell(C,'ActivationResults.xls')
% DataCollection = zeros(3,length(GraceA))
% DataCollection(1,1:length(BothActivated)) = BothActivated;
% DataCollection(2,1:length(NeitherActivated)) = NeitherActivated;
% DataCollection(3,1:length(OnlyGraceA_Activated)) = OnlyGraceA_Activated;

%       T = table(BothActivated, NeitherActivated, OnlyGraceA_Activated); 
%       filename = 'GraceAB_ResultsOfActivation.xlsx';
%     writetable(T,filename,'Sheet',1,'Range','B1')
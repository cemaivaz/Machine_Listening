%Functionthatlocatesthefundamentalfrequencyoftiaccordingtoinput
%thresholdsTThandATh.

%??Input??
%TTh:Harmonicthreshold
%ATh:Amplitudethreshold
%ti:An?by?2matrix[index_of_minima,value_of_minima]
%Row1MUSTbethecurrentpickedfundamentalfrequencyminima

%??Output??
%tg:Newminimaofthefundamentalfrequency
%ti:Remainingminimas

function[tg,ti]=findf0(TTh,ATh,ti)
N=[2 3 4 5 6];%Thresholdconstantstotestforharmonics

tg=ti(1,:);%Extractcurrentfundamentalfrequencyofminimas
ti(1,:)=[];%Deletesthecurrentfundamentalfrequencyofminimasfromt
ti=sortrows(ti,1);%Sortsrowsaccordingtoindexes(ascending)
ti=ti(end:-1:1,:);%Reversesti...sortmustbedesending
index=[];%Musttoallocatedtoavoidruntimewarnings

%Testforcandidatesofharmonicminimastothefundamentalfrequencyminima
for i=1:length(N)%Runelementsinn?times
x=tg(1)./ti(:,1);%Calculatestg/ti
newIndexes=find(abs(N(i)-x)<TTh);%Findsindexesofti,thatsatifiesthethresholdTTh
if length(newIndexes)>0
index=[index;newIndexes];%Addsthefoundindextoageneralindexvariable
end
end

%Testamongthecandidatesofaminimastothefundamentalfrequencyminima
%TestthecandidatesifthethredsholdiswithinthethresholdATh
newTgIndex=0;%Indexamongindexofnewfundamentalfrequencyminima
for i=1:length(index)%Runsonlyamongtheharmoniccandidates
if ti(index(i),2)-tg(2)<ATh%TestifAThthresholdifsatified
newTgIndex=i;%Saveindexof’index’ofnewfundamentalfrequencyminima
end
end
if newTgIndex>0%Ifanewfundamentalfrequencywasfound,setnewtg
tgTemp=ti(index(newTgIndex),:);%Makescopyofnewfundamentalfrequencyminima
ti(index(i),:)=tg;%Copiesoldfundamentalfrequencyminimaintoti
tg=tgTemp;%Copiesnewfundamentalfrequencyminimaintotg
ti=sortrows(ti,1);%Sortsrowsaccordingtoindexes(ascending)
ti=ti(end:-1:1,:);%Reversesti...sortmustbedesending
end

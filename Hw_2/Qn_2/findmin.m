%FunctionlocatesnMinsofthesmallestminimasofdata,whereminimas
%havetobeseparatedwithatleastspredsamples.Alsoaminimacan’t
%excidemaxMinValue!

%??Input??
%data:Datatolocateminimasin
%nMins:#ofminimastofindindata(ifpossible)
%spred:Minimum#ofsamplesbetweenminimaslocated
%maxMinValue:Maximumvalueofaminima

%??Output??
%ti:An?by?2matrixwherecolumn1isindexesofminimaand
%column2isvaluesofminima.
%minCount:Totalnumberoffoundminimaindata

function[ti minCount]=findmin(data,nMins,spred,maxMinValue)

[nm]=size(data);%Flipsvectorcorrect!
if m==1
data=data';
end

minimas=localminima(data)';%LocatesALLminimasindataasindexes
minimas(:,2)=data(minimas(:,1))';%Extendsminstoholddatavaluestoindexesofminimas
minCount=length(minimas);%Totalnumberoflocalminsfound

indexes=find(minimas(:,2)>=maxMinValue);%FindindexesofminimathatexcidesmaxMinValue
minimas(indexes,:)=[];%DeletesminimasthatexcidemaxMinValue
[m n]=size(minimas);%Newsizeofminsmatrix

if m>=1%Onlyrunsifdataareavailable
ti=zeros(nMins,2);%Allocatesnumberofsmallestminimaswanted
for i=1:nMins
[minima mIndex]=min(minimas(:,2));%Findsglobalminamongminimas
ti(i,:)=minimas(mIndex,:);%Copiesfoundglobalminimaintronewmatrice
indexes=find(abs(minimas(mIndex,1)-minimas(:,1))<spred);%Findsindexesofminimaswithin’spred’samples offoundminima
minimas([mIndexindexes'],:)=[];%Deletesminimawithin’spred’samplesoffoundminima
if length(minimas)<=0%Breakloopifnomoredata
ti(i+1:end,:)=[];%Deletesusedallocatedspace
break;
end
end
else
ti=zeros(0,2);%Allocatesa0?by?2matrix(toavoidruntimeerror)
end

%Thisfunctioncalculatesthetoppointofacurvefittedparable.
function[xx yy]=parreg(x,y)
%x:x?koordinatestodocurvefittingaround
%y:y?koordinatestodocurvefittingaround

[m n]=size(y);%Vectorxhastostand(am?by?1vector)
if n~=1
y=y';
end
[m n]=size(x);%Vectoryhastostand(am?by?1vector)
if n~=1
x=x';
end

%CurvefittingintheformAx=b
X=[ones(size(x)) x x.^2];%Xvalues
a=X\y;%Leastsqauressolutionofinclination
xx=-a(2)/(2*a(3));%Calculatesxoftoppointofparable
yy=[1 xx xx.^2]*a;%Calculatestoftoppointofparable

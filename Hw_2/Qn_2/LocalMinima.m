%mins=LocalMinima(x)
%
%findspositionsofallstrictlocalminimaininputarray

function mins=LocalMinima(x)

nPoints=length(x);
Middle=x(2:(nPoints-1));
Left=x(1:(nPoints-2));
Right=x(3:nPoints);
mins=1+find(Middle<Left&Middle<Right);

%WrittenbyKennethD.Harris
%ThissoftwareisreleasedundertheGNUGPL
%www.gnu.org/copyleft/gpl.html
%anycomments,orifyoumakeanyextensions
%letmeknowatharris@axon.rutgers.edu

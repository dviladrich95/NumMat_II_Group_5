% Local P1 mass matrix
function m=localmass2D(Fdet)

m=Fdet*[1 1/2 1/2;1/2 1 1/2;1/2 1/2 1]/12;

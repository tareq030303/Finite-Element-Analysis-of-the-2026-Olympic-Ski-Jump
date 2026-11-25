function [file_i,xy,nnod,sizee,idb,ndof,incidenze,l,gamma,m,EA,EJ,posiz,nbeam,pr]=loadstructure

file_i = 'input.inp';
fid = fopen(file_i, 'r');

% Lettura Nodi
findcard(fid, '*NODES');
nnod = 0;
while true
    riga = scom(fid);
    if contains(riga, '*ENDNODES')
        break;
    end
    dati = sscanf(riga, '%f');
    nnod = nnod + 1;
    vincoli(nnod, :) = dati(2:4)';
    xy(nnod, :) = dati(5:6)';
end

% Costruzione matrice idb (Gradi di libertà)
idb = zeros(nnod, 3);
ndof = 0;
for i = 1:nnod
    for j = 1:3
        if vincoli(i, j) == 0
            ndof = ndof + 1;
            idb(i, j) = ndof;
        end
    end
end
ndof_tot = ndof;
for i = 1:nnod
    for j = 1:3
        if vincoli(i, j) == 1
            ndof_tot = ndof_tot + 1;
            idb(i, j) = ndof_tot;
        end
    end
end

% Lettura Travi (*BEAMS)
findcard(fid, '*BEAMS');
nbeam = 0;
while true
    riga = scom(fid);
    if contains(riga, '*ENDBEAMS')
        break;
    end
    dati = sscanf(riga, '%f');
    nbeam = nbeam + 1;
    n1 = dati(2);
    n2 = dati(3);
    pr(nbeam, 1) = dati(4);

    incidenze(nbeam, :) = [idb(n1, :) idb(n2, :)];
    posiz(nbeam, :) = xy(n1, :);
    dx = xy(n2, 1) - xy(n1, 1);
    dy = xy(n2, 2) - xy(n1, 2);
    l(nbeam, 1) = sqrt(dx^2 + dy^2);
    gamma(nbeam, 1) = atan2(dy, dx);
end

% Lettura Proprietà (*PROPERTIES)
findcard(fid, '*PROPERTIES');
while true
    riga = scom(fid);
    if contains(riga, '*ENDPROPERTIES')
        break;
    end
    dati = sscanf(riga, '%f');
    id_prop = dati(1);
    prop_tab(id_prop, :) = dati(2:4)';
end

for k = 1:nbeam
    m(k, 1)  = prop_tab(pr(k), 1);
    EA(k, 1) = prop_tab(pr(k), 2);
    EJ(k, 1) = prop_tab(pr(k), 3);
end

sizee = max(max(xy) - min(xy));
fclose(fid);
end
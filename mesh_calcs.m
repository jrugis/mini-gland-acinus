%
% EXAMPLE mesh file calculations
%
% J.Rugis
% 05.08.21
%
%

% read the  mesh ply file
fname = 'AcinusSevenCells.ply';
%[acinus, duct, verts, faces, tets, lnodes, lradii, lsegs] = read_ply(fname);

% get the number of acinii
nacinii = size(acinus,2);

% get the cell data for the nth acinus
n = 1;
first_cell = acinus(n).icells + 1;
last_cell = acinus(n).icells + acinus(n).ncells;
averts = verts(first_cell : last_cell);
afaces = faces(first_cell : last_cell);
atets = tets(first_cell : last_cell);




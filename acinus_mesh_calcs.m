%
% EXAMPLE mesh file calculations for an acinus
%
% J.Rugis
% 05.08.21
%
%

function [averts, afaces, atets, face_areas, face_centers, tet_volumes, tet_centers] = acinus_mesh_calcs(acinus_n, fname)

  % read the mesh ply file
  [acinus, duct, verts, faces, tets, lnodes, lradii, lsegs] = read_ply(fname);

  % get the cell data for the nth acinus
  n = acinus_n;
  first_cell = acinus(n).icells + 1;
  last_cell = acinus(n).icells + acinus(n).ncells;
  averts = verts(first_cell : last_cell);
  afaces = faces(first_cell : last_cell);
  atets = tets(first_cell : last_cell);

  % get the lumen data for the nth acinus
  first_segment = acinus(n).ilsegs + 1;
  last_segment = acinus(n).ilsegs + acinus(n).nlsegs;
  alsegs = lsegs(first_segment : last_segment,:);

  % calculate surface triangle (face) attributes
  face_centers = cell(1,acinus(n).ncells);
  face_areas = cell(1,acinus(n).ncells);
  for i = 1:acinus(n).ncells
    nfaces = size(afaces{i},1);
    face_centers{i} = zeros(nfaces,3);
    face_areas{i} = zeros(nfaces,1);
    for j = 1: nfaces
      [face_areas{i}(j), face_centers{i}(j,:)] = calc_tri(afaces{i}(j,:),averts{i});
    end
  end
  
  % calculate tetrahedron attributes
  tet_centers = cell(1,acinus(n).ncells);
  tet_volumes = cell(1,acinus(n).ncells);
  for i = 1:acinus(n).ncells
    ntets = size(atets{i},1);
    tet_centers{i} = zeros(ntets,3);
    tet_volumes{i} = zeros(ntets,1);
    for j = 1: ntets
      [tet_volumes{i}(j), tet_centers{i}(j,:)] = calc_tet(atets{i}(j,:),averts{i});
    end
  end

  % calculate overall cell attributes

  %cell_volumes = zeros(acinus(n).ncells,1);
  % for cell volume use abs(sum( p1.Dot(p2.Cross(p3)) )) over surface tris
  % check against sum of tet volumes
  
  %cell_centroids = acinus(n).ncells;
  % for cell centroid use mean of surface tri centers

end

% calculate triangle area and center
function [tarea, tctr] = calc_tri(pts, verts)
  v1 = verts(pts(1)+1,:);
  v2 = verts(pts(2)+1,:);
  v3 = verts(pts(3)+1,:);
  tarea = norm(cross(v1-v2, v1-v3)) / 2.0;
  tctr = mean([v1; v2; v3]);
end

% calculate tetrahedron volume and center
function [tvol, tctr] = calc_tet(pts, verts)
  v1 = verts(pts(1)+1,:);
  v2 = verts(pts(2)+1,:);
  v3 = verts(pts(3)+1,:);
  v4 = verts(pts(4)+1,:);
  tvol = abs(det([v2 - v1; v3 - v1; v4 - v1])) / 6.0;
  tctr = mean([v1; v2; v3; v4]);
end


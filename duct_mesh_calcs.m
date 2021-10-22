%
% EXAMPLE mesh file calculations for a duct
%
% J.Rugis
% 06.08.21
%
% cell_geom is a cell array of struct containing the following features
% cell_raw    - api_idx      [1, n_api_mesh]
%             - bas_idx      [1, n_bas_mesh]
%             - lat_idx      [1, n_lat_mesh]
%             - face_coord   [9, n_mesh]
%             - face_area    [1, n_mesh]
%             - cell_type    0 = intercalated, 1 = striated               
%             - volume
%
% lumen_geom  - segment      [1, n_int+1]
%             - L
%             - start
%             - end
%             - length
%             - n_int
%             - radius
%             - volume

function [cell_geom, lumen_geom] = duct_mesh_calcs(mesh_file, version)

% read the mesh ply file
[acinus, duct, verts, faces, tets, lnodes, lradii, lsegs] = read_ply(mesh_file, version);

% get the cell data for the duct
n = 1; % there is only one duct in the current model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% striated cells only for now
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% intercalated cells
%first_icell = duct(n).iicells + 1;
%last_icell = duct(n).iicells + duct(n).nicells;
%ifaces = faces(first_icell : last_icell);

% striated cells
first_scell = duct(n).iscells + 1;
last_scell = duct(n).iscells + duct(n).nscells;
sfaces = faces(first_scell : last_scell);

n_cell = duct(n).nscells;   % number of striated cells
cell_geom = cell(n_cell,1); % pre-alocate matlab cell array
for i = 1:n_cell            % iterate over the cells
  sf = size(sfaces{i},1);
  face_coord = zeros(sf,9);
  face_area = zeros(sf,1);
  face_center = zeros(sf,3);
  face_type = zeros(sf,1);  % 1=apical, 2=basolateral, 3=basal
  for j = 1:sf              % iterate over the cell faces
    c1 = verts(sfaces{i}(j,1)+1,:);
    c2 = verts(sfaces{i}(j,2)+1,:);
    c3 = verts(sfaces{i}(j,3)+1,:);
    face_coord(j,:) = [c1, c2, c3];
    [face_area(j), face_center(j,:)] = calc_tri(face_coord(j,:));
    face_type(j) = get_face_type(face_center(j,:), 0.8, 12.5);
  end
  cell_center = calc_cell_center(face_center); % use the face centers
  
  %cell_volumes = zeros(acinus(n).ncells,1);
  % for cell volume use abs(sum( p1.Dot(p2.Cross(p3)) )) over surface tris
  % check against sum of tet volumes
  
  cell_raw.center = cell_center; % structure member assignment
  %cell_raw.api_idx = api_idx;
  %cell_raw.bas_idx = bas_idx;
  %cell_raw.lat_idx = lat_idx;
  cell_raw.face_coord = face_coord;
  cell_raw.face_area = face_area;
  %cell_raw.volume = volume;
  cell_geom{i} = cell_raw;       % save data for this cell
end

lumen_geom = cell(1);

  % get the lumen data for the nth acinus
  %first_segment = acinus(n).ilsegs + 1;
  %last_segment = acinus(n).ilsegs + acinus(n).nlsegs;
  %alsegs = lsegs(first_segment : last_segment,:);

end

% determine face type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: temporary simple version for now
%    Assumes duct aligned with z-axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ftype] = get_face_type(face_center, inner, outter)
  x = face_center(1);
  y = face_center(2);
  d = norm(x,y);
  if d-inner < 0.5      % apical
    ftype = 1;
  elseif outter-d < 1.5 % basal 
    ftype = 3;
  else                  % basolateral
    ftype = 2;
  end
end

% calculate triangle area and center
function [tarea, tctr] = calc_tri(verts)
  v1 = verts(1:3);
  v2 = verts(4:6);
  v3 = verts(7:9);
  tarea = norm(cross(v1-v2, v1-v3)) / 2.0;
  tctr = mean([v1; v2; v3]);
end

% calculate center of cell bounding box
function [center] = calc_cell_center(verts)
  cmax = max(verts,[],1);
  cmin = min(verts,[],1);
  center = mean([cmax; cmin]);
end


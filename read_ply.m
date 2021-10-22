%
% read in a striated cell mesh file
%
% J.Rugis
% 27.09.21
%
%

function [verts, faces, ] = read_ply(fname)
  % check ply file version
  pfile = fopen(fname,'r');
  V = get_version(pfile);
  if ~strcmp(V, '1.2')
    disp('ERROR: incorrect ply file version.');
    return;
  end
  
  fgetl(pfile);
  fgetl(pfile);
  
  % get the data counts
  pfile = fopen(fname,'r');
  nacinus = get_count(pfile, 'acinii');           
  nduct = get_count(pfile, 'duct');
  nlnode = get_count(pfile, 'lumen_node');
  nlseg = get_count(pfile, 'lumen_segment');
  ncell = get_count(pfile, 'cell');
  nvert = get_count(pfile, 'vertex');
  skip_header(pfile);

  % get acinii info
  acinus = struct([]);
  for i = 1:nacinus
    tokens = str2double(split(fgetl(pfile)));
    acinus(i).ncells = tokens(1);
    acinus(i).icells = tokens(2);
    acinus(i).nlsegs = tokens(3);
    acinus(i).ilsegs = tokens(4);
  end
  
  % get duct info
  duct = struct([]);
  for i = 1:nduct
    tokens = str2double(split(fgetl(pfile)));
    duct(i).nicells = tokens(1);
    duct(i).iicells = tokens(2);
    duct(i).nscells = tokens(3);
    duct(i).iscells = tokens(4);
    duct(i).nlsegs = tokens(5);
    duct(i).ilsegs = tokens(6);
  end
  
  % get lumen node data
  lnodes = zeros(nlnode,3);
  lradii = zeros(nlnode,1);
  for i = 1:nlnode
    tokens = str2double(split(fgetl(pfile)));
    lnodes(i,:) = tokens(1:3);
    lradii(i) = tokens(4);
  end
  
  % get lumen segment data
  lsegs = zeros(nlseg,2);
  for i = 1:nlseg
    tokens = str2double(split(fgetl(pfile)));
    lsegs(i,:) = tokens(1:2);
  end
  
  % get the cell info
  cells = struct([]);
  for i = 1:ncell
    tokens = str2double(split(fgetl(pfile)));
    cells(i).nverts = tokens(1);
    cells(i).nfaces = tokens(3);
    cells(i).ntets = tokens(5);
  end
  
  % get the vertex data
  verts = zeros(nvert,3);
  for i = 1:nvert
    verts(i,:) = str2double(split(fgetl(pfile)));
  end
  
  % get the face data
  faces = cell(1,ncell);
  for i = 1:ncell
    faces{i} = zeros(cells(i).nfaces,3);
    for j = 1:cells(i).nfaces
      faces{i}(j,:) = str2double(split(fgetl(pfile)));
    end
  end
  
  % get the tetrahdron data
  tets = cell(1,ncell);
  for i = 1:ncell
    tets{i} = zeros(cells(i).ntets,4);
    for j = 1:cells(i).ntets
      tets{i}(j,:) = str2double(split(fgetl(pfile)));
    end
  end
  
  fclose(pfile);
end

function [count] = get_count(pfile, value)
  while 1
    tokens = split(fgetl(pfile));
    if strcmp(tokens{1},'element') && strcmp(tokens{2},value)
      count = str2double(tokens{3});
      break;
    end
  end
end

function [] = skip_header(pfile)
  while 1
    if(strcmp(fgetl(pfile),'end_header'))
      break;
    end
  end
end

function v = get_version(pfile, version)
  while 1
    tokens = split(fgetl(pfile));
    if strcmp(tokens{1},'comment')
      v = tokens{5};
      break;
    end
  end

end

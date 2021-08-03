
clear all
close all
clc


load output_KCa3_VPLC0.0025.mat

%load new_saliva_output.mat
 
 
figure(1)
%trisurf(surftrilist,p(:,1),p(:,2),p(:,3),'FaceColor','cyan')

trisurf(apicaltrilist,p(:,1),p(:,2),p(:,3),'FaceColor','red','FaceAlpha',1)
hold on
trisurf(basaltrilist,p(:,1),p(:,2),p(:,3),'FaceColor','blue','FaceAlpha',0.1)

 rando  = randi(size(p,1),20,1);                                            % 20 random grid points 
 apical = find(dist_to_apical<0.5);                                         % apical grid points 
 middle = find(dist_to_apical>0.5 & dist_to_apical <2);                     % middle grid points 
 basal  = find(dist_to_apical>3);                                           % basal grid points

% pull out the points along a line through the cell

     % coordinates of mean apical and mean basal points
     aa = [-6.77,17.35 -0.23];                                  % just by trying a bunch of choices. Not a good method!
     bb = mean(p(basal,:));

     nn = bb-aa; nn = nn/norm(nn);                              % unit vector in direction of apical-to-basal line
     perp_dist=zeros(np,1);
     for i=1:np
        perp_dist(i) = norm(p(i,:)-aa - dot(p(i,:)-aa,nn)*nn);  % formula for dist of point to line
     end
     linescan = find(perp_dist<0.32);                           % indices of all grid points close enough to line

     % Now put blobs all along the line, just to see where they are
         X=p(linescan,1);
         Y=p(linescan,2);
         Z=p(linescan,3);
         scatter3(X,Y,Z,100,'g', ...
             'MarkerEdgeColor','k',...
             'MarkerFaceColor','g')
         
%      % Now put blobs in all apical region nodes, just to see where they are
%          X=p(apical,1);
%          Y=p(apical,2);
%          Z=p(apical,3);
%          scatter3(X,Y,Z,100,'g', ...
%              'MarkerEdgeColor','k',...
%              'MarkerFaceColor','b')
hold off

c_ls = (sol(linescan,:));                                       % the solution at each grid point along the line scan.
ip_ls = (sol(linescan+np,:));

figure(3)
plot(tim,c_ls)

figure(4)
plot(tim,mean(sol(basal,:)))
%ylim([0.0705 0.0712])

%% Plotting secretion stuff

        Nal 	= SSsol(1,:);
        Kl 		= SSsol(2,:);
        Cll     = SSsol(3,:);
        w       = SSsol(4,:);
        Na 		= SSsol(5,:);
        K 		= SSsol(6,:);
        Cl      = SSsol(7,:);
        HCO     = SSsol(8,:);
        H 		= SSsol(9,:);
        Va      = SSsol(10,:);
        Vb      = SSsol(11,:);
              
        Qa =  par.La*0.9 * ( 2 * ( Nal + Kl - Na - K - H ) - par.CO20 + par.Ul);  
        Qt =  par.Lt * ( 2 * ( Nal + Kl ) + par.Ul - par.Ie );
        Qtot=(Qa+Qt);
        
        figure(10)

        subplot(3,3,1)
        plot(tim,Qtot,'LineWidth',2)
        ylabel('fluid flow')
        
        subplot(3,3,2)
        plot(tim,mean(sol(apical,:)),'LineWidth',2)
        ylabel('mean apical calcium')

        subplot(3,3,3)
        plot(tim,w,'LineWidth',2)
        ylabel('cell volume')

        subplot(3,3,4)
        plot(tim,Nal,'LineWidth',2)
        ylabel('Na (lumen)')
        %ylim([0 140])

        subplot(3,3,5)
        plot(tim,Kl,'LineWidth',2)
        ylabel('K (lumen)')
        %ylim([0 140])

        subplot(3,3,6)
        plot(tim,Cll,'LineWidth',2)
        ylabel('Cl (lumen)')
        %ylim([0 140])

        subplot(3,3,7)
        plot(tim,Va,'LineWidth',2)
        ylabel('Va')
        %ylim([-60 -20])

        subplot(3,3,8)
        plot(tim,Vb,'LineWidth',2)
        %ylim([-60 -20])
        ylabel('Vb')

        subplot(3,3,9)
        plot(tim,Cl,'LineWidth',2)
        ylabel('Cl')

        set(gcf,'Position',[1400,100,1000,1200])
        
%         figure(11)
%         VtK = par.RTF * log( Kl / par.Ke );  
%         VtNa = par.RTF * log( Nal / par.Nae ); % tight junction K
%         Vt = Va - Vb;
%         JtK  = par.GtK * par.St * ( Vt - VtK ) / par.F;
%         JtNa = par.GtNa * par.St * ( Vt - VtNa ) / par.F;
%         plot(t,JtK,t,JtNa,'LineWidth',2)

%         figure(12)
%         plot(tim,Ip,'LineWidth',2)
%         ylabel('IP_3')  
% ----------------------------------------------------------


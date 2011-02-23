function [chosen,ACTIVITY,MEANDIST] = collicular_positions_two_thresholds(DATA,N,ELL_PARAMS,THRESH0,SIGMA,RANSTART) 
% to produce N evenly distributed nodes within an ellipse of collliculus
% which are not too close and where there is 
% sufficient signal in the ELEVATION and AZIMUTHAL components of amplitude
% AND at least 33% of the pixels within the circle of radius SIGMA 
%surrounding each node are active


% minimal spacing set automatically

%                 Loading the raw data

[ss Azim_amp Elev_amp Azim_phase Elev_phase MapAnsColor] = load_data(DATA);

colormap(MapAnsColor);

%   REMEMBER THAT ORDER OF INDICES IN AZIM_AMP etc IS (Y,X)


%                First of all define the elliptical area of colliculus
%                within which the mean amplitude is sought

RA = ELL_PARAMS(1);
RB = ELL_PARAMS(2);
ANG = ELL_PARAMS(3);
X0 = ELL_PARAMS(4);
Y0 = ELL_PARAMS(5);

ACC = 150;
[ELEV_AMP, AZIM_AMP]= select_pixels(ELL_PARAMS,Elev_amp,Azim_amp);

III = find(ELEV_AMP > 0);
p_i_e = length(III);

%---------------------------------------------------------------------------------------------

%  FOR BOTH  ELEV AND AZIM 
%    FIND MEAN ACTIVITY OVER ALL THE 250 X 250 pixels

global_mean_E = mean(mean(Elev_amp));
global_mean_A = mean(mean(Azim_amp));


%         THEN AVERAGE OVER NON ZERO VALUES OF ACTIVITY

E_AMP = reshape(ELEV_AMP,62500,1);
II = find(E_AMP >0);
mean_E = mean(E_AMP(II));


A_AMP = reshape(AZIM_AMP,62500,1);
II = find(A_AMP >0);
mean_A = mean(A_AMP(II));

%        SET THRESHOLDS CALCULATED ACCORDING TO MEANS WITHIN ELLIPSE

THRESH_E = THRESH0(1)*mean_E;
THRESH_A = THRESH0(2)*mean_A;


%-----------------------------------------------------------------------------
%        FIND PIXELS THAT HAVE ENOUGH AZIM and ELEV ACTIVITY

[ELEV_AMP, AZIM_AMP]= select_pixels(ELL_PARAMS,Elev_amp,Azim_amp);

ACTIVITY = ELEV_AMP > THRESH_E & AZIM_AMP > THRESH_A;

II = find (ACTIVITY ==1);
ellipse_active = length(II)

ELEV_AMP = ELEV_AMP.*ACTIVITY;

AZIM_AMP = AZIM_AMP.*ACTIVITY;

figure(99)
clf

%-------------------------------------------------------------------

%                     PLOT OUT ACTIVITY PATTERNS
subplot(3,2,3)

imagesc(Elev_amp);

hold on

[h X_ELL Y_ELL] = ellipse(RA,RB,ANG,X0,Y0,'k');

title({['Raw Elev Amp, mean: ',num2str(global_mean_E,3)];['Mean within ellipse: ',num2str(mean_E,3),' (',num2str(p_i_e),' pixels)']});


xlabel('X');
ylabel('Y');

axis ij 

axis image

hold off


subplot(3,2,4)
   
imagesc(Azim_amp);

hold on

[h X_ELL Y_ELL] = ellipse(RA,RB,ANG,X0,Y0,'k');

title({['Raw Azim Amp, mean: ',num2str(global_mean_A,3)];['Mean within ellipse: ',num2str(mean_A,3),' (',num2str(p_i_e),' pixels)']});

xlabel('X');
ylabel('Y');

axis image

axis ij

hold off


subplot(3,2,5)

imagesc(ELEV_AMP);

hold on

[h X_ELL Y_ELL] = ellipse(RA,RB,ANG,X0,Y0,'k');

title({['Elev amp of ',num2str(ellipse_active),' pixels within ellipse which are'];['above Elev, Amp thresholds: ',num2str(THRESH_E,3),' and ',num2str(THRESH_A,3)]});

axis ij

axis image

hold off

subplot(3,2,6)

imagesc(AZIM_AMP);

hold on

[h X_ELL Y_ELL] = ellipse(RA,RB,ANG,X0,Y0,'k');

title({['Azim amp of ',num2str(ellipse_active),'  pixels within ellipse which are'];['above Elev, Amp thresholds: ',num2str(THRESH_E,3),' and ',num2str(THRESH_A,3)]});

axis ij

axis image

hold off


% CHOOSE THE RECORDING POINTS

W1 = floor(min(X_ELL));
W2 = ceil(max(X_ELL));
W3 = floor(min(Y_ELL));
W4 = ceil(max(Y_ELL));

XLIST = [W1:W2];
XW = length(XLIST);

YLIST = [W3:W4];
YW = length(YLIST);

rand('twister', RANSTART);

MINSPACING = 0.75*sqrt(pi*RA*RB/N)

NLIM = N*100;

finished =0;
while finished == 0
      chosen = zeros(N,2);
      MINSPACING = 0.95*MINSPACING
      MSPACING = round(MINSPACING);
      CANDIDATE = ACTIVITY;

      n=1;
      ntry = 1;
      min_points=pi*SIGMA^2/3;    % 33 % of sampling area to be filled
      while n<=N & ntry <=NLIM & finished ==0
           SORTLISTX = randperm(XW);
            X = XLIST(SORTLISTX(1));
            SORTLISTY = randperm(YW);
            Y = YLIST(SORTLISTY(1));
            ntry=ntry+1;
            num_candidates = 0; %keep track of the number of ACTIVE points
%               num_poss = 0; %count all points within radius SIGMA
            for XX = X-SIGMA:X+SIGMA
                  for YY = Y-SIGMA:Y+SIGMA
                       if sqrt((XX-X)^2 + (YY-Y)^2) < SIGMA
                          num_candidates = num_candidates +CANDIDATE(YY,XX);
%                          num_poss = num_poss +1;
                       end
                  end
	        end
%            if CANDIDATE(Y,X) == 1 &  num_candidates > num_posS/
            if CANDIDATE(Y,X) == 1 &  num_candidates > min_points
                chosen(n,1) = X;
                chosen(n,2) = Y;
                num_can(n)=num_candidates;
                if n==N 
                    finished = 1; 
                end
                n=n+1;
                for XX = X-MSPACING:X+MSPACING
                    for YY = Y-MSPACING:Y+MSPACING
                        if sqrt((XX-X)^2 + (YY-Y)^2) < MINSPACING
                            CANDIDATE(YY,XX)=0;
                        end
                    end
                end
            end
      end
end

num_can

%                      reorder collicular points

[DUMMY ISORT] = sort(chosen(:,1));
chosen = chosen(ISORT,:)


%                 FIND MEAN DIST BETWEEN  NEAREST-NEIGHBOURS  ***NOT USED**

%       BY DOING TRIANGULATION ON COLLICULUS AND THEN MEASURING NN LENGTHS

%[ttriangles, tneighbours] = triangulation(chosen,ACC)
%numb_ttriangles = size(ttriangles,1);
%numb_points = size(chosen,1)

%ll = 1;
%for i=1:numb_points-1
%    for j=i+1:numb_points

%	 if tneighbours(i,j) == 1
%	   lengths(ll) = dist(chosen(i,:),chosen(j,:)');
%	   ll=ll+1;
%	end
%    end
%end
%MEANDIST = mean(lengths);
%STDDIST = std(lengths);
%disp(['Number of trials needed: ',num2str(ntry);]);


%            THIS CALCULATES THE MEAN MINIMUM DISTANCE BETWEEN NEIGHBOURS
%          THIS METHOD PREFERRED TO CALCULATING THE MEAN DISTANCE OVERALL 
%          TO EMPHASISE THAT THE PROXIMITY TO AT LEAST ONE NEIGHBOUR IS MOST IMPORTANT

dd = 10000*eye(N)+ (dist(chosen'));

MEANDIST = mean(min(dd))
STDDIST = std(min(dd));
disp(['Number of trials needed: ',num2str(ntry);]);



subplot(3,2,1)

hold on  

plot(chosen(:,1),chosen(:,2),'*');

title([num2str(N), ' chosen points within the ellipse.']);

xmax = W2;
xmin = W1;

ymax = W4;
ymin = W3;

rectangle('Position',[W1 W3  (W2-W1) (W4-W3)],'EdgeColor',[0 0 0]);

[h X_ELL Y_ELL] = ellipse(RA,RB,ANG,X0,Y0,'k');

hold off
axis([xmin xmax ymin ymax]);

axis equal

axis ij

axis image

subplot(3,2,2)


hold on

text(-0.1,1.0,[num2str(date),'        ',ss],'Color','b');

text(-0.1,0.9,['COLLICULAR\_POSITIONS\_TWO\_THRESHOLDS.M']);

text(-0.1,0.8,['At least 33% of pixels in the sampling circle must be active']);
text(-0.1,0.6,[num2str(N),' CHOSEN POINTS INSIDE AN ELLIPSE'],'Color','r');

text(-0.1,0.5,['MEAN, COEFFVAR OF MIN NN DISTS: ',num2str(MEANDIST,3),' pixels,',num2str(100*STDDIST/MEANDIST,2),' %'],'Color','r');

text(-0.1,0.4,['DATA, N, ELL\_PARAMS, THRESHS, SIGMA, RANSTART']);
text(-0.1,0.3,[num2str(DATA),',',num2str(N),',',num2str(RA),',',num2str(RB),',',num2str(ANG*180/pi,3),',',num2str(X0),',',num2str(Y0),',',num2str(THRESH0(1)),',',num2str(THRESH0(2)),',',num2str(SIGMA),',',num2str(RANSTART)]);

text(-0.1,0.2,['Number [total] active pixels in ellipse: ',num2str(ellipse_active),'[',num2str(p_i_e),']'],'Color','r');

text(-0.1,0.05,['Thresholds: fractions ',num2str(THRESH0(1)),',',num2str(THRESH0(2)),'  of mean activity levels']);
axis off

colormap(MapAnsColor);

orient tall

filename = [num2str(DATA),'_fred99.eps'];
print(99, '-depsc2', filename)

clear

InputDir = uigetdir
cd(InputDir)

%[fileName_in] = uigetfile();
[fileName_in] = uigetfile('*.tif');
Im_1 = imread(fileName_in);

Im_2 = imtool(Im_1);

% set size in pixels, define square region size for each sample in plate layout
SizeRegionPerWell_x = 50;
SizeRegionPerWell_y = 50;

% Im_1 = imread('1.tif');
% Im_1 = imread('0000072_02_800 16b.tif');

% convert to handle grayscale RGB input images
if size(Im_1, 3) > 1
   Im_1 = Im_1(:,:,1); 
end

% find im values for display
Im_1_PixInt_Min = min(Im_1(:));
Im_1_PixInt_Max = max(Im_1(:));

Im_1_rsrv = Im_1;




% user selects points
h_f1 = figure, imshow(Im_1,[]), title ('adjust image contrast then close contrast control window')

ImContrast_temp = imcontrast(h_f1);
ContrastMax = str2num(ImContrast_temp.Children(1).Children(3).Children.Children(2).Children.Children(2).Children(2).String)
ContrastMin = str2num(ImContrast_temp.Children(1).Children(3).Children.Children(2).Children.Children(2).Children(5).String)

uiwait(ImContrast_temp)

%%% select positions to orient sample boundaries within plate
uiwait(msgbox('select center of upper left well with mouse then hit enter', 'modal'))
[UL_x,UL_y] = getpts(h_f1)

uiwait(msgbox('select center of upper right well with mouse then hit enter', 'modal'))
[UR_x, UR_y] = getpts(h_f1)

uiwait(msgbox('select center of lower right well with mouse then hit enter', 'modal'))
[LR_x, LR_y] = getpts(h_f1)

% plot points for confirmation
Im_1_points = false(size(Im_1));
Im_1_points(uint16(UL_y),uint16(UL_x))=1;
Im_1_points(uint16(UR_y),uint16(UR_x))=1;
Im_1_points(uint16(LR_y),uint16(LR_x))=1;
% figure, imshow(Im_1_points,[]), title('selected corner points')

%select plate layout
PlateChoice_List = {'96 well','384 well','1536 well'};
[PlateChoice_Index,PlateChoice_tf] = listdlg('ListString',PlateChoice_List,'PromptString',{'Select your plate type.',''},'SelectionMode','single');

% NumTotalCol = 24;
% NumTotalRow = 16;

NumTotalCol_List = [12; 24; 48];

NumTotalCol_First = 1;
NumTotalCol_Last = NumTotalCol_List(PlateChoice_Index);

NumTotalRow_List = [8; 16; 32];

NumTotalRow_First = 1;
NumTotalRow_Last = NumTotalRow_List(PlateChoice_Index);


prompt = {'Enter first row coordinate selected:','Enter first column coordinate selected:','Enter last row coordinate selected:','Enter last column coordinate selected:'};
dlgtitle = 'Input';
dims = [1 35];
definput = {num2str(NumTotalRow_First), num2str(NumTotalCol_First), num2str(NumTotalRow_Last), num2str(NumTotalCol_Last)};
answer_PlateCoord = inputdlg(prompt,dlgtitle,dims,definput);

NumTotalRow_First = str2num(cell2mat(answer_PlateCoord(1)));
NumTotalCol_First = str2num(cell2mat(answer_PlateCoord(2)));
NumTotalRow_Last = str2num(cell2mat(answer_PlateCoord(3)));
NumTotalCol_Last = str2num(cell2mat(answer_PlateCoord(4)));


Mask_Outline_Single = false(SizeRegionPerWell_y+1,SizeRegionPerWell_x+1);
Mask_Outline_Single(:,1)=1;
Mask_Outline_Single(:,end)=1;
Mask_Outline_Single(1,:)=1;
Mask_Outline_Single(end,:)=1; 

    % for single time point analysis
    Count_Time = 0;
    Count_Time = Count_Time + 1;
    
    Count_Time
    
    Output_PlateLayout_1 = uint16(zeros(NumTotalRow_List(PlateChoice_Index),NumTotalCol_List(PlateChoice_Index))); 
    
    
 % myfilename = sprintf('%d.tif', Count_Time);
  
 % Im_1 = imread(myfilename);
  
    if size(Im_1,3) > 1
        
        Im_1 = Im_1(:,:,1);
    end

        Count_Output = 0;
        
        Im_1_points_dil = imdilate(Im_1_points, strel('disk',1));
% figure, imshow(Im_1_points_dil,[])

% y scale is backward becuase it is an image, flips sign of slope
Line_UL_UR_slope = (UL_y - UR_y) / (UL_x - UR_x);
Line_UL_UR_degrees = atand(Line_UL_UR_slope);
    
Im_2 = imrotate(Im_1, Line_UL_UR_degrees);
% figure, imshow(Im_2,[Im_1_PixInt_Min Im_1_PixInt_Max])

test1 = Im_2;

Im_2_points_dil = imrotate(Im_1_points_dil, Line_UL_UR_degrees);
% figure, imshow(Im_2_points_dil,[])

Im_2_points_dil_L = bwlabel(Im_2_points_dil);

Im_2_points_props = regionprops(Im_2_points_dil_L, 'Centroid');



% the below values should be left as floating point format... 
% warning: loss of precision. should convert to integer at last step

Rotated_WellCenters_XY = reshape([Im_2_points_props.Centroid],[2,size(Im_2_points_props,1)])'; % correct matrix fomat for paired xy points


Rotated_WellCenters_PerPlate_Xmin = min(Rotated_WellCenters_XY(:,1));
Rotated_WellCenters_PerPlate_Xmax = max(Rotated_WellCenters_XY(:,1));

Rotated_WellCenters_PerPlate_Xdist = Rotated_WellCenters_PerPlate_Xmax - Rotated_WellCenters_PerPlate_Xmin;

Rotated_WellDistance_X = Rotated_WellCenters_PerPlate_Xdist / ((NumTotalCol_Last - NumTotalCol_First));

Rotated_WellCenters_PerPlate_Ymin = min(Rotated_WellCenters_XY(:,2));
Rotated_WellCenters_PerPlate_Ymax = max(Rotated_WellCenters_XY(:,2));


Rotated_WellCenters_PerPlate_Ydist = Rotated_WellCenters_PerPlate_Ymax - Rotated_WellCenters_PerPlate_Ymin;

Rotated_WellDistance_Y = Rotated_WellCenters_PerPlate_Ydist / ((NumTotalRow_Last - NumTotalRow_First))


SizeRegionPerWell_x_HalfDist = (SizeRegionPerWell_x/2);
SizeRegionPerWell_y_HalfDist = (SizeRegionPerWell_y/2);

          
Mask_Outline_All = uint8(zeros(size(Im_2)));


CountRow = 0;


    for C1 = NumTotalRow_First:NumTotalRow_Last
        % C1 = 1
        % C1
        CountRow = CountRow + 1;
        
        CountCol = 0;
        for C2 = NumTotalCol_First:NumTotalCol_Last
          % C2 = 10
          % C2
          CountCol = CountCol + 1;
          
          

          SingleWell_Center_X = uint16(((CountCol-1)*Rotated_WellDistance_X) + Rotated_WellCenters_PerPlate_Xmin);
          SingleWell_Center_Y = uint16(((CountRow-1)*Rotated_WellDistance_Y) + Rotated_WellCenters_PerPlate_Ymin);
                
          % test1 = Im_2;
          test1(uint16(SingleWell_Center_Y), uint16(SingleWell_Center_X)) = 40000;

          % 2018b
          % h_ROI = images.roi.Rectangle(gca,'Position',[(SingleWell_Center_X - SizeRegionPerWell_x_HalfDist),(SingleWell_Center_Y - SizeRegionPerWell_y_HalfDist),SizeRegionPerWell_x,SizeRegionPerWell_y],'StripeColor','r');

          Im_2_SingleWell = Im_2(uint16((SingleWell_Center_Y - SizeRegionPerWell_y_HalfDist)):uint16((SingleWell_Center_Y + SizeRegionPerWell_y_HalfDist)), uint16((SingleWell_Center_X - SizeRegionPerWell_x_HalfDist)):uint16((SingleWell_Center_X + SizeRegionPerWell_x_HalfDist)));
          % figure, imshow(Im_2_SingleWell,[]), title(strcat('r', num2str(CountRow), ' ', 'c', num2str(CountCol)))

          Im_2_SingleWell_MedFilt = medfilt2(Im_2_SingleWell,[5 5], 'symmetric');
           % figure, imshow(Im_2_SingleWell_MedFilt,[])

          [Thrsh,Thrsh_EM] = multithresh(Im_2_SingleWell_MedFilt,1);
          
          % select lowest 5% of pixels, then measure median to estimate
          % local background
          BG_x1 = sort(Im_2_SingleWell_MedFilt);
          BG_x2 = BG_x1(1:uint16(size(BG_x1,1)/20));
          BG_x3 = median(BG_x2(:));

          Count_Output = Count_Output + 1;
          
          if Count_Time == 1 % set this on first time point only
              BG_Well_perTime(Count_Output,1) = BG_x3;
          end
          
           if Thrsh_EM > 0.3 
          Mask = imbinarize(Im_2_SingleWell_MedFilt, double(Thrsh)/256);
          % figure, imshow(Mask,[])

         
         
          Mask_L = bwlabel(Mask);
          % figure, imshow(Mask_L,[])
          
           Mask_L_Stats = regionprops(Mask_L,'Area', 'Circularity', 'Solidity');
          % Mask_L_Stats.Area, Mask_L_Stats.Circularity, Mask_L_Stats.Solidity
          
          % below line selects correct regions within wells for analysis
          % round well shapes of correct size, adjust values below if
          % needed
           ROI_idx = find([Mask_L_Stats.Area] > 50 & [Mask_L_Stats.Area] < 500 & [Mask_L_Stats.Solidity] > 0.8); 
             Mask_2 = ismember(Mask_L, ROI_idx);
%             Mask_2 = false(size(Mask));
%             Mask_2(Mask > 0) = 1;
            
           
          

          Mask_MeanInt_Uncorrected = mean(Im_2_SingleWell(Mask_2==1));
          BG_MeanInt = mean(Im_2_SingleWell(Mask_2==0));
          Mask_MeanInt_SubtractBG = Mask_MeanInt_Uncorrected - BG_MeanInt;

          Mask_Area = size(Mask_2(Mask_2==1),1);
         
        % Thrsh_EM_List(Count_Well_Per_Time,1) = Thrsh_EM;

         OutputNumbers(Count_Output,1) = C1; % well row coord
         OutputNumbers(Count_Output,2) = C2; % well column coord
         OutputNumbers(Count_Output,3) = Mask_Area;
         OutputNumbers(Count_Output,4) = Mask_MeanInt_Uncorrected;
         OutputNumbers(Count_Output,5) = BG_MeanInt;
         OutputNumbers(Count_Output,6) = Mask_MeanInt_SubtractBG;
         
         OutputNumbers(Count_Output,7) = Mask_MeanInt_SubtractBG * Mask_Area; % integrated intenity of spot, BG subtracted

         OutputNumbers_MeanInt_Uncorrected(Count_Output,Count_Time) = mean(Im_2_SingleWell(:));
         OutputNumbers_MeanInt_SubtractBG(Count_Output,Count_Time) = Mask_MeanInt_SubtractBG;
         OutputNumbers_MeanInt_SubtractBG(Count_Output,Count_Time) = mean(Im_2_SingleWell(:))-BG_x3;

        OutputNumbers_BG(Count_Output,Count_Time)=BG_x3;
        
         Mask_Outline_Single = bwperim(Mask_2);
         % Mask_Outline_Single = true(size(Im_2_SingleWell));
        
        Mask_Outline_Single_FullField = false(size(Im_2)); 
        Mask_Outline_Single_FullField(((SingleWell_Center_Y - SizeRegionPerWell_y_HalfDist)):uint16((SingleWell_Center_Y + SizeRegionPerWell_y_HalfDist)), uint16((SingleWell_Center_X - SizeRegionPerWell_x_HalfDist)):uint16((SingleWell_Center_X + SizeRegionPerWell_x_HalfDist))) = Mask_Outline_Single;
        
        Mask_Outline_All(Mask_Outline_Single_FullField == 1) = 200;
        
           else % conditional for threshold quality
          
              % conditional
              % zeros are deployed when no spot is reliably segmented    
             OutputNumbers(Count_Output,1) = C1; % well row coord
             OutputNumbers(Count_Output,2) = C2; % well column coord
             OutputNumbers(Count_Output,3) = 0 % Mask_Area;
             OutputNumbers(Count_Output,4) = 0 %Mask_MeanInt_Uncorrected;
             OutputNumbers(Count_Output,5) = 0 %BG_MeanInt;
             OutputNumbers(Count_Output,6) = 0 %Mask_MeanInt_SubtractBG;

             OutputNumbers(Count_Output,7) =  0 % Mask_MeanInt_SubtractBG * Mask_Area; % integrated intenity of spot, BG subtracted

             OutputNumbers_MeanInt_Uncorrected(Count_Output,Count_Time) = 0 % mean(Im_2_SingleWell(:));
             OutputNumbers_MeanInt_SubtractBG(Count_Output,Count_Time) = 0 % Mask_MeanInt_SubtractBG;
             OutputNumbers_MeanInt_SubtractBG(Count_Output,Count_Time) = 0 % mean(Im_2_SingleWell(:))-BG_x3;

            OutputNumbers_BG(Count_Output,Count_Time)= 0 %BG_x3;
              
          end % end conditional for threshold quality
        
        

        end % end C2 Column Loop
    end % C1 Row Loop 
    
    % figure, imshow(Mask_Outline_All,[])
    
    test1 = imdilate(test1, strel('disk',3));

    % imshow(test1,[min(Im_1(:)) max(Im_1(:))]), title(strcat('TimePoint ',num2str(Count_Time)))
    
    ScaleInt =  double(max(Im_2(:))-median(BG_Well_perTime(:)))/256;
    Im_2_NoBG = Im_2 - double(median(BG_Well_perTime(:)));
    Im_2_NoBG_8bit = uint8(Im_2_NoBG/ScaleInt);
    % figure, imshow(Im_2_NoBG_8bit,[])
    
    
    Im_2_RegionOverlay = uint8(Mask_Outline_All);
    Im_2_RegionOverlay(:,:,2)=Im_2_NoBG_8bit;
    Im_2_RegionOverlay(:,:,3) = 0;
      figure, imshow(Im_2_RegionOverlay,[])

 OutputNumbers_SubtractSingleBG = double(OutputNumbers_MeanInt_Uncorrected) - double(median(BG_Well_perTime(:)));
 
 % OutputList =  OutputNumbers(:,6); mean int minus local background
OutputList =OutputNumbers_SubtractSingleBG; % mean int minus average global bg

Concat_RowCol_and_Output = [OutputNumbers(:,1:2),OutputList];

% save output as excel file in same directory as input files
Filename_Out = strcat(fileName_in,'_MatlabResults_1.xlsx');
writematrix(Concat_RowCol_and_Output, Filename_Out,'FileType','spreadsheet')



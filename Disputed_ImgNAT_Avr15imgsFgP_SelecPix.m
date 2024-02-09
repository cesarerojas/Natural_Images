function [] = Disputed_ImgNAT_Avr15imgsFgP_SelecPix(RESULTS_folder, numCaptDevice, numInitialImg, width, height)
%Function that analyses natural images.
% -------------------------------------------------------------------------
% Copyright (c) 2024 Instituto Politécnico Nacional (IPN), México.
% All Rights Reserved.
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Permission to use, copy, modify, and distribute this software for
% educational, research, and non-profit purposes is hereby granted, without 
% fee or written agreement is hereby granted, provided that this copyright 
% notice appears in all copies. IPN does not warrant that the operation of the 
% program will be uninterrupted or error-free. The end user understands that 
% the program has been developed for research purposes and is advised not to
% rely exclusively on the program for any reason. In no even shall IPN be 
% liable to any party for any direct, indirect, special, incidental, or 
% consequential damages, including lost profits, arising out of the use of this
% software. IPN disclaims all warranties, and has no obligation to provide
% maintenance, support, updates, improvements, or modifications.
% -------------------------------------------------------------------------
% Version: 20240208
% -------------------------------------------------------------------------
% Authors:
%      César Enrique Rojas-López,     IPN-ESIME Culhuacan.
%      Omar Jiménez-Ramírez,          IPN-ESIME Culhuacan.
%      Luis Niño-de-Rivera-Oyarzabal, IPN-ESIME Culhuacan.
%      Leonardo Palacios-Luengas,     UAM-Iztapalapa.   
%      Rubén Vázquez-Medina,          IPN-CICATA Querétaro.

% Contact: rvazquez@ipn.mx | February 2024
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% This function estimates the Mahalanobis distance between the intrisic signal left by the capture device 
% in natural digital images and different camera fingerprints
%
% INPUTS:
%     RESULTS_folder     - Full path where the results will be stored
%     path_DisputedImg   - Full path to where the disputed images will be read.
%                          This variable will be requested 10 times since 10 disputed images 
%                          will be analyzed for each capture device. It is
%                          assumed that each folder of disputed image will store the d_k images of all the capture devices.
%                          Example:
%                             01 iPhone_SE2020_1_nat_41.JPG
%                             02 iPhone_XR_nat_41.JPG
%                             03 Motorola_G4Plus_nat_41.jpg
%                             04 Samsung_GalaxyA01_nat_41.jpg
%                             05 Samsung_GalaxyNote9_nat_41.jpg
%                             06 Motorola_G20_nat_41.jpg
%                             07 iPhone_SE2020_2_nat_41.JPG
%                             08 Huawei_Y9-2019_nat_41.jpg
%     numCaptDevice      - Number of capture devices considered in the analysis
%     numInitialImg      - Consecutive number from which the disputed images 
%                          will be considered in each capture device folder
%     width              - Width of the image clipping extracted from reference flat images,
%                          considering as reference point the image centroid.
%     height             - Heigth of the image clipping extracted from reference flat images,
%                          considering as reference point the image centroid.
%
% OUTPUTS:
%     AverageImgs        - Matrix with a number of rows equal to width x height and a number 
%                          of columns equal to 16. The odd columns contain the pixel intensity 
%                          of the green layer of disputed image and the even columns contain the PRNU.
%     total              - mat file with an 8 x 16 matrix. In each pair of columns, the first columnn refers 
%                          to the number of times the inequality DM(T_{i,j},H_{01-15,k}) ≤ DM(T_{i,j},H_{01-15,l}
%                          was satisfied, and the second column refers to when it was not satisfied, where T_{i,j}
%                          is the intrinsic signal left by the j-th capture device on the disputed image number i, 
%                          and H_{01-15,l} is the fingerprint for l-th capture device using fifteen reference flat image. 
%                          l is the first index, which varies by affecting the rows 
%                          in the tables, and k is the second index, which varies by affecting each pair of columns
%                          Thus, this file shows the 8 results of a disputed image when it is assumed 
%                          to belong to a different capture device.
%                           
%
% REQUIRED FUNCTION:
%    Average_15RefImgs.m - It calculates the camera fingerprint using fifteen reference flat images.
%
%
% PROCEDURE:
% 1.- With the function "Average_15RefImgs" 15 reference images of each capture device 
%     are read to obtain their average of the icero and PRNU. 
%     The number of the capture device is sent, as well as the width and height required.
% 2.- Mean, sigma and delta are obtained from the first average capture device to calculate upper and lower limits. 
% 3.- The first disputed image is read out and the vectors of icero and PRNU are obtained, T.
% 4.- Icero and PRNU pixels of the disputed image are limited.
% 5.- The pixel icero and PRNU of the averaged reference images are read, R. (see 1)
% 6.- T and H are made the same size.
% 7.- The Mahalanobis Distance, DM, is obtained for each capture device. 
% 8.- The Mahalanobis distance is calculated step by step. 
% 9.- The results are saved.

close all 
clc

AverageImgs = Average_15RefImgs(numCaptDevice, width, height);

while numInitialImg <= 50  % Only ten images are considered in the analysis.
    icero_y_PRNU = 1;
    DC_References_Limits = 1;
    fprintf('Full path to the disputed image number %d from the eight capture devices: ', numInitialImg)
    path_DisputedImg = input ('', 's');
    while (DC_References_Limits <= numCaptDevice)
        avrg = mean(AverageImgs(:,icero_y_PRNU));
        sigma  = std(AverageImgs(:,icero_y_PRNU));
        delta = sigma*2; % 95.4% of the area is considered under the probability distribution function.
        lower_limit = avrg - delta;
        upper_limit = avrg + delta;
        clear avrg sigma delta width height

        %***DISPUTED IMAGES***
        cd (path_DisputedImg);
        disputedFile = ls(path_DisputedImg);
        numberDisputesFiles = size(disputedFile,1);
        column_counts = 1;

        for i = 3 : numberDisputesFiles
            
            filenameCurrentDisputed = (disputedFile(i,:));
            disputedImage = imreadort(filenameCurrentDisputed);
            I = uint8(double(disputedImage(:,:,2)));

            CV_eta  = NoiseExtractFromImage(I,2);
            Icero_D = double(I)-CV_eta;
            PRNU_D  = WienerInDFT(CV_eta,std2(CV_eta ));

            vector_Icero_D = Icero_D(:);
            vector_PRNU_D  = PRNU_D(:);
            clear CV_eta Icero_D PRNU_D

            %SELECTS PIXELS OF DISPUTED IMAGE THAT COMPLY WITH LIMITS
            V_Icero_and_PRNU = zeros(size(vector_Icero_D,1),1);
            for n=1:size(vector_Icero_D,1)
               if ((vector_Icero_D(n)>=(lower_limit)) && (vector_Icero_D(n)<=(upper_limit)))
                   V_Icero_and_PRNU (n,2) = vector_PRNU_D(n);
                   V_Icero_and_PRNU (n,1) = vector_Icero_D(n);
               end
            end
            clear vector_Icero_D vector_PRNU_D

            %DELETES ROWS THAT DID NOT COMPLY WITH SELECTION LIMITS
            V_Icero_and_PRNU_WITH_ZEROS = V_Icero_and_PRNU;
            indexDelete = find(V_Icero_and_PRNU(:,1)~=0);
            if isempty(indexDelete) == 0
                V_Icero_and_PRNU_WITH_ZEROS = V_Icero_and_PRNU((indexDelete(:)),1:2);
                V_Icero_interest = V_Icero_and_PRNU_WITH_ZEROS(:,1);
                V_PRNU_interest  = V_Icero_and_PRNU_WITH_ZEROS(:,2);

                T = [V_Icero_interest, V_PRNU_interest];
                clear V_Icero_and_PRNU V_Icero_and_PRNU_WITH_ZEROS 
            end

            k = 1;
            icero_y_PRNU2 = 1;
            while (k <= numCaptDevice)
                vector_Icero_Ref (:,1) = AverageImgs(:,icero_y_PRNU2);
                vector_PRNU_Ref (:,1) = AverageImgs(:,icero_y_PRNU2+1);
                H = [vector_Icero_Ref, vector_PRNU_Ref];
                    
                %CROP VECTOR OF SELECTED PIX OF REFERENCE IMAGE TO VECTOR SIZE OF SELECTED PIX OF THE DISPUTED IMAGE
                if length(V_Icero_interest) < length(vector_Icero_Ref)
                    vector_Icero_Ref2 = vector_Icero_Ref(1:length(V_Icero_interest),:);
                    vector_PRNU_Ref2  = vector_PRNU_Ref (1:length(V_Icero_interest),:);
                    H = [vector_Icero_Ref2, vector_PRNU_Ref2];
                    T = [V_Icero_interest, V_PRNU_interest];
                else
                    V_Icero_interes2 = V_Icero_interest(1:length(vector_Icero_Ref),:);
                    V_PRNU_interes2  = V_PRNU_interest(1:length(vector_Icero_Ref),:);
                    T = [V_Icero_interes2, V_PRNU_interes2];
                end
                
                if (k == 1)
                    DM  = zeros (length(T),size(disputedFile,1)-2);
                end
                dif = T-H; COVARIANCE = cov(T-H); COV_INVERSE = inv(COVARIANCE);
                
                for n = 1 : size(T,1)
                    DM(n,k) = sqrt((dif(n,:) * COV_INVERSE) * dif(n,:)');
                end                

                k = k + 1;
                icero_y_PRNU2 = icero_y_PRNU2 + 2;
                clear H
            end
            
            if column_counts == 1
                total = zeros ((numCaptDevice),(numCaptDevice*2));
            end
            C = zeros (size(DM,1),size(DM,2)*3);
            fil = 1;
            for column = 1 : size(DM,2)
                for row = 1 : size(DM,1)
                    if (DM(row,i-2)<DM(row,column))
                        C(row,fil)=1;
                    else
                        if (DM(row,i-2)>DM(row,column))
                            C(row,fil+1)=1;
                        else
                            C(row,fil+2)=1;
                        end
                    end
                end
                fil = fil + 3;
            end

            %GETS COMPARISON TOTALS
            col = 1;

            for h = 1 : size(C,2)/3
                total(h,column_counts) = sum(C(:,col))+sum(C(:,col+2));
                total(h,column_counts+1) = sum(C(:,col+1));
                col = col +3;
            end
            column_counts = column_counts +2; 

            clear DM        
            cd(path_DisputedImg);
        end
        filename  = sprintf('%s.mat', [RESULTS_folder,'\RESUME D',int2str(numInitialImg),', - Selecc_PIX - 2 sigma - Average_15Img assuming belonging DC0',int2str(DC_References_Limits),', - DM 1 a 1 ']);
        save (filename, 'total');
        column_counts = 1;
        clear total 
        DC_References_Limits = DC_References_Limits + 1;
        icero_y_PRNU = icero_y_PRNU + 2;
    end
    numInitialImg = numInitialImg + 1;
end
end
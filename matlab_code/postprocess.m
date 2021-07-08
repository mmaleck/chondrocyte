% plot the solution

fileID = fopen('IV_I_K_Ca.txt','r');
formatSpec = '%f';
sizeA = [2 Inf];
I_K_Ca = fscanf(fileID,formatSpec, sizeA);
fclose(fileID);

subplot(4,4,1)
plot(I_K_Ca(1,:), I_K_Ca(2,:))
title('I_{KCa}')

fileID = fopen('IV_act_DR.txt','r');
formatSpec = '%f';
sizeA = [2 Inf];
act_DR = fscanf(fileID,formatSpec, sizeA);
fclose(fileID);

subplot(4,4,2)
plot(actDR(1,:), act_DR(2,:))
title('act_{DR}')

fileID = fopen('IV_I_bq.txt','r');
formatSpec = '%f';
sizeA = [2 Inf];
IV_I_bq = fscanf(fileID,formatSpec, sizeA);
fclose(fileID);

subplot(4,4,3)
plot(IV_I_bq(1,:), IV_I_bq(2,:))
title('I_{bq}')

fileID = fopen('IV_I_K_DR.txt','r');
formatSpec = '%f';
sizeA = [2 Inf];
IV_I_K_DR = fscanf(fileID,formatSpec, sizeA);
fclose(fileID);

subplot(4,4,3)
plot(IV_I_bq(1,:), IV_I_bq(2,:))
title('I_{bq}')


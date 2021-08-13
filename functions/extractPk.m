function [Eavg,Iavg,E_nullavg,I_nullavg] = extractPk(tempMat,adjInd_E,adjInd_I,tempMat_base,E_nullLoc,I_nullLoc) 
% This function accepts an input matrix (observations x time) and arrays of
% max and min peak locations to better estimate max and min amplitudes by
% taking average of -0.5ms to 0.5ms around the peak. 

for sweep_i = 1:size(tempMat,1)
    % estimate peak size post-stimulation period
    clearvars *pkBaseidx     
    % baseline current before peak
    if adjInd_E(sweep_i) < adjInd_I(sweep_i)   % if E precedes I
        pkBaseidx=[adjInd_E(sweep_i)-130,adjInd_E(sweep_i)-30]; % 13ms-3ms prior to peak detection

        if min(pkBaseidx)<=350 && min(pkBaseidx)>=300
            Ebase_beforepk(sweep_i) = mean(tempMat(sweep_i,(pkBaseidx(1)-100):(pkBaseidx(2)-100)),2); 
            Ibase_beforepk(sweep_i) = mean(tempMat(sweep_i,(pkBaseidx(1)-100):(pkBaseidx(2)-100)),2); 
        else
            Ebase_beforepk(sweep_i) = mean(tempMat(sweep_i,pkBaseidx(1):pkBaseidx(2)),2);  
            Ibase_beforepk(sweep_i) = mean(tempMat(sweep_i,pkBaseidx(1):pkBaseidx(2)),2);  
        end

    else % if I precedes E
        IpkBaseidx=[adjInd_I(sweep_i)-130,adjInd_I(sweep_i)-30]; % 13ms-3ms prior to peak detection
        EpkBaseidx=[adjInd_E(sweep_i)-130,adjInd_E(sweep_i)-30];
        if min(IpkBaseidx)<=350 && min(IpkBaseidx)>=300
            Ibase_beforepk(sweep_i) = mean(tempMat(sweep_i,(IpkBaseidx(1)-100):(IpkBaseidx(2)-100)),2); 
        else
            Ibase_beforepk(sweep_i) = mean(tempMat(sweep_i,IpkBaseidx(1):IpkBaseidx(2)),2); 
        end

        if min(EpkBaseidx)<=350 && min(EpkBaseidx)>=300
            Ebase_beforepk(sweep_i) = mean(tempMat(sweep_i,(EpkBaseidx(1)-100):(EpkBaseidx(2)-100)),2); 
        else
            Ebase_beforepk(sweep_i) = mean(tempMat(sweep_i,EpkBaseidx(1):EpkBaseidx(2)),2); 
        end
    end

    Eavg(sweep_i) = mean(tempMat(sweep_i,(adjInd_E(sweep_i)-5):(adjInd_E(sweep_i)+5)),2)-...
        Ebase_beforepk(sweep_i);
    Iavg(sweep_i) = mean(tempMat(sweep_i,(adjInd_I(sweep_i)-5):(adjInd_I(sweep_i)+5)),2)-...
        Ibase_beforepk(sweep_i);

    % do the same for pre-stimulation baseline period
    % first 10ms used for baseline estimate
    base_beforepk_null(sweep_i) = mean(tempMat_base(sweep_i,(E_nullLoc(sweep_i)-130+130):(E_nullLoc(sweep_i)-30+130)),2);  % remember that tempMat_base has additional 130 points in the beginning
    E_nullavg(sweep_i) = mean(tempMat_base(sweep_i,(E_nullLoc(sweep_i)-5+130):(E_nullLoc(sweep_i)+5+130)),2)-...
        base_beforepk_null(sweep_i);
    I_nullavg(sweep_i) = mean(tempMat_base(sweep_i,(I_nullLoc(sweep_i)-5+130):(I_nullLoc(sweep_i)+5+130)),2)-...
        base_beforepk_null(sweep_i);

end
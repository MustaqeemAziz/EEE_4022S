AlignedMatrix = RangeAlignmentStartEnd('Tigress_INB_1_P728_B1_STC88_HRR','velocity estimation',1042,1121,3);
ISARimage = AutofocusMultiple(AlignedMatrix, 'phase difference');
function AlignedMatrix = RangeAlignmentStartEnd(Filename, Technique, StartProfile, EndProfile, OrderOfPolyFit)

HRR_Profile = RadarDataset(StartProfile:EndProfile,:);

figure;
imagesc(20*log10(abs(HRR_Profile)));
title("Unaligned Data");
colormap('jet');
colorbar;

NumBursts = size(HRR_Profile,1);    % Number of bursts in measured data
NumCells = size(HRR_Profile,2);     % Number of range bins in measured data
HRR_reference = HRR_Profile(1,:);

%

Technique to choose from

if strcmpi(Technique,'dominant scatterer')
    
    n_vector =  fftshift(-NumCells/2:1:(NumCells/2-1));
    CorrelationMatrix = zeros(NumBursts, 2*NumCells-1);
    AlignedMatrixlocal = zeros(NumBursts, NumCells);
    CorrelationMatrix(1,:) = xcorr(abs(HRR_reference), abs(HRR_Profile(1,:)));
    AlignedMatrixlocal(1,:) = abs(HRR_Profile(1,:));
    FilterWeight = zeros(1,NumBursts);
    FilterWeight(1) = 0;
    
    for a = 2:1:NumBursts
        CorrelationMatrix(a,:) = xcorr(abs(HRR_reference), abs(HRR_Profile(a,:)));
        [max_numlocal, max_indlocal] = max(CorrelationMatrix(a,:).');
        BinShiftsReqlocal = max_indlocal - NumCells;
        AlignedMatrixlocal(a,:) = abs(ifft(exp(-1i*2*pi*BinShiftsReqlocal*n_vector/NumCells).*(fft(HRR_Profile(a,:),[],2)),[], 2));
        % Use for mean averaging
        HRR_reference = mean((AlignedMatrixlocal(1:1:a,:)));
        % Use for exponential averaging
        FilterWeightInitial = FilterWeight(a-1);
        FilterWeight(a) = 1/a + (1-1/a)*FilterWeightInitial;
        HRR_reference = mean(abs(AlignedMatrixlocal(1:a,:)).*(FilterWeight(:,1:a).'));
        
    end
    
    [max_num, max_ind] = max(CorrelationMatrix.');
    BinShiftsReq = max_ind - NumCells;
    
    NegIdx = find(BinShiftsReq < 0);
    BinShiftsReq(NegIdx) = BinShiftsReq(NegIdx) + NumCells;
    
    Burst_num = (0:1:NumBursts-1);
    Smooth_coeff = polyfit(Burst_num, BinShiftsReq, OrderOfPolyFit);
    Smooth_curve = polyval(Smooth_coeff, Burst_num);
    
    figure;
    plot(Burst_num, BinShiftsReq, '-*', Burst_num, Smooth_curve, '--')
    title('Required Bin Shifts','fontsize',20);
    xlim([0 NumBursts]);
    xlabel('Range Profile Number','fontsize',20);
    ylabel('Bin Shift','fontsize',20);
    fig = gcf;
    fig.Color = [1 1 1];
    set(gca,'fontsize',20)
    
    PhaseShiftRange_New = zeros(NumBursts,NumCells);
    
    for b = 1:1:NumBursts
        PhaseShiftRange_New(b,:) = exp(-1i*2*pi*Smooth_curve(b)*n_vector/NumCells);
    end
    
    AlignedMatrix{1} = ifft(PhaseShiftRange_New.*(fft(HRR_Profile,[],2)),[], 2);
    
elseif strcmpi(Technique,'envelope correlation')
    
    BinShifts = -0.5:0.0025:0.5;                         % Bin shifts to consider in the autocorrelation
    BinShifts_Initial = BinShifts;
    BinShiftFound = zeros(1,NumBursts);
    a = length(BinShifts);
    n_vector =  fftshift(-NumCells/2:1:(NumCells/2-1));
    AlignedMatrixlocal = zeros(NumBursts,NumCells);
    AlignedMatrixlocal(1,:) = HRR_reference;
    
    FilterWeight = zeros(1,NumBursts);
    FilterWeight(1) = 0;
    
    for b = 2:1:NumBursts
        ProfileToAlign = HRR_Profile(b,:);
        
        ShiftMatrix = zeros(a,NumCells);
        
        for c = 1:1:a
            ShiftMatrix(c,:) = exp(-1i*2*pi*BinShifts(c)*n_vector/NumCells);
        end
        
        ShiftedMthProfile = ifft(ShiftMatrix.*(fft(ProfileToAlign,[],2)),[], 2);
        
        figure;
        imagesc(20*log10(abs(ProfileToAlign)));
        title("ProfileToAlign");
        colormap('jet');
        colorbar;
        figure;
        imagesc(20*log10(abs(ShiftedMthProfile)));
        title("ShiftedMthProfile");
        colormap('jet');
        colorbar;
        
        Reference_ProfileMatrix = repmat(HRR_reference,a,1);
        CorrelationValuesMatrix = abs(Reference_ProfileMatrix).*abs(ShiftedMthProfile);
        CorrelationVector = sum(CorrelationValuesMatrix,2);
        
        [max_num, max_ind] = max(CorrelationVector.');
        BinShiftFound(b) = BinShifts(max_ind);
        
        AlignedProfile = ShiftedMthProfile(max_ind,:);
        AlignedMatrixlocal(b,:) = AlignedProfile;
        
        BinShifts = BinShifts_Initial + BinShiftFound(b);
        % Use for mean averaging
        HRR_reference = mean(abs(AlignedMatrixlocal(1:b,:)));
        % Use for exponential averaging
        FilterWeightInitial = FilterWeight(b-1);
        FilterWeight(b) = 1/b + (1-1/b)*FilterWeightInitial;
        HRR_reference = mean(abs(AlignedMatrixlocal(1:b,:)).*(FilterWeight(:,1:b).'));
    end
    
    Burst_num = (1:1:NumBursts);
    Smooth_coeff = polyfit(Burst_num, BinShiftFound, OrderOfPolyFit);
    Smooth_curve = polyval(Smooth_coeff, Burst_num);
    
    figure;
    plot(Burst_num, BinShiftFound, '-*', Burst_num, Smooth_curve, '--')
    title('Required Bin Shifts','fontsize',20);
    xlim([0 NumBursts]);
    xlabel('Range Profile Number','fontsize',20);
    ylabel('Bin Shift','fontsize',20);
    fig = gcf;
    fig.Color = [1 1 1];
    set(gca,'fontsize',20)
    
    PhaseShiftRange_New = zeros(NumBursts,NumCells);
    
    for d = 1:1:NumBursts
        PhaseShiftRange_New(d,:) = exp(-1i*2*pi*(Smooth_curve(d)+10)*n_vector/NumCells);
    end
    
    AlignedMatrix{1} = ifft(PhaseShiftRange_New.*(fft(HRR_Profile,[],2)),[], 2);
    
    
elseif strcmpi(Technique,'Velocity estimation')
    
    n_vector =  fftshift(-NumCells/2:1:(NumCells/2-1)); % look at website to understand the order
    
    a = 1;
    b_initial = 20;
    b = b_initial;
    c = ceil(NumBursts/b);
    MultipleMatrix = cell(c,1);
    
    for d = 1:1:c
        MultipleMatrix{d} = abs(HRR_Profile(a:1:b,:));
        a = a + b_initial;
        b = b + b_initial;
        if b > NumBursts
            b = NumBursts;
        end
    end
    
    AccumulatedProfileMatrix = zeros(c,NumCells);
    AccumulatedProfiles = sum(MultipleMatrix{1});
    LogAccumulatedProfiles = log(AccumulatedProfiles);
    AccumulatedProfileMatrix(1,:) = LogAccumulatedProfiles - mean(LogAccumulatedProfiles);
    Reference_AccumulatedProfile = AccumulatedProfileMatrix(1,:);
    
    CorrelationMatrix_AccumulatedProfile = zeros(c, 2*NumCells-1);
    CorrelationMatrix_AccumulatedProfile(1,:) = xcorr(Reference_AccumulatedProfile, AccumulatedProfileMatrix(1,:));
    
    
    for e = 2:1:length(MultipleMatrix)
        AccumulatedProfiles = sum(MultipleMatrix{e});
        LogAccumulatedProfiles = log(AccumulatedProfiles);
        AccumulatedProfileMatrix(e,:) = LogAccumulatedProfiles - mean(LogAccumulatedProfiles);
        CorrelationMatrix_AccumulatedProfile(e,:) = xcorr(Reference_AccumulatedProfile, AccumulatedProfileMatrix(e,:));
    end
    
    [max_num, max_ind] = max(CorrelationMatrix_AccumulatedProfile.');
    Velocity_Offset = max_ind - NumCells;
    
    NegIdx = find(Velocity_Offset < 0);
    Velocity_Offset(NegIdx) = Velocity_Offset(NegIdx) + NumCells;
    
    Velocity_xaxis = (0:1:c-1);
    Velocity_coeff = polyfit(Velocity_xaxis, Velocity_Offset, 1);
    Velocity_curve = polyval(Velocity_coeff, Velocity_xaxis);
    
    figure;
    plot(Velocity_xaxis, Velocity_Offset, '-*', Velocity_xaxis, Velocity_curve, ':')
    title('Velocity_Offset');
    
    VelocityaxisMultiplier = (NumBursts-1)/(length(Velocity_Offset)-1);
    BinShiftsReqaxis = (VelocityaxisMultiplier+1:VelocityaxisMultiplier:NumBursts);
    
    BinShiftsReq = zeros(1,NumBursts);
    f = 1;
    g = 1;
    BinShiftsReq(f:1:g) = Velocity_Offset(1);
    
    for h = 2:1:c
        f = g + 1;
        g = round(BinShiftsReqaxis(h-1));
        if g > NumBursts
            g = NumBursts;
        end
        BinShiftsReq(f:1:g) = Velocity_Offset(h);
    end
    
    Burst_num = (0:1:NumBursts-1);
    Smooth_coeff = polyfit(Burst_num, BinShiftsReq, 1);
    Smooth_curve = polyval(Smooth_coeff, Burst_num);
    
    figure;
    plot(Burst_num, BinShiftsReq, '-*', Burst_num, Smooth_curve, '--')
    title('Required Bin Shifts','fontsize',20);
    xlim([0 NumBursts]);
    xlabel('Range Profile Number','fontsize',20);
    ylabel('Bin Shift','fontsize',20);
    fig = gcf;
    fig.Color = [1 1 1];
    set(gca,'fontsize',20)
    
    PhaseShiftRange_New = zeros(NumBursts,NumCells);
    
    for g = 1:1:NumBursts
        PhaseShiftRange_New(g,:) = exp(-1i*2*pi*(Smooth_curve(g))*n_vector/NumCells);
    end
    
    
    AlignedMatrix{1} = ifft(PhaseShiftRange_New.*(fft(HRR_Profile,[],2)),[], 2);
    
    
end

AlignedMatrix1 = AlignedMatrix{1};
AlignedMatrix1 = (normr(abs(AlignedMatrix1))).';
BinEntropy = zeros(1,NumCells);
for z = 1:1:NumCells
    
    BinEntropy(z) = mean(abs(AlignedMatrix1(:,z)));
    
end

Entropy = -1*mean(BinEntropy.*log(BinEntropy));


figure;
imagesc(20*log10(abs(AlignedMatrix{1})));
title('HRR Profile with Alignment','fontsize',20);
colormap('jet');
colorbar;
xlabel('Range Bin Number','fontsize',20);
ylabel('Range Profile Number','fontsize',20);
fig = gcf;
fig.Color = [1 1 1];
dim = [0.15 0.6 0.3 0.3];
str = {'Entropy: ',Entropy};
annotation('textbox',dim,'BackgroundColor',[0.94 0.97 1],'String',str,'Fontsize',28,'FitBoxToText','on');
set(gca,'fontsize',28)

end
function ISARimage = AutofocusMultiple(AlignedMatrix, Technique)

C = 3e8;
%% Use for creating multiple CPI's
AlignedMatrix = cell2mat(AlignedMatrix);
NumBursts = size(AlignedMatrix,1);    % Number of bursts in measured data

y = 1;
x_initial = 512;
x = x_initial;
w = 2*ceil(NumBursts/x)-1;
MultipleMatrixCell = cell(w,1);

for v = 1:1:w
    MultipleMatrixCell{v} = num2cell(AlignedMatrix(y:1:x,:));
    y = y + x_initial/2;
    x = x + x_initial/2;
    if x > NumBursts
        x = NumBursts;
    end
end

CPI_Length = length(MultipleMatrixCell);
if strcmpi(Technique,'phase difference')
    
    for CPI_n = 1:1:CPI_Length
        
        AlignedMatrixLocal = AlignedMatrix(MultipleMatrixCell{CPI_n});
        NumCellslocal = size(AlignedMatrixLocal,2);
        NumBurstslocal = size(AlignedMatrixLocal,1);
        
        RangeAxis = (0:1:(NumCellslocal-1))*C/(2*NumCellslocal*4e6);
        DopplerFreq_Hz = ((-NumBurstslocal/2):1:(NumBurstslocal/2-1))*154/NumBurstslocal;
        
        AlignedMatrixLocalAmplitude = abs(AlignedMatrixLocal);
        RangeBinVariance = var(AlignedMatrixLocalAmplitude,[],1);
        PowerMatrix = AlignedMatrixLocalAmplitude.^2;
        RangeBinPower = sum(PowerMatrix,1);
        ThresholdPower = (1/NumCellslocal)*sum(RangeBinPower);
        DSV = max(RangeBinVariance);
        
        b = length(RangeBinVariance);
        
        for c = 1:1:b
            if RangeBinPower(c) > ThresholdPower
                if RangeBinVariance(c) < DSV
                    DSV = RangeBinVariance(c);
                    DS = c;
                end
            end
        end
        
        PhaseShift = phase(AlignedMatrixLocal(1,DS));
        
        PhaseVector = zeros(size(AlignedMatrixLocal,1),1);
        
        for d = 1:1:size(AlignedMatrixLocal,1)
            PhaseVector(d) = phase(AlignedMatrixLocal(d,DS)) - PhaseShift;
        end
        
        PhaseShiftMultiplier = exp(-1i*PhaseVector);
        AutoFocusMatrix = AlignedMatrixLocal.*PhaseShiftMultiplier;
        
        WindowVector = blackman(NumBurstslocal);
        WindowMatrix = repmat(WindowVector,1, NumCellslocal);
        
        ISARimage = fftshift(fft(AutoFocusMatrix.*WindowMatrix, [], 1),1); % include windowing to reduce sidelobes in Doppler dimension
        
        ISARimage1 = normc((abs(ISARimage)));
        
        BinEntropyISAR = zeros(1,NumCellslocal);
        for z = 1:1:NumCellslocal
            BinEntropyISAR(z) = mean((ISARimage1(z,:)));
        end
        figure;
        plot(BinEntropyISAR);
        EntropyISAR = -1*mean(BinEntropyISAR.*log(BinEntropyISAR));
        figure;
        imagesc(RangeAxis,DopplerFreq_Hz,20*log10(abs(ISARimage)));
        title('ISAR image','fontsize',20);
        colormap('jet');
        colorbar;
        xlabel('Down-range (m)','fontsize',20);
        ylabel('Doppler (Hz)','fontsize',20);
        fig = gcf;
        fig.Color = [1 1 1];
        dim = [0.15 0.6 0.3 0.3];
        str = {'Entropy: ',EntropyISAR};
        annotation('textbox',dim,'BackgroundColor',[0.94 0.97 1],'String',str,'Fontsize',20,'FitBoxToText','on');
        set(gca,'fontsize',20)
        axis xy;
        
    end
    
    
elseif strcmpi(Technique,'phase conjugate')
    
    for CPI_n = 1:1:CPI_Length
        
        AlignedMatrixLocal = AlignedMatrix(MultipleMatrixCell{CPI_n}));
        AlignedMatrixLocal=AlignedMatrix{CPI_n};
        AlignedMatrixLocal = AlignedMatrixLocal(185:587,:);
        NumCellslocal = size(AlignedMatrixLocal,2);
        NumBurstslocal = size(AlignedMatrixLocal,1);
        
        RangeAxis = (0:1:(NumCellslocal-1))*C/(2*NumCellslocal*4e6);
        DopplerFreq_Hz = ((-NumBurstslocal/2):1:(NumBurstslocal/2-1))*154/NumBurstslocal;
        
        AlignedMatrixLocalAmplitude = abs(AlignedMatrixLocal);
        RangeBinVariance = var(AlignedMatrixLocalAmplitude,[],1);
        PowerMatrix = AlignedMatrixLocalAmplitude.^2;
        RangeBinPower = sum(PowerMatrix,1);
        ThresholdPower = (1/NumCellslocal)*sum(PowerMatrix(:));
        DSV = max(RangeBinVariance);
        
        b = length(RangeBinVariance);
        
        for d = 1:1:b
            if RangeBinPower(d) > ThresholdPower
                if RangeBinVariance(d) < DSV
                    DSV = RangeBinVariance(d);
                    DS = d;
                end
            end
        end
        
        DSRangeBin = conj(AlignedMatrixLocal(:,DS));
        DSMatrix = repmat(DSRangeBin,1,NumCellslocal);
        AutoFocusMatrix = AlignedMatrixLocal.*DSMatrix;
        
        WindowVector = blackman(NumBurstslocal);
        WindowMatrix = repmat(WindowVector,1, NumCellslocal);
        
        ISARimage = fftshift(fft(AutoFocusMatrix.*WindowMatrix, [], 1),1); % include windowing to reduce sidelobes in Doppler dimension
        
        
        
        ISARimage1 = normc((abs(ISARimage)));
        
        BinEntropyISAR = zeros(1,NumCellslocal);
        for z = 1:1:NumCellslocal
            BinEntropyISAR(z) = mean((ISARimage1(z,:)));
        end
        figure;
        plot(BinEntropyISAR);
        EntropyISAR = -1*mean(BinEntropyISAR.*log(BinEntropyISAR));
        figure;
        imagesc(RangeAxis,DopplerFreq_Hz,20*log10(abs(ISARimage)));
        title('ISAR image','fontsize',20);
        colormap('jet');
        colorbar;
        xlabel('Down-range (m)','fontsize',20);
        ylabel('Doppler (Hz)','fontsize',20);
        fig = gcf;
        fig.Color = [1 1 1];
        dim = [0.15 0.6 0.3 0.3];
        str = {'Entropy: ',EntropyISAR};
        annotation('textbox',dim,'BackgroundColor',[0.94 0.97 1],'String',str,'Fontsize',20,'FitBoxToText','on');
        set(gca,'fontsize',20)
        axis xy;
        
    end
    
end

end







cd(FOV1.defaultDir);
syncObj = bSyncIm(FOV1);

%% Select Slice and run Patches
for nSlice = 1:4
    acqBlocks = [1 syncObj.sliceFrames(1,nSlice)];
    for blockNum = 2:size(syncObj.sliceFrames,1)
        acqBlocks(blockNum,:) = ...
            [1+syncObj.sliceFrames(blockNum-1,nSlice), syncObj.sliceFrames(blockNum,nSlice)];
    end
    mapFile = sprintf('%sFOV1_memmap_slice%0.2d.mat',[cd filesep],nSlice);
    data = matfile(mapFile);
    saveFile = sprintf('%spatchResults_slice%0.2d_v0820',FOV1.defaultDir,nSlice);
    runSNCpatchNMF(data,saveFile,acqBlocks)
end
Acquisition2P([],@SClk);
cd(FOV1.defaultDir)
motionCorrect(FOV1)
writeDir = 'E:\memmaps\19-160709';
for i=1:4
    indexMovie(FOV1,i,1,writeDir);
end
FOV1.syncInfo = bSyncIm(FOV1);
fprintf('%d frames from Acq Obj \n',sum(FOV1.correctedMovies.slice(end).channel.size(:,3))),
FOV1.save,

%% Select Slice and run Patches
for nSlice = 1:4
    extractSourcesNMF(FOV1,nSlice);
end
FOV1.save

%% Deconvolve and save clustered sources
cd(FOV1.defaultDir),
[dF,dF_denoised,dF_deconv,...
    traceBs,traceGs,traceSNs,traceSnScales,A,b,f] = deconv_NMF(FOV1);
save('cellData_0921','dF','dF_denoised','dF_deconv',...
    'traceBs','traceGs','traceSNs','traceSnScales','A','b','f'),

function [data,timesteps] =merge_tiles(param,outdir,model,scenario,startyear,endyear,uppery,lowery,leftx,rightx,tilesize,res)

rows_tiles=tilesize/res;
rows_global=(uppery-lowery)/res;
cols=abs(leftx-rightx)/res;

for upperyy=uppery:-tilesize:lowery+tilesize
  loweryy=upperyy-tilesize;
  yrow=floor((uppery-upperyy)/res)+1;
  filename=[outdir,filesep,model,'_',num2str(upperyy),'_',num2str(loweryy),'_',num2str(leftx),'_',num2str(rightx),filesep,param,'_',model,'_',scenario,'_',num2str(startyear),'-',num2str(endyear)];
  [data_tile timesteps]=readras(filename,1,1000);
  data(yrow:yrow+rows_tiles-1,1:cols,:)=data_tile(1:rows_tiles,1:cols,:);
end

end

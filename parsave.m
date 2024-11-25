function parsave(fname,val,output,timesteps,data_cm_var)
  save(fname, 'val','output','timesteps','data_cm_var','-v7.3')
end
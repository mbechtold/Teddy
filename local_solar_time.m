%calculates local solar time offset to UTC for each pixel center
function lst_offset=local_solar_time(res)
s=15/res;
a=1;
for n=-12:12
  lst_offset(a:a+s-1)=repmat(n,[1,s]);
  a=a+s;
end
lst_offset(1:s/2)=[];
lst_offset(end-s/2+1:end)=[];
return;
end
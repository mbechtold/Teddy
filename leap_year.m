function b=leap_year(year)
if (rem(year,400)==0 | rem(year,4)==0 & rem(year,100)~=0)
  b=1;
else
  b=0;
end
b=logical(b);
return
end
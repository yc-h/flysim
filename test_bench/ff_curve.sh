#!/bin/bash
rm FI.dat
declare -i N
N=0
end=300
prestrength=0
initialstrength=0
strength=0
prestrength=$initialstrength
#-------------------------------------
while [ "  $N " -lt " $end " ]
do
i=`expr $N + 1 `
echo $N,$strength
./flysim.out -s rough > log
rm log

strength=`gawk -v pres=$prestrength 'BEGIN{stre=0;}
{ 
	stre=pres+30;
       }
END{print stre;}' network.pro`


sed -i 6s/$prestrength/$strength/g network.pro
prestrength=$strength


  mean=`gawk -v N=$N 'BEGIN{n=0;f=0;y=0;}
  {
    if($1>2)
    {
      n++;
      f=f+$2;
      y=y+$3;
    }
    else { }
  }

END{print (N+1)*30,f/n,y/n;}' Freq_0.dat>>FI.dat`


N=`expr $N + 1 `
done
echo "job complete"


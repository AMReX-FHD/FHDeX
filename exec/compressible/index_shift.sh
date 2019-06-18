#!/bin/bash

nvars=4

declare -a varnames=("cu" "fluxx" "fluxy" "fluxz")
declare -a indxnames=("i" "j" "k")
declare -A nodality
declare -a loopnode
declare -i ndims=3

declare -i loopcnt=0
declare -i inloop=0

nodekey="-1"

# for ((i=1;i<=num_rows;i++)) do
#     for ((j=1;j<=num_columns;j++)) do
#         nodality[$i,$j]=$RANDOM
#     done
# done

for ((d=0;d<3;d++)) do
loopnode[$d+1]=0
done

nodality[1,1]=0
nodality[1,2]=0
nodality[1,3]=0

nodality[2,1]=1
nodality[2,2]=0
nodality[2,3]=0

nodality[3,1]=0
nodality[3,2]=1
nodality[3,3]=0

nodality[4,1]=0
nodality[4,2]=0
nodality[4,3]=1

## now loop through the above array
for ((i=1;i<=nvars;i++)) do
    for ((j=1;j<=3;j++)) do
        echo ${nodality[$i,$j]}
    done
done

# Note: must be on same level as spatial loop
while IFS= read -r line; do
    echo "Text read from file: $line"

    # Detect beginning of loops
    if echo "$line" | grep -q "do .* ="; then
	echo "loop detected";
	loopcnt+=1
	# echo "$loopcnt";

	if [ "$loopcnt" -eq "$ndims" ] && [ "$inloop" -eq "0" ] 
	then
	    echo "entering loop";
	    inloop=1
	fi
	
	for ((d=0;d<3;d++)) do
	# echo "Hack: ${indxnames[$d]}";

	if echo "$line" | grep -q "${indxnames[$d]}.*="; then
	    echo "beginning of ${indxnames[$d]} loop detected";

	    if echo "$line" | grep -q -- "-1"; then
		echo "nodal in ${indxnames[$d]} direction";

	    	loopnode[$d+1]=1
	    fi
	fi
	done
    fi
    
    # Detect end of loops
    if echo "$line" | grep -q "end.*do"; then
	echo "endloop detected";
	loopcnt+=-1
	# echo "$loopcnt";
	if [ "$loopcnt" -eq "0" ]
	then
	    echo "end of spatial loop detected";
            inloop=0

	    for ((d=0;d<3;d++)) do
	    loopnode[$d+1]=0
	    done
	fi
    fi
    
    # Replace string occurances
    if [ "$inloop" -eq "1" ] 
    then
    	echo "inside loop";
	for ((d=0;d<3;d++)) do
	echo "Nodality in ${indxnames[$d]} = ${loopnode[$d+1]}";
	done

    	for ((n=0;n<nvars;n++)) do
    	if echo "$line" | grep -q "${varnames[$n]}"; then
    	    echo "fixing nodality for ${varnames[$n]}";
	    
    	    for ((d=0;d<3;d++)) do
	    echo "comparing loop nodality ${loopnode[$((d+1))]} to mf nodality ${nodality[$((n+1)),$((d+1))]}";

    	    if [ "${loopnode[$((d+1))]}" -eq "1" ] && [ "${nodality[$((n+1)),$((d+1))]}" -eq "0" ]
	    then
		echo "testing cases...";

    	    	case $d in

    	    	    0)
    	    		echo "cc in x"
    	    		line=$(echo "$line" | sed "s/${varnames[$n]}(i+1,\(.*\))/${varnames[$n]}(i,\1)/g")
    	    		line=$(echo "$line" | sed "s/${varnames[$n]}(i,\(.*\))/${varnames[$n]}(i-1,\1)/g")
    	    		;;

    	    	    1)
    	    		echo "cc in y"
    	    		;;

    	    	    2)
    	    		echo "cc in z"
    	    		;;

    	    	esac
		
		echo "Replaced line:       $line"

    	    fi
    	    done
    	fi
    	done
    fi

done < "$1"

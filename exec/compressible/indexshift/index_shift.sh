#!/bin/bash

outputFilename="output.f90"

nvars=23

# foobar in place of zeta - redundant with "eta" (both share nodality)

declare -a varnames=("prim" "con" "fluxx" "fluxy" "fluxz" "eta" "foobar" "kappa" "chi" "Dij"
    "visccorn" "cornux" "cornuy" "cornuz" "cornvx" "cornvy" "cornvz" "cornwx" "cornwy" "cornwz"
    "ranfluxx" "ranfluxy" "ranfluxz" "rancorn")

declare -a indxnames=("i" "j" "k")
declare -A nodality
declare -a loopnode
declare -i ndims=3

declare -i loopcnt=0
declare -i inloop=0

d05="....."
d06="......"
d07="......."
d08="........"
d09="........."
d10=".........."
d11="..........."
d12="............"

regx="\(\|..\|....\)"
# regxnd="\(\|.\|...\)"
regxnd="\(\|.\|...\|....\|$d05\|$d06\|$d07\|$d08\|$d09\|$d10\|$d11\|$d12\)"

# for ((i=1;i<=num_rows;i++)) do
#     for ((j=1;j<=num_columns;j++)) do
#         nodality[$i,$j]=$RANDOM
#     done
# done

for ((d=0;d<3;d++)) do
    loopnode[$d]=0
done

nodality[0,0]=0; nodality[0,1]=0; nodality[0,2]=0
nodality[1,0]=0; nodality[1,1]=0; nodality[1,2]=0

nodality[2,0]=1; nodality[2,1]=0; nodality[2,2]=0
nodality[3,0]=0; nodality[3,1]=1; nodality[3,2]=0
nodality[4,0]=0; nodality[4,1]=0; nodality[4,2]=1

nodality[5,0]=0; nodality[5,1]=0; nodality[5,2]=0
nodality[6,0]=0; nodality[6,1]=0; nodality[6,2]=0
nodality[7,0]=0; nodality[7,1]=0; nodality[7,2]=0
nodality[8,0]=0; nodality[8,1]=0; nodality[8,2]=0
nodality[9,0]=0; nodality[9,1]=0; nodality[9,2]=0

nodality[10,0]=1; nodality[10,1]=1; nodality[10,2]=1
nodality[11,0]=1; nodality[11,1]=1; nodality[11,2]=1
nodality[12,0]=1; nodality[12,1]=1; nodality[12,2]=1
nodality[13,0]=1; nodality[13,1]=1; nodality[13,2]=1
nodality[14,0]=1; nodality[14,1]=1; nodality[14,2]=1
nodality[15,0]=1; nodality[15,1]=1; nodality[15,2]=1
nodality[16,0]=1; nodality[16,1]=1; nodality[16,2]=1
nodality[17,0]=1; nodality[17,1]=1; nodality[17,2]=1
nodality[18,0]=1; nodality[18,1]=1; nodality[18,2]=1
nodality[19,0]=1; nodality[19,1]=1; nodality[19,2]=1

nodality[20,0]=1; nodality[20,1]=0; nodality[20,2]=0
nodality[21,0]=0; nodality[21,1]=1; nodality[21,2]=0
nodality[22,0]=0; nodality[22,1]=0; nodality[22,2]=1

nodality[23,0]=1; nodality[23,1]=1; nodality[23,2]=1


## now loop through the above array
for ((i=0;i<nvars;i++)) do
    for ((j=0;j<3;j++)) do
        echo ${nodality[$i,$j]}
    done
done

rm ${outputFilename}

# Note: must be on same level as spatial loop
# while IFS= read -r line; do
while IFS= read -r line || [ -n "$line" ]; do
    echo "Text read from file: $line"

    # Detect beginning of loops
    if echo "$line" | grep -q "do .*="; then
        echo "   loop detected";
        loopcnt+=1
        echo "   loop count = $loopcnt";

        if [ "$loopcnt" -eq "$ndims" ] && [ "$inloop" -eq "0" ]
        then
            echo "   entering loop";
            inloop=1
        fi

        for ((d=0;d<3;d++)) do
            # echo "Hack: ${indxnames[$d]}";

            if echo "$line" | grep -q "${indxnames[$d]}.*="; then
                echo "   beginning of ${indxnames[$d]} loop detected";

                if echo "$line" | grep -q -- "-1"; then
                    echo "   nodal in ${indxnames[$d]} direction";

                    loopnode[$d]=1
                fi
            fi
        done
    fi

    # Detect end of loops
    if echo "$line" | grep -q "end.*do"; then
        echo "   endloop detected";
        loopcnt+=-1
        echo "   loop count = $loopcnt";

        if [ "$loopcnt" -eq "0" ]
        then
            echo "  end of spatial loop detected";
            inloop=0

            for ((d=0;d<3;d++)) do
                loopnode[$d]=0
            done
        fi
    fi

    # Replace string occurances
    if [ "$inloop" -eq "1" ]
    then
        echo "   inside loop";
        # for ((d=0;d<3;d++)) do
        # echo "   nodality in ${indxnames[$d]} = ${loopnode[$d]}";
        # done

        linetemp="$line"

        for ((n=0;n<nvars;n++)) do
            if echo "$line" | grep -q "${varnames[$n]}("; then
                echo "   checking nodality for ${varnames[$n]}";

                for ((d=0;d<3;d++)) do
                    echo "   nodality: loop ${loopnode[$((d))]}, mf ${nodality[$((n)),$((d))]}";

                    if [ "${loopnode[$((d))]}" -eq "1" ] && [ "${nodality[$((n)),$((d))]}" -eq "0" ]
                    then
                        echo "   loop nodal, mf cell centered...";

                        # Note: ordering of sed statements matters
                        case $d in

                            0)
                                echo "   cc in x"
                                linetemp=$(echo "$linetemp" | sed -e "s/${varnames[$n]}(i-1,${regx}${regx}${regxnd})/${varnames[$n]}(i-2,\1\2\3)/g")
                                linetemp=$(echo "$linetemp" | sed -e "s/${varnames[$n]}(i,${regx}${regx}${regxnd})/${varnames[$n]}(i-1,\1\2\3)/g")
                                linetemp=$(echo "$linetemp" | sed -e "s/${varnames[$n]}(i+1,${regx}${regx}${regxnd})/${varnames[$n]}(i,\1\2\3)/g")
                                linetemp=$(echo "$linetemp" | sed -e "s/${varnames[$n]}(i+2,${regx}${regx}${regxnd})/${varnames[$n]}(i+1,\1\2\3)/g")
                                ;;

                            1)
                                echo "   cc in y"
                                linetemp=$(echo "$linetemp" | sed -e "s/${varnames[$n]}(${regx}j-1,${regx}${regxnd})/${varnames[$n]}(\1j-2,\2\3)/g")
                                linetemp=$(echo "$linetemp" | sed -e "s/${varnames[$n]}(${regx}j,${regx}${regxnd})/${varnames[$n]}(\1j-1,\2\3)/g")
                                linetemp=$(echo "$linetemp" | sed -e "s/${varnames[$n]}(${regx}j+1,${regx}${regxnd})/${varnames[$n]}(\1j,\2\3)/g")
                                linetemp=$(echo "$linetemp" | sed -e "s/${varnames[$n]}(${regx}j+2,${regx}${regxnd})/${varnames[$n]}(\1j+1,\2\3)/g")
                                ;;

                            2)
                                echo "   cc in z"

                                # with 4th index
                                linetemp=$(echo "$linetemp" | sed -e "s/${varnames[$n]}(${regx}${regx}k-1,${regxnd})/${varnames[$n]}(\1\2k-2,\3)/g")
                                linetemp=$(echo "$linetemp" | sed -e "s/${varnames[$n]}(${regx}${regx}k,${regxnd})/${varnames[$n]}(\1\2k-1,\3)/g")
                                linetemp=$(echo "$linetemp" | sed -e "s/${varnames[$n]}(${regx}${regx}k+1,${regxnd})/${varnames[$n]}(\1\2k,\3)/g")
                                linetemp=$(echo "$linetemp" | sed -e "s/${varnames[$n]}(${regx}${regx}k+2,${regxnd})/${varnames[$n]}(\1\2k+1,\3)/g")

                                # without 4th index
                                linetemp=$(echo "$linetemp" | sed -e "s/${varnames[$n]}(${regx}${regx}k-1)/${varnames[$n]}(\1\2k-2)/g")
                                linetemp=$(echo "$linetemp" | sed -e "s/${varnames[$n]}(${regx}${regx}k)/${varnames[$n]}(\1\2k-1)/g")
                                linetemp=$(echo "$linetemp" | sed -e "s/${varnames[$n]}(${regx}${regx}k+1)/${varnames[$n]}(\1\2k)/g")
                                linetemp=$(echo "$linetemp" | sed -e "s/${varnames[$n]}(${regx}${regx}k+2)/${varnames[$n]}(\1\2k+1)/g")
                                ;;

                        esac

                        echo "Original line:       $line"
                        line="$linetemp"
                        echo "Modified line:       $line"

                        # echo "Suggested update:    $linetemp"
                        # echo "do you wish to change line?"
                        # select modchoice in "yes" "no"; do
                        #     case $modchoice in
                        #         yes ) line="$linetemp"; break;;
                        #         no ) echo "Line not changed:    $line"; exit;;
                        #     esac
                        # done
                        # # echo "Modified line:       $line"

                    fi

                    if [ "${loopnode[$((d))]}" -eq "0" ] && [ "${nodality[$((n)),$((d))]}" -eq "1" ]
                    then
                        echo "   loop cell centered, mf nodal...";

                        # Note: ordering of sed statements matters
                        case $d in

                            0)
                                echo "   nodal in x"
                                linetemp=$(echo "$linetemp" | sed -e "s/${varnames[$n]}(i+1,${regx}${regx}${regxnd})/${varnames[$n]}(i+2,\1\2\3)/g")
                                linetemp=$(echo "$linetemp" | sed -e "s/${varnames[$n]}(i,${regx}${regx}${regxnd})/${varnames[$n]}(i+1,\1\2\3)/g")
                                linetemp=$(echo "$linetemp" | sed -e "s/${varnames[$n]}(i-1,${regx}${regx}${regxnd})/${varnames[$n]}(i,\1\2\3)/g")
                                linetemp=$(echo "$linetemp" | sed -e "s/${varnames[$n]}(i-2,${regx}${regx}${regxnd})/${varnames[$n]}(i-1,\1\2\3)/g")
                                ;;

                            1)
                                echo "   nodal in y"
                                linetemp=$(echo "$linetemp" | sed -e "s/${varnames[$n]}(${regx}j+1,${regx}${regxnd})/${varnames[$n]}(\1j+2,\2\3)/g")
                                linetemp=$(echo "$linetemp" | sed -e "s/${varnames[$n]}(${regx}j,${regx}${regxnd})/${varnames[$n]}(\1j+1,\2\3)/g")
                                linetemp=$(echo "$linetemp" | sed -e "s/${varnames[$n]}(${regx}j-1,${regx}${regxnd})/${varnames[$n]}(\1j,\2\3)/g")
                                linetemp=$(echo "$linetemp" | sed -e "s/${varnames[$n]}(${regx}j-2,${regx}${regxnd})/${varnames[$n]}(\1j-1,\2\3)/g")
                                ;;

                            2)
                                echo "   nodal in z"

                                # with 4th index
                                linetemp=$(echo "$linetemp" | sed -e "s/${varnames[$n]}(${regx}${regx}k+1,${regxnd})/${varnames[$n]}(\1\2k+2,\3)/g")
                                linetemp=$(echo "$linetemp" | sed -e "s/${varnames[$n]}(${regx}${regx}k,${regxnd})/${varnames[$n]}(\1\2k+1,\3)/g")
                                linetemp=$(echo "$linetemp" | sed -e "s/${varnames[$n]}(${regx}${regx}k-1,${regxnd})/${varnames[$n]}(\1\2k,\3)/g")
                                linetemp=$(echo "$linetemp" | sed -e "s/${varnames[$n]}(${regx}${regx}k-2,${regxnd})/${varnames[$n]}(\1\2k-1,\3)/g")

                                # without 4th index
                                linetemp=$(echo "$linetemp" | sed -e "s/${varnames[$n]}(${regx}${regx}k+1)/${varnames[$n]}(\1\2k+2)/g")
                                linetemp=$(echo "$linetemp" | sed -e "s/${varnames[$n]}(${regx}${regx}k)/${varnames[$n]}(\1\2k+1)/g")
                                linetemp=$(echo "$linetemp" | sed -e "s/${varnames[$n]}(${regx}${regx}k-1)/${varnames[$n]}(\1\2k)/g")
                                linetemp=$(echo "$linetemp" | sed -e "s/${varnames[$n]}(${regx}${regx}k-2)/${varnames[$n]}(\1\2k-1)/g")
                                ;;

                        esac

                        echo "Original line:       $line"
                        line="$linetemp"
                        echo "Modified line:       $line"

                    fi
                done
            fi
        done
    fi

    # Write output to file
    echo "$line" >> ${outputFilename}

done < "$1"

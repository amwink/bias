#!/bin/bash

CD=$PWD;
DD=/scratch2/ibrouwer2/appsms/appsms_import/sourcedata;
DF=${CD}/dicoms_list.txt;
OF=${CD}/dicoms_prot.txt;

# find one dicoms per (non-empty, file-containing) directory
if [[ ! -f $DF ]]; then
    find $DD -type d -exec sh -c 'find "{}" -maxdepth 1 -type f -name \*dcm\* | sort | head -n 3 | tail -n 1' ";" > ${CD}/dicoms_list.txt;
fi

# get all the parameters used in The Word Document
if [[ ! -f $OF ]]; then

    printf "%16s\t%6s\t%-8s\t%-8s\t%-8s\t%-8s\t%s\n" intent.dir est.prot. TE TI TR seq.guess fname > $OF;
        
    while read dcmname; do
    
	intent=$( basename ${dcmname%/*} ); # this directory specifies the intent
	intent=${intent#*_}
	intent=${intent#*_}
	
	dump=$(dcmdump $dcmname);
	
	TE=$(grep EchoTime <<< $dump | awk '{print $3}' | tr "[]" "  " | xargs);
	if [[ "$TE" = "" ]]; then TE=0; fi;
	
	TI=$(grep InversionTime <<< $dump | awk '{print $3}' | tr "[]" "  " | xargs);
	if [[ "$TI" = "" ]]; then TI=0; fi;
	
	TR=$(grep RepetitionTime <<< $dump | awk '{print $3}' | tr "[]" "  " | xargs);
	if [[ "$TR" = "" ]]; then TR=0; fi;
	
	# TR may be in s not ms
	if (( $(echo "$TR < 10" | bc -l) )); then
	    TR=$(echo "$TR * 1000" | bc -l);
	fi
	
	SE=$(grep ScanningSequence <<< $dump | awk '{print $3}' | tr "[]" "  " | xargs);
	if [[ "$SE" = "" ]]; then SE="n/a"; fi;
	
	#assuming: if no "SE" in this string, then GE
	if [[ ${SE//SE/} != ${SE} ]]; then
	    SE="SE";
	else
	    if [[ ${SE//GR/} != ${SE} ]]; then
		SE="GR";
	    else
		SE="NA";
	    fi;
	fi;
	
	if (( $(echo "$TR > 0" | bc -l) )); then
	    
	    if (( $(echo "$TR < 500" | bc -l) )); then
		prot="SWI"
	    else
		if (( $(echo "$TI > 0" | bc -l) )); then
		    if [[ $SE == "SE" ]]; then
			prot="FLAIR"
		    else
			if [[ $SE == "GR" ]]; then
			    prot="MPRAGE"
			else
			    prot="NA"
			fi
		    fi
		else
		    if [[ $SE == "SE" ]]; then
			prot="DTI"
		    else
			if [[ $SE == "GR" ]]; then
			    prot="fMRI"
			else
			    prot="NA"
			fi
		    fi
		fi    
	    fi
	    
	else
	    
	    prot="NA";
	    
	fi
	
	printf "%16s\t%6s\t%-8.2f\t%-8.2f\t%-8.2f\t%-8s\t%s\n" "${intent:0:15}" "$prot" "$TE" "$TI" "$TR" "$SE" $dcmname;

	printf "." >&2
	
    done < $DF >> $OF;
    
fi 

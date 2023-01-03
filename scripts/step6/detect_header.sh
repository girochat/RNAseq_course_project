#!/usr/bin/env bash

# path to the project directory
HOME=/data/courses/rnaseq_course/lncRNAs/Project2/users/grochat

# extract the file ID for the output file
if [ $2 = "gtf" ] 
then
	file_ID=$(basename $1 .gtf)
else
	file_ID=$(basename $1 .bed)
fi

# get all possible headers (1st column) of the bed/gtf file given in argument
awk 'BEGIN{
	i=1
}{
	header=$1
	if (!(header in headers)){
		headers[header]=1
		headers_value[i]=header
		i++
	}
}END{
	print length(headers) 
	for (value in headers_value){
		print headers_value[value]
	} 
}' $1 > $HOME/scripts/${file_ID}_headers.txt

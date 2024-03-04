while read Have_DEG_Group No_DEG_Group
do
	#statements
	#echo "Groups containing differential genes are being processed"
	mkdir -p ./summary/Have_DEG_Group/$Have_DEG_Group
	find -name "*$Have_DEG_Group*" -exec cp {} ./summary/Have_DEG_Group/$Have_DEG_Group/ \;
	#echo "Groups that do not contain differential genes are being processed"
	mkdir -p ./summary/No_DEG_Group/$No_DEG_Group
	find -name "*$No_DEG_Group*" -exec cp {} ./summary/No_DEG_Group/$No_DEG_Group/ \;
	find -name "NA" -exec rm -r {} \; 
done < deg_group.txt

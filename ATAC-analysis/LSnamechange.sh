#a ton of annoying manipulating the file names to get them to look the same in preparation for ChromVAR analysis

cd LSA2
for oldname in *rmdup*; do
	newname=`basename $oldname | sed 's/A_/-LSA2-A_/' | sed 's/B_/-LSA2-B_/' | sed 's/C_/-LSA2-C_/' | sed 's/Q/0/'`
	mv $oldname $newname
done
cd ../LSA3
for oldname in *rmdup*; do
	newname=`basename $oldname | sed 's/LSA3-PS//' | sed 's/LSA3-C/C/' | sed 's/-1_/-LSA3-A_/' | sed 's/-2_/-LSA3-B_/' | sed 's/-3_/-LSA3-C_/'`
	mv $oldname $newname
done
cd ../LSA5
for oldname in *rmdup*; do
	newname=`basename $oldname | sed 's/D//' | sed 's/yc//' | sed 's/-1_/-LSA5-A_/' | sed 's/-2_/-LSA5-B_/' | sed 's/-3_/-LSA5-C_/'`
	mv $oldname $newname
done
cd ..

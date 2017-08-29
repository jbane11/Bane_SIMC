### Shell script for both making all the links and running th monte carlo

if [[ "$#" -eq "1" ]] 
then 
	RUN=$1
fi
if [[ "$#" -eq "0" ]] 
then 
	echo "Please enter run number " 
	read -e RUN
fi

[ -f mcrun_info/input_run$RUN.dat ] && echo "Linking the correct mc info file" || echo "Please create input_run$RUN.dat in /mcrun_info" 
ln -sf mcrun_info/input_run$RUN.dat input.dat

[ -f datarun_info/run$RUN.in ] && echo "Linking the correct data run file" || echo "Please create run$RUN.in in /datarun_info" 
ln -sf datarun_info/run$RUN.in reconmc.in

[ -f rzdat/run$RUN.rzdat ] && echo "Linking the correct rzdat run file" || echo "Please generate the rzdat file needed"
ln -sf rzdat/run$RUN.rzdat run.rzdat


if [ ! -f "rzdat/run$RUN.rzdat" ]
then
	echo "Exiting, somthing is wrong"
	exit
fi



sh reconmc.sh "$RUN"
echo " "
echo "Just finish run # $RUN!!!" 
echo " "
echo " Let's make those ntuples root files and move them to the proper places!!!"
echo " "
h2root mc\ \ "$RUN.rzdat"
echo " "
mv mc\ \ "$RUN.root" "rootfiles/mc$RUN.root"
echo " "
mv mc\ \ "$RUN.rzdat" "mc_rzdat/mc$RUN.rzdat"


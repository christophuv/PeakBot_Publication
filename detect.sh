
#echo Cleaning old models and logs
#rm -r logs/*
#echo

cd /home/users/cbueschl/LCHRMS-GPU/peakbot_example

echo Cleaning old figures
rm -r *.png
echo

echo Detecting with PeakBot
echo

python detectPB.py
echo
echo

cd /home/users/cbueschl/LCHRMS-GPU
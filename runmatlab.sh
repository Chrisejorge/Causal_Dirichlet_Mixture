#! /bin/bash
# Import environmental variable
# . /etc/profile
# . ~/.bash_profile
# Call MATLAB with the appropriate input and output,
# make it immune to hangups and quits using ''nohup'',
# and run it in the background.

nohup matlab -nodisplay -nodesktop -nosplash < $1 > $2 &
echo $! > save_pid.txt &


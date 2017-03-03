# Before running this script, go into the pulser main menu and under menu 6 (Operating Mode) change "PB5 set LOCAL" to "PB5 set REMOTE"
echo "-+-+-+- Starting matchstick.sh -+-+-+-"
# Requests for user input to create matchsticks
#read -p "Please enter the lowest voltage matchstick to be created (in V) - minimum is 1V: " lowerLimit
#read -p "Please enter the highest voltage matchstick to be created (in V) - maximum is 10V: " upperLimit
#read -p "Please enter the step size between matchsticks: " stepSize # make sure you have an integer number of steps between the lower and upper limits!
## matchsticks start at 0.2V pulser voltage up to 1.0V in 0.1V steps, and then 1V to 6V in 1V steps. This is hardcoded below. If you want this to change then edit the script below.
read -p "Please enter ten times the statistics required in each matchstick peak: " stats
read -p "Please enter the rate you want the pulser to run at: " rate
# calculate pause time required to accumulate required statistics
((pauseTime = stats/rate))
#echo $pauseTime
#read -p "Please enter the attenuation factor you want the pulser to use (either 1, 2, 5, 10, 20, 50 or 100): " atten
## default attenuation factor is 20X - if this needs to be changed, change it below.
read -p "Please enter which polarity pulse you want the pulser to output (0 = positive, 1 = negative): " polarityChoice
echo "Other relevant pulser parameters, and their default values, are:"
echo "Delay = 250 nanoseconds"
echo "Rise Time = 50 nanoseconds"
echo "Fall Time = 10 microseconds"
echo "If you want to change these default parameters please edit them in the 'matchsticks.sh' script." # change them below, when they are passed into /dev/USB1
echo
# Start sending commands to the pulser to set it up
echo "Setting up pulser..."
printf "set trigger mode internal\r" > /dev/ttyUSB1
sleep 1
printf "set tail pulse 1\r" > /dev/ttyUSB1
sleep 1
printf "clamp baseline 0\r" > /dev/ttyUSB1
sleep 1
printf "set delay 250\r" > /dev/ttyUSB1
sleep 1
printf "set rise time 100\r" > /dev/ttyUSB1
sleep 1
printf "set fall time 100000\r" > /dev/ttyUSB1
sleep 1
printf "set attenuation 20\r" > /dev/ttyUSB1
sleep 1
printf "set rep rate $rate\r" > /dev/ttyUSB1
sleep 1
printf "set polarity positive $polarityChoice\r" > /dev/ttyUSB1
sleep 1
echo "Done!"
# Send commands to pulser to make the matchsticks
echo "Making matchsticks..."
# 0.03V matchstick
voltage=$(echo 'scale=1;3/100' | bc )
echo $voltage
printf "set amplitude $voltage\r" > /dev/ttyUSB1
sleep 1
printf "set pulse on 1\r" > /dev/ttyUSB1
sleep $pauseTime
printf "set pulse on 0\r" > /dev/ttyUSB1
sleep 1
# 0.07V matchstick
voltage=$(echo 'scale=1;7/100' | bc )
echo $voltage
printf "set amplitude $voltage\r" > /dev/ttyUSB1
sleep 1
printf "set pulse on 1\r" > /dev/ttyUSB1
sleep $pauseTime
printf "set pulse on 0\r" > /dev/ttyUSB1
sleep 1
# 0.1V matchstick
voltage=$(echo 'scale=1;10/100' | bc )
echo $voltage
printf "set amplitude $voltage\r" > /dev/ttyUSB1
sleep 1
printf "set pulse on 1\r" > /dev/ttyUSB1
sleep $pauseTime
printf "set pulse on 0\r" > /dev/ttyUSB1
sleep 1
# 0.2V matchstick
voltage=$(echo 'scale=1;20/100' | bc )
echo $voltage
printf "set amplitude $voltage\r" > /dev/ttyUSB1
sleep 1
printf "set pulse on 1\r" > /dev/ttyUSB1
sleep $pauseTime
printf "set pulse on 0\r" > /dev/ttyUSB1
sleep 1
# 0.3V matchstick
voltage=$(echo 'scale=1;30/100' | bc )
echo $voltage
printf "set amplitude $voltage\r" > /dev/ttyUSB1
sleep 1
printf "set pulse on 1\r" > /dev/ttyUSB1
sleep $pauseTime
printf "set pulse on 0\r" > /dev/ttyUSB1
sleep 1
# 0.4V matchstick
voltage=$(echo 'scale=1;40/100' | bc )
echo $voltage
printf "set amplitude $voltage\r" > /dev/ttyUSB1
sleep 1
printf "set pulse on 1\r" > /dev/ttyUSB1
sleep $pauseTime
printf "set pulse on 0\r" > /dev/ttyUSB1
sleep 1
# 0.5V matchstick
voltage=$(echo 'scale=1;50/100' | bc )
echo $voltage
printf "set amplitude $voltage\r" > /dev/ttyUSB1
sleep 1
printf "set pulse on 1\r" > /dev/ttyUSB1
sleep $pauseTime
printf "set pulse on 0\r" > /dev/ttyUSB1
sleep 1
# 0.6V matchstick
voltage=$(echo 'scale=1;60/100' | bc )
echo $voltage
printf "set amplitude $voltage\r" > /dev/ttyUSB1
sleep 1
printf "set pulse on 1\r" > /dev/ttyUSB1
sleep $pauseTime
printf "set pulse on 0\r" > /dev/ttyUSB1
sleep 1
# 0.7V matchstick
voltage=$(echo 'scale=1;70/100' | bc )
echo $voltage
printf "set amplitude $voltage\r" > /dev/ttyUSB1
sleep 1
printf "set pulse on 1\r" > /dev/ttyUSB1
sleep $pauseTime
printf "set pulse on 0\r" > /dev/ttyUSB1
sleep 1
# 0.8V matchstick
voltage=$(echo 'scale=1;80/100' | bc )
echo $voltage
printf "set amplitude $voltage\r" > /dev/ttyUSB1
sleep 1
printf "set pulse on 1\r" > /dev/ttyUSB1
sleep $pauseTime
printf "set pulse on 0\r" > /dev/ttyUSB1
sleep 1
# 0.9V matchstick
voltage=$(echo 'scale=1;90/100' | bc )
echo $voltage
printf "set amplitude $voltage\r" > /dev/ttyUSB1
sleep 1
printf "set pulse on 1\r" > /dev/ttyUSB1
sleep $pauseTime
printf "set pulse on 0\r" > /dev/ttyUSB1
sleep 1
# 1.0V matchstick
voltage=$(echo 'scale=1;100/100' | bc )
echo $voltage
printf "set amplitude $voltage\r" > /dev/ttyUSB1
sleep 1
printf "set pulse on 1\r" > /dev/ttyUSB1
# doubling the statistics in this peak so it can be more easily identified in matchsticks spectra
sleep $pauseTime
sleep $pauseTime
printf "set pulse on 0\r" > /dev/ttyUSB1
sleep 1
# 1.5V matchstick
voltage=$(echo 'scale=1;150/100' | bc )
echo $voltage
printf "set amplitude $voltage\r" > /dev/ttyUSB1
sleep 1
printf "set pulse on 1\r" > /dev/ttyUSB1
sleep $pauseTime
printf "set pulse on 0\r" > /dev/ttyUSB1
sleep 1
# 2.0V matchstick
voltage=$(echo 'scale=1;200/100' | bc )
echo $voltage
printf "set amplitude $voltage\r" > /dev/ttyUSB1
sleep 1
printf "set pulse on 1\r" > /dev/ttyUSB1
sleep $pauseTime
printf "set pulse on 0\r" > /dev/ttyUSB1
sleep 1
# 2.5V matchstick
voltage=$(echo 'scale=1;250/100' | bc )
echo $voltage
printf "set amplitude $voltage\r" > /dev/ttyUSB1
sleep 1
printf "set pulse on 1\r" > /dev/ttyUSB1
sleep $pauseTime
printf "set pulse on 0\r" > /dev/ttyUSB1
sleep 1
# 3.0V matchstick
voltage=$(echo 'scale=1;300/100' | bc )
echo $voltage
printf "set amplitude $voltage\r" > /dev/ttyUSB1
sleep 1
printf "set pulse on 1\r" > /dev/ttyUSB1
sleep $pauseTime
printf "set pulse on 0\r" > /dev/ttyUSB1
sleep 1
# 3.5V matchstick
voltage=$(echo 'scale=1;350/100' | bc )
echo $voltage
printf "set amplitude $voltage\r" > /dev/ttyUSB1
sleep 1
printf "set pulse on 1\r" > /dev/ttyUSB1
sleep $pauseTime
printf "set pulse on 0\r" > /dev/ttyUSB1
sleep 1
# 4.0V matchstick
voltage=$(echo 'scale=1;400/100' | bc )
echo $voltage
printf "set amplitude $voltage\r" > /dev/ttyUSB1
sleep 1
printf "set pulse on 1\r" > /dev/ttyUSB1
sleep $pauseTime
printf "set pulse on 0\r" > /dev/ttyUSB1
sleep 1
# 4.5V matchstick
voltage=$(echo 'scale=1;450/100' | bc )
echo $voltage
printf "set amplitude $voltage\r" > /dev/ttyUSB1
sleep 1
printf "set pulse on 1\r" > /dev/ttyUSB1
sleep $pauseTime
printf "set pulse on 0\r" > /dev/ttyUSB1
sleep 1
# 5.0V matchstick
voltage=$(echo 'scale=1;500/100' | bc )
echo $voltage
printf "set amplitude $voltage\r" > /dev/ttyUSB1
sleep 1
printf "set pulse on 1\r" > /dev/ttyUSB1
sleep $pauseTime
printf "set pulse on 0\r" > /dev/ttyUSB1
sleep 1
# 5.5V matchstick
voltage=$(echo 'scale=1;550/100' | bc )
echo $voltage
printf "set amplitude $voltage\r" > /dev/ttyUSB1
sleep 1
printf "set pulse on 1\r" > /dev/ttyUSB1
sleep $pauseTime
printf "set pulse on 0\r" > /dev/ttyUSB1
sleep 1
# 6.0V matchstick
voltage=$(echo 'scale=1;600/100' | bc )
echo $voltage
printf "set amplitude $voltage\r" > /dev/ttyUSB1
sleep 1
printf "set pulse on 1\r" > /dev/ttyUSB1
sleep $pauseTime
printf "set pulse on 0\r" > /dev/ttyUSB1
sleep 1
echo "Done!"
exit

cc=gcc
${cc} inter.c -lm -lpthread -fopenmp -o main
./main $1 $2 $3

if [ $1 == 'cubic' ]
then
	python3 drawcubic.py
else
	python3 draw.py $1
fi
rm main

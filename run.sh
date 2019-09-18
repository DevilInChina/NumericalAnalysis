para='gau doo cro cho jac gs sor cg'
File='3matrices'
gcc main.c -o main -lm -lpthread
mat=$(ls ${File})
for name in ${para}
do
	for matrix in ${mat}
	do
		./main $name ${File}/$matrix
	done
done

CC = icpc  -g
CFLAGS =# -O3 

C = $(CC) $(CFLAGS)
EXEC = BSwhguan

.PHONY : clean

$(EXEC):
	$(C) -c $(EXEC).cpp
	$(C) $(EXEC).o -o $(EXEC)
BS:
	./$(EXEC)  T1whguan.txt 1
clean:
	rm -f *.o
	rm -f $(EXEC)
	echo "clean obj"

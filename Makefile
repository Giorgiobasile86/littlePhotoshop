main_iplib: main_iplib.o bmp.o ip_lib.o
	gcc bmp.o main_iplib.o ip_lib.o -Wall -lm -g -O1 -fsanitize=address -fsanitize=undefined -Wextra -o main_iplib
bmp.o: bmp.c bmp.h ip_lib.h
	gcc -c -Wall -std=gnu89 bmp.c -o bmp.o
ip_lib.o: ip_lib.c ip_lib.h bmp.h
	gcc -c -Wall --ansi --pedantic -std=gnu89 ip_lib.c -o ip_lib.o
main_iplib.o: main_iplib.c bmp.h ip_lib.h
	gcc -Wall --ansi --pedantic -std=gnu89 -c main_iplib.c -o main_iplib.o

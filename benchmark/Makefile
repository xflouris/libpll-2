#    Copyright (C) 2015 Diego Darriba
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as
#    published by the Free Software Foundation, either version 3 of the
#    License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#    Contact: Diego Darriba <Diego.Darriba@h-its.org>,
#    Exelixis Lab, Heidelberg Instutute for Theoretical Studies
#    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany

LIB_DIR = ../src/.libs
ABS_LIB_DIR=`readlink -m "../src/.libs"`

CC = gcc
CFLAGS = -I../test/src -L$(LIB_DIR) -g -O3 -Wall -std=c99
CLIBS = -lpll -lm

CFILES = $(shell find src -name '*.c')

OBJFILES = $(patsubst src/benchmarks/%.c, obj/%, $(CFILES))

OBJCOMMON = ../test/src/common.o

all: $(OBJCOMMON) $(OBJFILES)

obj/test-cache: src/benchmarks/test-cache.c obj/test_avx2.o obj/test_avx512.o ../test/src/common.c $(DEPS)
	@mkdir -p "$(@D)"
	$(CC) "-Wl,-rpath,$(ABS_LIB_DIR)" $(CPPFLAGS) $(CFLAGS) -o $@ $< obj/test_avx2.o obj/test_avx512.o $(OBJCOMMON) $(CLIBS) $(LDFLAGS)

obj/%: src/benchmarks/%.c ../test/src/common.c $(DEPS)
	@mkdir -p "$(@D)"
	$(CC) "-Wl,-rpath,$(ABS_LIB_DIR)" $(CPPFLAGS) $(CFLAGS) -o $@ $< $(OBJCOMMON) $(CLIBS) $(LDFLAGS)

obj/test_avx2.o: src/lib/test_avx2.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -mavx2 -mfma -c -o $@ $< $(CLIBS) $(LDFLAGS)

obj/test_avx512.o: src/lib/test_avx512.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -mavx512f -c -o $@ $< $(CLIBS) $(LDFLAGS)

%.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c -o $@ $< $(CLIBS) $(LDFLAGS)
	
clean:
	rm -rf obj result/* $(OBJCOMMON)

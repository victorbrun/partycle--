COMPILER 	= clang++
FLAGS 		= -Wall -o2
SOURCES		= main.cpp superellipsoid.cpp coordinate_indexer.cpp
HEADERS		= superellipsoid.h coordinate_indexer.h
OUTPUT 		= partycle--

compile: ${SOURCES} ${HEADERS}
	${COMPILER} ${FLAGS} ${SOURCES} -o ${OUTPUT}

run: ${OUTPUT}
	./${OUTPUT}

compile_run: compile run



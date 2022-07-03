CC = g++ -std=c++11
CFLAGS = -I/.
#LIB_PATH = vsc-hard-mounts/leuven-apps/thinking/thinking/2014a/software/GSL/1.16-foss-2014a/lib
LIB_PATH = vsc-hard-mounts/leuven-apps/thinking/2018a/software/GSL/2.4-GCCcore-6.4.0/lib
LIB = -lgsl -lgslcblas

DEPS1 = error.h dictionary.h meteorology.h openField.h UopenField.h KopenField.h math.h forest.h Uforest.h Kforest.h KTaylorOF.h KTaylorForest.h #windPowerLaw.h eddyDiffTaylor.h
DEPS2 = error.h dictionary.h meteorology.h openField.h UopenField.h KopenField.h math.h forest.h Uforest.h Kforest.h KTaylorOF.h KTaylorForest.h particle.h radiation.h #windPowerLaw.h eddyDiffTaylor.h
DEPS3 = error.h dictionary.h concentration.h radiation.h
DEPS4 = error.h dictionary.h meteorology.h particle.h concentration.h openField.h UopenField.h KopenField.h math.h forest.h Uforest.h Kforest.h KTaylorOF.h KTaylorForest.h
DEPS6 = error.h dictionary.h radiation.h

ODIR = include

OBJ_ = genMeteo.o error.o dictionary.o meteorology.o openField.o UopenField.o KopenField.o math.o forest.o Uforest.o Kforest.o KTaylorOF.o KTaylorForest.o# windPowerLaw.o eddyDiffTaylor.o
OBJ = $(patsubst %, $(ODIR)/%, $(OBJ_))

OBJ2_ = paths.o error.o dictionary.o meteorology.o openField.o UopenField.o KopenField.o math.o forest.o Uforest.o Kforest.o KTaylorOF.o KTaylorForest.o particle.o radiation.o
OBJ2 = $(patsubst %, $(ODIR)/%, $(OBJ2_))

OBJ3_ = reconstruct.o error.o dictionary.o concentration.o radiation.o
OBJ3 = $(patsubst %, $(ODIR)/%, $(OBJ3_))

OBJ4_ = simConc.o error.o dictionary.o meteorology.o particle.o concentration.o openField.o UopenField.o KopenField.o math.o forest.o Uforest.o Kforest.o KTaylorOF.o KTaylorForest.o
OBJ4 = $(patsubst %, $(ODIR)/%, $(OBJ4_))

OBJ6_ = posToDose.o error.o dictionary.o radiation.o
OBJ6 = $(patsubst %, $(ODIR)/%, $(OBJ6_))

$(ODIR)/%.o: %.cpp $(DEPS2)
	$(CC) -c -o $@ $< $(CFLAGS)

default: conc
 
EXE = genMeteo
EXE2 = simDose
EXE3 = reconstruct
EXE4 = simConc
EXE5 = posToDose

meteo: $(OBJ)
	$(CC) -L/$(LIB_PATH) -o $(EXE) $^ $(LIB)

dose: $(OBJ2)
	$(CC) -L/$(LIB_PATH) -o $(EXE2) $^ $(LIB)

reconstruct: $(OBJ3)
	$(CC) -o $(EXE3) $^

conc: $(OBJ4)
	$(CC) -L/$(LIB_PATH) -o $(EXE4) $^ $(LIB)

posdose: $(OBJ5)
	$(CC) -o $(EXE5) $^

clean: 
	rm include/*.o; rm output/*; rm $(EXE); rm $(EXE2); rm $(EXE3); rm $(EXE4);
